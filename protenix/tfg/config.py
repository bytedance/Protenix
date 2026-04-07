# Copyright 2024 ByteDance and/or its affiliates.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Configuration utilities for TFG (Training-Free Guidance).

This module converts a user-facing config mapping (e.g. loaded from JSON/YAML)
into strongly-typed runtime objects used by the TFG sampler.

What this file defines
----------------------
- `Schedule`: small callable objects that map a scalar time `t` in `[0, 1]` to a
  float (used for weights and other time-varying hyper-parameters).
- `Term`: a wrapper around a potential implementation from
  `protenix.tfg.potentials`, with an evaluation interval, a weight schedule, and
  (optionally scheduled) parameters.
- `TFGConfig`: the top-level configuration consumed by `TFGEngine`.

"""

import math
from dataclasses import dataclass
from typing import Any, Iterable, Mapping

import torch

from protenix.tfg import potentials


class Schedule:
    """
    Lightweight schedule interface.

    A schedule maps a scalar time `t` (typically normalized to `[0, 1]`) to a
    Python `float`.
    """

    def __call__(self, t: float) -> float:
        """Evaluate the schedule at time `t`.

        Args:
            t: Normalized time, typically in `[0, 1]`.

        Returns:
            A Python `float`.
        """
        raise NotImplementedError


@dataclass(frozen=True)
class Constant(Schedule):
    """A schedule that always returns the same value."""

    value: float

    def __call__(self, t: float) -> float:
        """Return `self.value` (ignores `t`)."""
        return float(self.value)


@dataclass(frozen=True)
class ExponentialInterpolation(Schedule):
    """
    Exponential interpolation between `start` and `end`.

    `alpha == 0` reduces to linear interpolation.
    """

    start: float
    end: float
    alpha: float = 0.0

    def __call__(self, t: float) -> float:
        """Interpolate between `start` and `end`.

        The interpolation is linear when `alpha == 0`, otherwise it applies an
        exponential warping so the change can be concentrated near `t=0` or
        `t=1`.

        Args:
            t: Normalized time in `[0, 1]`.

        Returns:
            Interpolated scalar.
        """
        if self.alpha == 0.0:
            return float(self.start + (self.end - self.start) * t)
        num = math.exp(self.alpha * t) - 1.0
        den = math.exp(self.alpha) - 1.0
        return float(self.start + (self.end - self.start) * (num / den))


def schedule_from_cfg(obj: Any) -> Schedule:
    """Parse a schedule from a config object.

    Supported forms:
    - `Schedule`: returned as-is.
    - `float`/`int`: treated as a constant schedule.
    - `Mapping`: schedule specification with a `type` key.

      Currently supported schedule types:
      - `{"type": "const", "value": 0.1}`
      - `{"type": "exp_interpolation", "start": 1.0, "end": 0.0, "alpha": 3.0}`
    """

    if isinstance(obj, Schedule):
        return obj
    if isinstance(obj, (int, float)):
        return Constant(float(obj))
    if isinstance(obj, Mapping):
        cfg = dict(obj)
        if "type" not in cfg:
            raise KeyError(
                "Schedule config must contain key 'type'. "
                "Examples: {'type': 'const', 'value': 1.0}"
            )
        t = str(cfg["type"]).lower()
        if t == "const":
            return Constant(float(cfg["value"]))
        if t == "exp_interpolation":
            return ExponentialInterpolation(
                start=float(cfg["start"]),
                end=float(cfg["end"]),
                alpha=float(cfg.get("alpha", 0.0)),
            )
        raise ValueError(f"Unknown schedule type: {t}")

    raise TypeError(f"Unsupported schedule config type: {type(obj)}")


_REQUIRED_FEATURES: dict[str, set[str]] = {
    # Potential class name -> required feature keys.
    #
    # These keys must exist in the `input_feature_dict` passed into
    # `TFGEngine.step(...)`. The list is best-effort and is intended for
    # fail-fast validation rather than exhaustive schema enforcement.
    "PairwiseDistancePotential": {
        "pairwise_distance_index",
        "pairwise_distance_is_bond",
        "pairwise_distance_is_angle",
        "pairwise_distance_upper_bound",
        "pairwise_distance_lower_bound",
        "ref_element",
    },
    "InterchainBondPotential": {"interchain_bond_index"},
    "VinaStericPotential": {
        "asym_id",
        "atom_to_token_idx",
        "ref_element",
        "interchain_bond_index",
    },
    "SymmetricChainPotential": {
        "asym_id",
        "atom_to_token_idx",
        "symmetric_chain_index",
    },
    "StereoBondPotential": {"stereo_bond_index", "stereo_bond_orientation"},
    "ChiralAtomPotential": {
        "chiral_index",
        "chiral_orientation",
        "asym_id",
        "atom_to_token_idx",
    },
    "PlanarImproperPotential": {
        "planar_improper_index",
        "planar_improper_is_carbonyl",
    },
    "LinearBondPotential": {"linear_triple_bond_index"},
    "ExperimentalTorsionPotential": {
        "experimental_torsion_index",
        "experimental_torsion_force_constant",
        "experimental_torsion_sign",
    },
}


@dataclass
class Term:
    """A single TFG energy term.

    A term wraps:
    - a potential implementation (from `protenix.tfg.potentials`),
    - an evaluation `interval` in inner-loop steps,
    - a time-varying `weight` schedule (multiplies energy and gradient),
    - optional (possibly scheduled) potential parameters.
    """

    # Name of the potential class, and the lookup key in `potentials.CLASS_REGISTRY`.
    name: str
    # Evaluate this term every `interval` diffusion steps. `interval <= 0` disables.
    interval: int
    # Multiplicative weight schedule applied to both energy and gradient.
    weight: Schedule
    # Potential parameters; each value may be a constant or a `Schedule`.
    param_templates: dict[str, Any]
    # Instantiated potential object implementing `energy(...)` and optionally
    # `energy_and_grad(...)` / `project(...)`.
    _potential: Any
    # Whether this term is allowed to participate in the projection loop.
    enable_projection: bool = True

    def required_features(self) -> set[str]:
        """Return required feature keys for this term (best-effort)."""
        return set(_REQUIRED_FEATURES.get(self.name, set()))

    def active(self, step_i: int) -> bool:
        """Whether this term should be evaluated at `step_i`."""
        return self.interval > 0 and (step_i % self.interval == 0)

    def _params_at(self, t: float) -> dict[str, Any]:
        """Materialize parameters for time `t`.

        Any parameter value that is a `Schedule` is evaluated at `t`; other
        values are passed through unchanged.
        """
        out = {}
        for k, v in self.param_templates.items():
            out[k] = v(t) if isinstance(v, Schedule) else v
        return out

    def energy(
        self, coords: torch.Tensor, feats: Mapping[str, Any], t: float
    ) -> torch.Tensor:
        """Compute weighted energy for a single term.

        Args:
            coords: Coordinates, shape `[*batch, N_atom, 3]`.
            feats: Feature dict for this sample.
            t: Normalized time in `[0, 1]`.

        Returns:
            Energy tensor with shape `[*batch]`.
        """
        w = float(self.weight(t))
        if w == 0.0:
            return torch.zeros(
                coords.shape[:-2], device=coords.device, dtype=coords.dtype
            )
        e = self._potential.energy(coords, feats, self._params_at(t))
        return e * w

    def energy_and_grad(
        self, coords: torch.Tensor, feats: Mapping[str, Any], t: float
    ) -> tuple[torch.Tensor, torch.Tensor]:
        """Compute weighted energy and gradient for a single term.

        Args:
            coords: Coordinates, shape `[*batch, N_atom, 3]`.
            feats: Feature dict for this sample.
            t: Normalized time in `[0, 1]`.

        Returns:
            (energy, grad)
            - energy: shape `[*batch]`
            - grad: shape `[*batch, N_atom, 3]`
        """
        w = float(self.weight(t))
        if w == 0.0:
            zero_e = torch.zeros(
                coords.shape[:-2], device=coords.device, dtype=coords.dtype
            )
            return zero_e, torch.zeros_like(coords)
        e, g = self._potential.energy_and_grad(coords, feats, self._params_at(t))
        return e * w, g * w

    def project(
        self, coords: torch.Tensor, feats: Mapping[str, Any], t: float
    ) -> torch.Tensor:
        """Optional projection update.

        Some potentials implement an auxiliary `project()` method (e.g. to enforce
        constraints). If not available, returns zeros.
        """
        if not hasattr(self._potential, "project"):
            # For potentials without a projection operator, we return zeros so
            # the caller can safely accumulate deltas.
            return torch.zeros_like(coords)

        return self._potential.project(coords, feats, self._params_at(t))


def _build_terms(term_cfg: Mapping[str, Any] | None) -> list[Term]:
    """Build `Term` objects from a config section.

    `term_cfg` must be a mapping:
    `{term_name: {interval, weight, ...potential_params}}`.
    """
    terms = []
    if term_cfg is None:
        return terms

    if not isinstance(term_cfg, Mapping):
        raise TypeError(
            f"terms must be a mapping of term_name -> term_config, got {type(term_cfg)}"
        )

    items = term_cfg.items()

    for name, raw in items:
        if raw is None:
            raw = {}
        raw = dict(raw)
        interval = int(raw.pop("interval", 1))
        weight_cfg = raw.pop("weight", 0.0)
        enable_projection = bool(raw.pop("enable_projection", True))
        weight = schedule_from_cfg(weight_cfg)
        param_templates = {}
        for k, v in raw.items():
            if isinstance(v, dict) and "type" in v:
                # Schedule spec.
                try:
                    param_templates[k] = schedule_from_cfg(v)
                    continue
                except Exception:
                    # Not a valid schedule spec; keep raw dict.
                    pass
            param_templates[k] = v

        if name not in potentials.CLASS_REGISTRY:
            raise KeyError(
                f"Unknown potential '{name}'. Available: {sorted(potentials.CLASS_REGISTRY.keys())}"
            )
        pot_cls = potentials.CLASS_REGISTRY[name]
        # Instantiate the potential. Parameters are passed at call-time so the
        # same instance can be reused across diffusion steps.
        pot = pot_cls()
        terms.append(
            Term(
                name=name,
                interval=interval,
                weight=weight,
                param_templates=param_templates,
                _potential=pot,
                enable_projection=enable_projection,
            )
        )
    return terms


def validate_features(feats: Mapping[str, Any], terms: Iterable[Term]) -> None:
    """
    Validate that `feats` contains required keys for all configured terms.
    """
    missing: dict[str, list[str]] = {}
    for term in terms:
        need = term.required_features()
        if not need:
            continue
        miss = [k for k in need if k not in feats]
        if miss:
            missing[term.name] = miss
    if missing:
        lines = ["TFG is missing required input features:"]
        for name, miss in missing.items():
            lines.append(f"- {name}: {miss}")
        raise KeyError("\n".join(lines))


@dataclass(frozen=True)
class TFGConfig:
    """Top-level configuration consumed by `TFGEngine`.

    Fields are intentionally plain Python types so they are easy to construct
    from JSON/YAML, serialize, and log.
    """

    # Whether TFG is enabled at all.
    enable: bool
    # Strength of denoiser-path guidance (gradient through the denoiser).
    rho: float
    # Step size for direct refinement on predicted x0.
    mu: float
    # MC noise std used by the logp estimator.
    eps_std: float
    # Number of MC samples.
    eps_batch: int
    # Outer loop count inside one diffusion step.
    outer_steps: int
    # Inner refinement loop count on x0 inside one outer step.
    inner_steps: int
    # Number of projection iterations.
    projection_outer_steps: int
    projection_inner_steps: int
    # Configured terms (potentials) to apply.
    terms: tuple[Term, ...]
    # Debug switch: log per-term energies at the last refinement step.
    log_last_step_energy: bool = False


def parse_tfg_config(guidance_cfg: Mapping[str, Any] | None) -> TFGConfig:
    """
    Parse a TFG config mapping into a `TFGConfig`.

    Expected top-level keys (all optional):
        - `enable`: bool
        - `rho`: float
        - `mu`: float
        - `mc`: mapping with `std` and `batch`
        - `steps`: mapping with `tfg_outer`, `tfg_inner`, `projection`
        - `terms`: mapping `{PotentialName: {interval, weight, ...params}}`
        - `log_last_step_energy`: bool

    This parser performs basic validation and normalization:
        - clamps batch sizes / step counts to be non-negative where applicable,
        - ensures `eps_batch >= 1`,
        - errors out on unknown top-level keys to avoid silent misconfiguration.
    """
    if guidance_cfg is None:
        # No guidance config means guidance is completely disabled.
        return TFGConfig(
            enable=False,
            rho=0.0,
            mu=0.0,
            eps_std=0.0,
            eps_batch=1,
            outer_steps=1,
            inner_steps=0,
            projection_outer_steps=0,
            projection_inner_steps=0,
            terms=(),
        )

    cfg = dict(guidance_cfg)
    allowed_top = {
        "enable",
        "rho",
        "mu",
        "mc",
        "steps",
        "terms",
        "log_last_step_energy",
    }
    extra = set(cfg.keys()) - allowed_top
    if extra:
        raise KeyError(f"Unsupported keys in TFG config: {sorted(extra)}")

    enable = bool(cfg.get("enable", False))
    rho = float(cfg.get("rho", 0.0))
    mu = float(cfg.get("mu", 0.0))

    mc = dict(cfg.get("mc", {}))
    eps_std = float(mc.get("std", 0.0))
    eps_batch = int(mc.get("batch", 1))
    eps_batch = max(1, eps_batch)

    steps = dict(cfg.get("steps", {}))

    outer_steps = max(1, int(steps.get("tfg_outer", 1)))
    inner_steps = max(0, int(steps.get("tfg_inner", 10)))
    projection_outer_steps = max(0, int(steps.get("projection_outer", 2)))
    projection_inner_steps = max(0, int(steps.get("projection_inner", 10)))
    term_cfg = cfg.get("terms", {})
    terms = tuple(_build_terms(term_cfg))
    if enable and len(terms) == 0:
        raise ValueError("TFG is enabled but no terms are configured")

    return TFGConfig(
        enable=enable,
        rho=rho,
        mu=mu,
        eps_std=eps_std,
        eps_batch=eps_batch,
        outer_steps=outer_steps,
        inner_steps=inner_steps,
        projection_inner_steps=projection_inner_steps,
        projection_outer_steps=projection_outer_steps,
        terms=terms,
        log_last_step_energy=bool(cfg.get("log_last_step_energy", False)),
    )
