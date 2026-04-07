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

"""TFG (Training-Free Guidance) energy terms (Potentials).

This module defines differentiable geometry/chemistry constraints implemented as
*potentials* (penalty energies) operating on atomic coordinates.

Core conventions:

- `coords`: atomic coordinates, typically shaped `[..., N, 3]`.
- `feats`: a dict containing indices/bounds/labels needed by constraints.
- `energy(...)`: returns a per-sample scalar energy.
- `energy_and_grad(...)`: returns `(energy, dE/dcoords)`.
- `project(...)`: returns a coordinate update `delta_x` (same shape as `coords`).

Internally we compute values and Jacobians for basic geometric primitives
(pairwise distance, angle, dihedral), then map `dE/d(value)` back to
`dE/d(coords)` via the chain rule.
"""

import math
from typing import Any, Optional

import torch

from protenix.data.constants import rdkit_vdws

CLASS_REGISTRY: dict[str, type] = {}

# Cache VDW radii tables to avoid per-call allocation.
# Note: keep the canonical table in fp32 for numerical stability; cast when used.
_VDW_RADII_128_CPU = torch.zeros(128, dtype=torch.float32)
_VDW_RADII_128_CPU[:118] = torch.as_tensor(rdkit_vdws, dtype=torch.float32)
_VDW_RADII_128_CACHE: dict[tuple[str, int | None], torch.Tensor] = {}


def _get_vdw_radii_128(device: torch.device) -> torch.Tensor:
    """Return a cached `[128]` fp32 VDW radii table on `device`."""
    key = (device.type, device.index)
    out = _VDW_RADII_128_CACHE.get(key)
    if out is None or out.device != device:
        out = _VDW_RADII_128_CPU.to(device=device)
        _VDW_RADII_128_CACHE[key] = out
    return out


def register(cls):
    """Decorator: register a potential class by its name.

    This enables constructing potentials dynamically from configuration.
    """
    CLASS_REGISTRY[cls.__name__] = cls
    return cls


class Potential:
    """Base class of all potentials.

    Subclasses implement `_eval(...)` to compute energy, and may implement
    `_project(...)` to provide a projection-like coordinate correction.
    """

    def __init__(
        self,
        default_params: Optional[dict[str, Any]] = None,
    ):
        self._default_params = default_params

    def _resolve_params(self, params: Optional[dict[str, Any]]) -> dict[str, Any]:
        """Merge default parameters with runtime parameters."""
        res = {}
        if self._default_params is not None:
            res.update(self._default_params)
        if params is not None:
            res.update(params)
        return res

    def _eval(
        self,
        coords: torch.Tensor,
        feats: dict[str, Any],
        params: dict[str, Any],
        need_grad: bool,
    ) -> torch.Tensor | tuple[torch.Tensor, torch.Tensor]:
        """Compute energy and optionally gradients.

        Args:
            coords: Atomic coordinates, `[..., N, 3]`.
            feats: Feature dict.
            params: Resolved params (defaults + overrides).
            need_grad: Whether to return coordinate gradients.

        Returns:
            If `need_grad=False`: energy tensor of shape `coords.shape[:-2]` (or scalar).
            If `need_grad=True`: `(energy, grad)` where `grad` has the same shape as `coords`.
        """
        raise NotImplementedError

    def energy(
        self,
        coords: torch.Tensor,
        feats: dict[str, Any],
        params: Optional[dict[str, Any]] = None,
    ) -> torch.Tensor:
        """Compute energy only (no gradients)."""
        p = self._resolve_params(params)
        return self._eval(coords, feats, p, False)

    def energy_and_grad(
        self,
        coords: torch.Tensor,
        feats: dict[str, Any],
        params: Optional[dict[str, Any]] = None,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        """Compute energy and coordinate gradients."""
        p = self._resolve_params(params)
        return self._eval(coords, feats, p, True)

    def project(
        self,
        coords: torch.Tensor,
        feats: dict[str, Any],
        params: Optional[dict[str, Any]] = None,
    ) -> torch.Tensor:
        """Compute a projection-like coordinate update `delta_x`."""
        p = self._resolve_params(params)
        return self._project(coords, feats, p)

    def _project(
        self, coords: torch.Tensor, feats: dict[str, Any], params: dict[str, Any]
    ) -> torch.Tensor:
        """Subclass hook for projection. Default: no-op (zeros)."""
        return torch.zeros_like(coords)


def _to_tensor(x, device, dtype):
    """Convert a Python scalar/array to a tensor on the given device/dtype."""
    return torch.as_tensor(x, device=device, dtype=dtype)


def _distance_value_and_grad(
    coords: torch.Tensor, index: torch.Tensor, need_grad: bool
) -> tuple[torch.Tensor, Optional[torch.Tensor]]:
    """Pairwise distances and (optional) Jacobians w.r.t. coordinates.

    Args:
        coords: `[..., N, 3]`.
        index: Atom-pair indices `[2, M]` (each column is (i, j)).
        need_grad: Whether to compute `d(dist)/d(coords)`.

    Returns:
        value: Distances `[..., M]`.
        grad_value: If requested, Jacobian `[..., 2, M, 3]` where the `2` axis
            corresponds to atoms (i, j). Otherwise None.
    """
    r = coords.index_select(-2, index[0]) - coords.index_select(-2, index[1])
    norm = torch.linalg.norm(r, dim=-1).clamp_min(1e-8)
    if not need_grad:
        return norm, None
    r_hat = r / norm.unsqueeze(-1)
    grad = torch.stack((r_hat, -r_hat), dim=-3)
    return norm, grad


def _angle_value_and_grad(
    coords: torch.Tensor, index: torch.Tensor, need_grad: bool
) -> tuple[torch.Tensor, Optional[torch.Tensor]]:
    """Angles for triples (i-j-k) and (optional) Jacobians.

    The angle is defined at vertex `j` between vectors `(i - j)` and `(k - j)`.
    The output is in radians (typically in [0, pi]).

    Shapes:
        - `index`: `[3, M]`
        - return angle: `[..., M]`
        - return grad (if requested): `[..., 3, M, 3]` for atoms (i, j, k)
    """
    r_ji = coords.index_select(-2, index[0]) - coords.index_select(-2, index[1])
    r_jk = coords.index_select(-2, index[2]) - coords.index_select(-2, index[1])
    cross = torch.linalg.cross(r_ji, r_jk, dim=-1)
    dot = torch.linalg.vecdot(r_ji, r_jk) + 1e-8
    angle = torch.atan2(torch.linalg.vector_norm(cross, dim=-1), dot)

    if not need_grad:
        return angle, None

    r2_ji = torch.linalg.vecdot(r_ji, r_ji) + 1e-8
    r2_jk = torch.linalg.vecdot(r_jk, r_jk) + 1e-8
    rp = torch.linalg.vector_norm(cross, dim=-1) + 1e-8

    grad_i = torch.linalg.cross(r_ji, cross, dim=-1) / (r2_ji * rp).unsqueeze(-1)
    grad_k = torch.linalg.cross(cross, r_jk, dim=-1) / (r2_jk * rp).unsqueeze(-1)
    grad_j = -grad_i - grad_k
    grad = torch.stack((grad_i, grad_j, grad_k), dim=-3)
    return angle, grad


def _dihedral_value_and_grad(
    coords: torch.Tensor, index: torch.Tensor, need_grad: bool
) -> tuple[torch.Tensor, Optional[torch.Tensor]]:
    """Dihedral (torsion) angles and (optional) Jacobians.

    For quadruples (i-j-k-l), the dihedral is the signed angle between planes
    (i-j-k) and (j-k-l), in radians.

    Shapes:
        - `index`: `[4, M]`
        - return phi: `[..., M]`
        - return grad (if requested): `[..., 4, M, 3]` for atoms (i, j, k, l)
    """
    r_ij = coords.index_select(-2, index[1]) - coords.index_select(-2, index[0])
    r_kj = coords.index_select(-2, index[1]) - coords.index_select(-2, index[2])
    r_kl = coords.index_select(-2, index[3]) - coords.index_select(-2, index[2])

    m = torch.linalg.cross(r_ij, r_kj, dim=-1)
    n = torch.linalg.cross(r_kj, r_kl, dim=-1)
    w = torch.linalg.cross(m, n, dim=-1)
    wlen = torch.linalg.vector_norm(w, dim=-1)
    s = torch.linalg.vecdot(m, n) + 1e-8
    phi = torch.atan2(wlen, s)
    ipr = torch.linalg.vecdot(r_ij, n)
    phi = -phi * torch.sign(ipr)

    if not need_grad:
        return phi, None

    iprm = torch.linalg.vecdot(m, m) + 1e-8
    iprn = torch.linalg.vecdot(n, n) + 1e-8
    nrkj2 = torch.linalg.vecdot(r_kj, r_kj)
    nrkj_1 = torch.rsqrt(nrkj2 + 1e-8)
    nrkj_2 = torch.square(nrkj_1)
    nrkj = nrkj2 * nrkj_1
    a = -nrkj / iprm
    f_i = -a.unsqueeze(-1) * m
    b = nrkj / iprn
    f_l = -b.unsqueeze(-1) * n
    p = torch.linalg.vecdot(r_ij, r_kj) * nrkj_2
    q = torch.linalg.vecdot(r_kl, r_kj) * nrkj_2
    uvec = p.unsqueeze(-1) * f_i
    vvec = q.unsqueeze(-1) * f_l
    svec = uvec - vvec
    f_j = f_i - svec
    f_k = f_l + svec
    grad = torch.stack((f_i, -f_j, -f_k, f_l), dim=-3)
    return phi, grad


def _abs_dihedral_value_and_grad(
    coords: torch.Tensor, index: torch.Tensor, need_grad: bool
) -> tuple[torch.Tensor, Optional[torch.Tensor]]:
    """Absolute dihedral value |phi| and (optional) Jacobians.

    Useful when only the magnitude matters (planarity, cis/trans) and the sign does not.
    """
    phi, grad = _dihedral_value_and_grad(coords, index, need_grad)
    if not need_grad:
        return torch.abs(phi), None
    sign = torch.where(phi < 0, -1.0, 1.0)
    grad = grad * sign[..., None, :, None]
    return torch.abs(phi), grad


def _planar_improper_value_and_grad(
    coords: torch.Tensor,
    index: torch.Tensor,
    need_grad: bool,
    zero_tol: float = 1e-10,
) -> tuple[torch.Tensor, Optional[torch.Tensor]]:
    # index: [4, M]
    # index[0]: nb1, index[1]: nb2, index[2]: center, index[3]: nb3
    p1 = coords.index_select(-2, index[0])  # iPoint
    p3 = coords.index_select(-2, index[1])  # kPoint
    p2 = coords.index_select(-2, index[2])  # jPoint (Center)
    p4 = coords.index_select(-2, index[3])  # lPoint
    r_ji = p1 - p2
    r_jk = p3 - p2
    r_jl = p4 - p2
    l2_ji = torch.linalg.vecdot(r_ji, r_ji)
    l2_jk = torch.linalg.vecdot(r_jk, r_jk)
    l2_jl = torch.linalg.vecdot(r_jl, r_jl)
    dist_mask = (l2_ji > zero_tol) & (l2_jk > zero_tol) & (l2_jl > zero_tol)
    d_ji = torch.sqrt(l2_ji).clamp_min(1e-8)
    d_jk = torch.sqrt(l2_jk).clamp_min(1e-8)
    d_jl = torch.sqrt(l2_jl).clamp_min(1e-8)

    r_ji_n = r_ji / d_ji.unsqueeze(-1)
    r_jk_n = r_jk / d_jk.unsqueeze(-1)
    r_jl_n = r_jl / d_jl.unsqueeze(-1)
    # n = -rJI x rJK
    n = torch.linalg.cross(-r_ji_n, r_jk_n, dim=-1)
    l2_n = torch.linalg.vecdot(n, n)
    mask = dist_mask & (l2_n > zero_tol)
    n_normalized = n / torch.sqrt(l2_n).unsqueeze(-1).clamp_min(1e-8)
    cos_y = torch.where(mask, torch.linalg.vecdot(n_normalized, r_jl_n), 0.0).clamp(
        -1.0, 1.0
    )
    sin_y_sq = (1.0 - cos_y * cos_y).clamp_min(0.0)
    sin_y = torch.sqrt(sin_y_sq).clamp_min(1e-8)
    energy = 1.0 - sin_y
    if not need_grad:
        return energy, None
    cos_theta = torch.linalg.vecdot(r_ji_n, r_jk_n).clamp(-1.0, 1.0)
    sin_theta_sq = (1.0 - cos_theta * cos_theta).clamp_min(1e-12)
    sin_theta = torch.sqrt(sin_theta_sq).clamp_min(1e-8)
    t1 = torch.linalg.cross(r_jl_n, r_jk_n, dim=-1)
    t2 = torch.linalg.cross(r_ji_n, r_jl_n, dim=-1)
    t3 = torch.linalg.cross(r_jk_n, r_ji_n, dim=-1)
    term1 = sin_y * sin_theta
    term2 = cos_y / (sin_y * sin_theta_sq)

    tg1 = (
        t1 / term1.unsqueeze(-1)
        - (r_ji_n - r_jk_n * cos_theta.unsqueeze(-1)) * term2.unsqueeze(-1)
    ) / d_ji.unsqueeze(-1)
    tg3 = (
        t2 / term1.unsqueeze(-1)
        - (r_jk_n - r_ji_n * cos_theta.unsqueeze(-1)) * term2.unsqueeze(-1)
    ) / d_jk.unsqueeze(-1)
    tg4 = (
        t3 / term1.unsqueeze(-1) - r_jl_n * (cos_y / sin_y).unsqueeze(-1)
    ) / d_jl.unsqueeze(-1)
    tg2 = -(tg1 + tg3 + tg4)
    # dE/dY = -cos_y (where sin_y = sin(Y), cos_y = cos(Y))
    # grad_E = (dE/dY) * grad_Y = (dE/dY) * (-tg) = -(-cos_y) * tg = cos_y * tg
    # cos_y shape: [..., M]
    # grad shape: [..., 4, M, 3]
    grad = torch.stack((tg1, tg3, tg2, tg4), dim=-3) * cos_y.unsqueeze(-2).unsqueeze(-1)
    grad = torch.where(mask.unsqueeze(-2).unsqueeze(-1), grad, 0.0)
    return energy, grad


def _flat_bottom_linear(
    value: torch.Tensor,
    k: torch.Tensor,
    lower: Optional[torch.Tensor],
    upper: Optional[torch.Tensor],
) -> tuple[torch.Tensor, torch.Tensor]:
    """
    Flat-bottom linear potential.

    When the value is outside the interval [lower, upper], the energy grows
    linearly with the violation:

        E = k * (max(0, lower - value) + max(0, value - upper))

    Returns:
        energy: Energy tensor with the same shape as `value`.
        grad: Gradient dE/dvalue with the same shape as `value`.
    """
    energy = torch.zeros_like(value)
    grad = torch.zeros_like(value)

    if lower is not None:
        diff_lb = lower - value
        mask_lb = diff_lb > 0
        energy += torch.where(mask_lb, k * diff_lb, 0.0)
        grad -= torch.where(mask_lb, k, 0.0)

    if upper is not None:
        diff_ub = value - upper
        mask_ub = diff_ub > 0
        energy += torch.where(mask_ub, k * diff_ub, 0.0)
        grad += torch.where(mask_ub, k, 0.0)

    return energy, grad


def _flat_bottom_parabolic(
    value: torch.Tensor,
    k: torch.Tensor,
    lower: Optional[torch.Tensor],
    upper: Optional[torch.Tensor],
) -> tuple[torch.Tensor, torch.Tensor]:
    """
    Flat-bottom parabolic potential.
    When the value is outside the interval [lower, upper], the energy grows
    quadratically with the violation:
        E = 0.5 * k * (max(0, lower - value)^2 + max(0, value - upper)^2)
    Returns:
        energy: Energy tensor with the same shape as `value`.
        grad: Gradient dE/dvalue with the same shape as `value`.
    """
    energy = torch.zeros_like(value)
    grad = torch.zeros_like(value)
    if lower is not None:
        diff_lb = lower - value
        mask_lb = diff_lb > 0
        energy += torch.where(mask_lb, 0.5 * k * (diff_lb**2), 0.0)
        grad -= torch.where(mask_lb, k * diff_lb, 0.0)
    if upper is not None:
        diff_ub = value - upper
        mask_ub = diff_ub > 0
        energy += torch.where(mask_ub, 0.5 * k * (diff_ub**2), 0.0)
        grad += torch.where(mask_ub, k * diff_ub, 0.0)
    return energy, grad


def _aggregate_atom_gradients(
    coords: torch.Tensor,
    indices: torch.Tensor,
    grad_values: torch.Tensor,
    dE_dvals: torch.Tensor,
) -> torch.Tensor:
    """Map constraint-space gradients to per-atom coordinate gradients.

    Given `dE/dv` and `dv/dx`, compose `dE/dx = (dE/dv) * (dv/dx)` and accumulate
    contributions back to a dense `[..., N, 3]` tensor via `scatter_add_`.
    """
    # Compose gradients: dE/dx = dE/dv * dv/dx.
    # Broadcast: [..., 1, M, 1] * [..., K, M, 3] -> [..., K, M, 3]
    g_combined = dE_dvals.unsqueeze(-2).unsqueeze(-1) * grad_values

    g_flat = g_combined.flatten(start_dim=-3, end_dim=-2)  # [..., K*M, 3]
    idx_flat = indices.flatten()  # [K*M]

    if g_flat.ndim > coords.ndim:
        extra_dims = g_flat.ndim - coords.ndim
        g_flat = g_flat.sum(dim=list(range(extra_dims)))

    batch_shape = coords.shape[:-2]
    target_idx = idx_flat.view(*([1] * len(batch_shape)), -1, 1).expand(
        *batch_shape, -1, 3
    )

    # Accumulate into per-atom coordinate gradients.
    out = torch.zeros_like(coords)
    out.scatter_add_(dim=-2, index=target_idx, src=g_flat)
    return out


def _solve_constraint_projection(
    coords: torch.Tensor,
    index: torch.Tensor,
    value: torch.Tensor,
    grad_value: torch.Tensor,
    mask: torch.Tensor,
    eps: float = 1e-8,
) -> torch.Tensor:
    """Solve a masked linearized constraint projection and return `delta_x`.

    For active constraints we use the minimum-norm solution:
        dx = -J^T (J J^T)^{-1} v
    where `v` is the constraint violation and `J = dv/dx` is the Jacobian.
    """
    # Shapes:
    #   coords:      [*B, N, 3]
    #   index:       [K, M]
    #   value:       [*B, M]
    #   grad_value:  [*B, K, M, 3]  (dv/dx)
    #   mask:        [*B, M] (bool)
    device, dtype = coords.device, coords.dtype
    batch_shape = coords.shape[:-2]
    n_atom = coords.shape[-2]
    k = int(index.shape[0])
    m_total = int(index.shape[1])

    if value.shape[-1] != m_total:
        raise ValueError(
            f"_solve_constraint_projection: value.shape[-1]={value.shape[-1]} != index.shape[1]={m_total}"
        )
    if mask.shape != value.shape:
        raise ValueError(
            f"_solve_constraint_projection: mask.shape={tuple(mask.shape)} != value.shape={tuple(value.shape)}"
        )
    if grad_value.shape[-3:] != (k, m_total, 3):
        raise ValueError(
            f"_solve_constraint_projection: grad_value.shape[-3:]={tuple(grad_value.shape[-3:])} "
            f"!= (K,M,3)=({k},{m_total},3)"
        )

    # Flatten batch dims -> loop per batch (each batch has a different active mask).
    b = math.prod(batch_shape) if len(batch_shape) > 0 else 1
    coords_b = coords.reshape(b, n_atom, 3)
    value_b = value.reshape(b, m_total)
    mask_b = mask.reshape(b, m_total)
    grad_b = grad_value.reshape(b, k, m_total, 3)

    delta = torch.zeros_like(coords_b)
    atom_arange_cache: dict[int, torch.Tensor] = {}

    # Linearized projection: dx = -J^T (J J^T)^{-1} v.
    for bi in range(b):
        active = mask_b[bi]
        n_active = int(active.sum().item())
        if n_active == 0:
            continue

        v = value_b[bi, active]  # [m]
        g = grad_b[bi, :, active, :]  # [K, m, 3]
        idx = index[:, active]  # [K, m]

        # Build Jacobian rows in atom space: [m, N, 3] -> [m, N*3].
        grad_atom = torch.zeros((n_active, n_atom, 3), device=device, dtype=dtype)
        if n_active not in atom_arange_cache:
            atom_arange_cache[n_active] = torch.arange(n_active, device=device)
        row = atom_arange_cache[n_active]
        for ki in range(k):
            grad_atom[row, idx[ki]] += g[ki]
        grad_flat = grad_atom.reshape(n_active, n_atom * 3)

        # Solve (J J^T + eps I) lambda = v in fp32 for numerical stability if needed.
        work_dtype = (
            torch.float32 if dtype in (torch.float16, torch.bfloat16) else dtype
        )
        gf = grad_flat.to(work_dtype)
        vv = v.to(work_dtype)
        denom = gf @ gf.transpose(-1, -2)
        denom = denom + (
            float(eps) * torch.eye(n_active, device=device, dtype=work_dtype)
        )
        lam = torch.linalg.solve(denom, vv)  # [m]
        dx_flat = -(gf.transpose(-1, -2) @ lam)  # [N*3]
        delta[bi] = dx_flat.to(dtype).reshape(n_atom, 3)

    return delta.reshape(*coords.shape)


def _zeros_energy(coords: torch.Tensor) -> torch.Tensor:
    """Return a zero energy tensor matching the batch shape."""
    return torch.zeros(coords.shape[:-2], device=coords.device, dtype=coords.dtype)


def _zeros_energy_and_grad(coords: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
    """Return `(0 energy, 0 gradients)` with shapes matching `coords`."""
    return _zeros_energy(coords), torch.zeros_like(coords)


def _sum_energy(e: torch.Tensor) -> torch.Tensor:
    """Sum per-constraint energies along the last axis into a per-sample scalar."""
    if e.ndim == 0:
        return e
    return e.sum(dim=-1)


# -----------------------------------------------------------------------------
# Concrete terms
# -----------------------------------------------------------------------------


@register
class InterchainBondPotential(Potential):
    """Upper-bound distance constraint for inter-chain covalent bonds.

    Expected `feats`:
        - `interchain_bond_index`: `[2, M]` atom-pair indices.

    Param:
        - `buffer`: maximum allowed distance for these bonds (same unit as `coords`).
    """

    def __init__(self, default_params: Optional[dict[str, Any]] = None):
        defaults = {"buffer": 2.0}
        if default_params is not None:
            defaults.update(default_params)
        super().__init__(defaults)

    def _eval(self, coords, feats, params, need_grad: bool):
        """Compute energy (and optionally gradients) for inter-chain bonds."""
        idx = feats["interchain_bond_index"]
        if idx.numel() == 0:
            return (
                _zeros_energy_and_grad(coords) if need_grad else _zeros_energy(coords)
            )

        value, grad_value = _distance_value_and_grad(coords, idx, need_grad)

        k = torch.ones_like(value)
        upper = torch.full(
            (idx.shape[1],),
            float(params["buffer"]),
            device=coords.device,
            dtype=value.dtype,
        )
        e, dE = _flat_bottom_linear(value, k, None, upper)
        if not need_grad:
            return _sum_energy(e)
        grad_atom = _aggregate_atom_gradients(coords, idx, grad_value, dE)
        return _sum_energy(e), grad_atom


@register
class PairwiseDistancePotential(Potential):
    """Pairwise distance bounds with category-specific buffers and VDW corrections.

    This term enforces lower/upper bounds for a list of atom pairs. Different
    buffers are used depending on whether a pair is a bond, angle-derived, or a
    generic non-bonded (clash) constraint.

    Expected `feats`:
        - `pairwise_distance_index`: `[2, M]`
        - `pairwise_distance_lower_bound`: `[M]`
        - `pairwise_distance_upper_bound`: `[M]`
        - `pairwise_distance_is_bond`: `[M]` (0/1)
        - `pairwise_distance_is_angle`: `[M]` (0/1)
        - `ref_element`: `[N, 118]` one-hot element type used for VDW radii

    Params:
        - `bond_buffer`, `angle_buffer`, `clash_buffer`: relative slack factors.
    """

    def __init__(self, default_params: Optional[dict[str, Any]] = None):
        defaults = {"bond_buffer": 0.05, "angle_buffer": 0.05, "clash_buffer": 0.05}
        if default_params is not None:
            defaults.update(default_params)
        super().__init__(defaults)

        # Cache static per-sample tensors (depends only on feats).
        self._cache_key = None
        self._cached_state_code = None  # int64 [M]
        self._cached_vdw_limit_f32 = None  # float32 [M]

        # Cache small scale tables keyed by (device, dtype, buffers).
        self._scale_table_cache: dict[tuple, tuple[torch.Tensor, torch.Tensor]] = {}

    @staticmethod
    def _cache_key_from_feats(feats: dict[str, Any]) -> tuple:
        idx = feats["pairwise_distance_index"]
        lb = feats["pairwise_distance_lower_bound"]
        ub = feats["pairwise_distance_upper_bound"]
        is_bond = feats["pairwise_distance_is_bond"]
        is_angle = feats["pairwise_distance_is_angle"]
        ref_element = feats["ref_element"]

        def _tid(x: torch.Tensor) -> tuple:
            return (x.data_ptr(), tuple(x.shape), str(x.dtype))

        return (
            str(idx.device),
            _tid(idx),
            _tid(lb),
            _tid(ub),
            _tid(is_bond),
            _tid(is_angle),
            _tid(ref_element),
        )

    def _get_scale_tables(
        self,
        *,
        device: torch.device,
        lb_dtype: torch.dtype,
        ub_dtype: torch.dtype,
        b_buf: float,
        a_buf: float,
        c_buf: float,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        min_ba = min(b_buf, a_buf)
        key = (
            device.type,
            device.index,
            str(lb_dtype),
            str(ub_dtype),
            b_buf,
            a_buf,
            c_buf,
        )
        out = self._scale_table_cache.get(key)
        if out is not None:
            return out

        l_table = torch.tensor(
            [c_buf, b_buf, a_buf, min_ba], device=device, dtype=lb_dtype
        )
        u_table = torch.tensor(
            [0.0, b_buf, a_buf, min_ba], device=device, dtype=ub_dtype
        )
        self._scale_table_cache[key] = (l_table, u_table)
        return l_table, u_table

    def _get_distance_bounds(self, feats, params):
        """Compute final (lower, upper) bounds per pair.

        Notes:
            - Clash constraints get `upper = +inf` (i.e., only a lower bound).
            - VDW radii are used to clamp bounds to avoid unrealistic contacts.
        """
        idx = feats["pairwise_distance_index"]
        if idx.numel() == 0:
            return idx, None, None

        lb_base = feats["pairwise_distance_lower_bound"]
        ub_base = feats["pairwise_distance_upper_bound"]

        cache_key = self._cache_key_from_feats(feats)
        if cache_key != self._cache_key:
            # Rebuild static caches.
            is_bond = feats["pairwise_distance_is_bond"].to(torch.int64)
            is_angle = feats["pairwise_distance_is_angle"].to(torch.int64)

            # Combined state code: 0=clash, 1=bond, 2=angle, 3=both (bond & angle).
            self._cached_state_code = is_bond + (is_angle << 1)

            # Precompute VDW-based per-pair limit in fp32 (independent of buffers).
            vdw_radii = _get_vdw_radii_128(idx.device)
            element_idx = feats["ref_element"].argmax(dim=-1)
            atom_radii = vdw_radii.index_select(
                0, element_idx.to(torch.long)
            )  # [N] fp32
            vdw_pair_sum = atom_radii[idx[0]] + atom_radii[idx[1]]
            self._cached_vdw_limit_f32 = 0.35 + 0.5 * vdw_pair_sum

            self._cache_key = cache_key

        state_code = self._cached_state_code
        vdw_limit = self._cached_vdw_limit_f32.to(dtype=lb_base.dtype)

        b_buf = float(params["bond_buffer"])
        a_buf = float(params["angle_buffer"])
        c_buf = float(params["clash_buffer"])

        l_table, u_table = self._get_scale_tables(
            device=idx.device,
            lb_dtype=lb_base.dtype,
            ub_dtype=ub_base.dtype,
            b_buf=b_buf,
            a_buf=a_buf,
            c_buf=c_buf,
        )

        l_scale = 1.0 - l_table[state_code]
        u_scale = 1.0 + u_table[state_code]

        lower = lb_base * l_scale
        upper = ub_base * u_scale
        # Clash: no upper bound.
        upper = torch.where(state_code == 0, float("inf"), upper)

        # Van der Waals radii (VDW) adjustment.
        bond_mask = (state_code & 1).bool()
        lower = torch.where(~bond_mask, torch.maximum(lower, vdw_limit), lower)
        upper = torch.where(bond_mask, torch.minimum(upper, vdw_limit), upper)

        return idx, lower, upper

    def _eval(self, coords, feats, params, need_grad: bool):
        """Compute energies (and optionally gradients) for all constrained pairs."""
        idx, lower, upper = self._get_distance_bounds(feats, params)
        if idx.numel() == 0:
            return (
                _zeros_energy_and_grad(coords) if need_grad else _zeros_energy(coords)
            )

        value, grad_value = _distance_value_and_grad(coords, idx, need_grad)

        k = torch.ones_like(value)
        e, dE = _flat_bottom_parabolic(value, k, lower, upper)
        if not need_grad:
            return _sum_energy(e)
        grad_atom = _aggregate_atom_gradients(coords, idx, grad_value, dE)
        return _sum_energy(e), grad_atom

    def _project_masked(self, coords, idx, lower, upper, idx_mask):
        idx = idx[..., idx_mask]
        if idx.numel() == 0:
            return torch.zeros_like(coords)
        lower = lower[idx_mask]
        upper = upper[idx_mask]
        value, grad_value = _distance_value_and_grad(coords, idx, True)
        mask_lb = value < lower
        mask_ub = value > upper
        mask = mask_lb | mask_ub
        if mask.sum() == 0:
            return torch.zeros_like(coords)
        v = value - torch.where(mask_lb, lower, upper)
        return _solve_constraint_projection(coords, idx, v, grad_value, mask)

    def _project(self, coords, feats, params):
        """Project on bond-distance and angle-distance pairs back into their [lower, upper] interval."""
        idx = feats["pairwise_distance_index"]
        idx, lower, upper = self._get_distance_bounds(feats, params)
        bond_mask = feats["pairwise_distance_is_bond"].bool()
        angle_mask = feats["pairwise_distance_is_angle"].bool()
        delta_x_angle = self._project_masked(
            coords,
            idx,
            lower,
            upper,
            angle_mask,
        )
        delta_x_bond = self._project_masked(
            coords + delta_x_angle, idx, lower, upper, bond_mask
        )
        return delta_x_bond + delta_x_angle


@register
class StereoBondPotential(Potential):
    """Stereo double-bond (cis/trans) constraint using |dihedral|.

    Expected `feats`:
        - `stereo_bond_index`: `[4, M]` torsion quadruples.
        - `stereo_bond_orientation`: `[M]` float/bool-like.
          If `> 0.5`, the preferred state is trans (in terms of the selected torsion quadruple), |phi| near pi.
          Otherwise, the preferred state is cis (in terms of the selected torsion quadruple), |phi| near 0.

    Param:
        - `buffer`: tolerance in radians.
    """

    def __init__(self, default_params: Optional[dict[str, Any]] = None):
        defaults = {"buffer": 0.52360}
        if default_params is not None:
            defaults.update(default_params)
        super().__init__(defaults)

    def _eval(self, coords, feats, params, need_grad: bool):
        """Compute the stereo-bond penalty energy (and optionally gradients)."""
        idx = feats["stereo_bond_index"]
        orient = feats["stereo_bond_orientation"].to(coords.dtype)
        if idx.numel() == 0:
            return (
                _zeros_energy_and_grad(coords) if need_grad else _zeros_energy(coords)
            )

        value, grad_value = _abs_dihedral_value_and_grad(coords, idx, need_grad)
        lower = torch.where(
            orient > 0.5,
            _to_tensor(torch.pi - params["buffer"], value.device, value.dtype),
            _to_tensor(float("-inf"), value.device, value.dtype),
        )
        upper = torch.where(
            orient > 0.5,
            _to_tensor(float("inf"), value.device, value.dtype),
            _to_tensor(params["buffer"], value.device, value.dtype),
        )

        k = torch.ones_like(value)
        e, dE = _flat_bottom_linear(value, k, lower, upper)
        e_sum = _sum_energy(e)
        if not need_grad:
            return e_sum
        grad_atom = _aggregate_atom_gradients(coords, idx, grad_value, dE)
        return e_sum, grad_atom


@register
class ChiralAtomPotential(Potential):
    """Chirality constraint for chiral centers.

    The chirality is represented as a torsion (dihedral) angle `phi` with an
    orientation label `orient` (typically +1 or -1). We enforce:

        orient > 0:  phi >=  buffer
        orient <= 0: phi <= -buffer

    The `_project` method provides a feasibility repair step using a linearized
    constraint projection.
    """

    def __init__(self, default_params: Optional[dict[str, Any]] = None):
        defaults = {"buffer": 0.34906, "scale_x": True}
        if default_params is not None:
            defaults.update(default_params)
        super().__init__(defaults)

    def _eval(self, coords, feats, params, need_grad: bool):
        """Compute chirality penalty energy (and optionally gradients)."""
        idx = feats["chiral_index"]
        if idx.numel() == 0:
            return (
                _zeros_energy_and_grad(coords) if need_grad else _zeros_energy(coords)
            )

        orient = feats["chiral_orientation"]
        value, grad_value = _dihedral_value_and_grad(coords, idx, need_grad)

        lower = torch.where(
            orient > 0,
            _to_tensor(params["buffer"], value.device, value.dtype),
            _to_tensor(float("-inf"), value.device, value.dtype),
        )
        upper = torch.where(
            orient > 0,
            _to_tensor(float("inf"), value.device, value.dtype),
            -_to_tensor(params["buffer"], value.device, value.dtype),
        )
        k = torch.ones_like(value)
        e, dE = _flat_bottom_linear(value, k, lower, upper)
        e_sum = _sum_energy(e)
        if not need_grad:
            return e_sum
        grad_atom = _aggregate_atom_gradients(coords, idx, grad_value, dE)
        return e_sum, grad_atom

    def _project(self, coords, feats, params):
        """Project violated chirality constraints using a linearized solver.

        Returns:
            `delta_x` with the same shape as `coords`.
        """
        idx = feats["chiral_index"]
        if idx.numel() == 0:
            return torch.zeros_like(coords)
        orient = feats["chiral_orientation"]
        buffer = float(params["buffer"])
        value, grad_value = _dihedral_value_and_grad(coords, idx, True)
        value = value * orient - buffer
        grad_value = grad_value * orient.reshape(
            *([1] * len(coords.shape[:-2])), 1, -1, 1
        )
        mask = value < 0
        if mask.sum() == 0:
            return torch.zeros_like(coords)
        delta_x = _solve_constraint_projection(coords, idx, value, grad_value, mask)

        if bool(params.get("scale_x", True)):
            # Apply a chain-wise rescaling so that the radius of gyration of
            # affected atoms remain unchanged
            device = coords.device
            dtype = coords.dtype
            batch_shape = coords.shape[:-2]
            n_atom = coords.shape[-2]
            mask_atom = torch.zeros(n_atom, dtype=bool, device=device)
            mask_atom[idx] = True
            atom_chain_id_full = feats["asym_id"][..., feats["atom_to_token_idx"]]
            atom_chain_id = atom_chain_id_full[mask_atom]
            n_chain = int(atom_chain_id.max()) + 1
            masked_delta_x = delta_x[..., mask_atom, :]
            masked_coords = coords[..., mask_atom, :]
            atom_chain_id_reshaped = atom_chain_id.reshape(
                *([1] * len(batch_shape)), -1, 1
            )
            atom_count = torch.zeros(
                *batch_shape, n_chain, 1, device=device, dtype=dtype
            ).scatter_add_(
                -2,
                atom_chain_id_reshaped.expand(*batch_shape, -1, 1),
                torch.ones_like(masked_coords[..., :1]),
            )
            coords_sum = torch.zeros(
                *batch_shape, n_chain, 3, device=device, dtype=dtype
            ).scatter_add_(
                -2, atom_chain_id_reshaped.expand(*batch_shape, -1, 3), masked_coords
            )
            center = coords_sum / atom_count.clamp_min(1)
            center = center[..., atom_chain_id, :]
            rg_sum = torch.zeros(
                *batch_shape, n_chain, 1, device=device, dtype=dtype
            ).scatter_add_(
                -2,
                atom_chain_id_reshaped.expand(*batch_shape, -1, 1),
                ((masked_coords - center) ** 2).sum(dim=-1, keepdim=True),
            )
            rg = rg_sum / atom_count.clamp_min(1)
            rg = rg[..., atom_chain_id, :]

            new_x = masked_coords + masked_delta_x
            new_sum = torch.zeros(
                *batch_shape, n_chain, 3, device=device, dtype=dtype
            ).scatter_add_(
                -2, atom_chain_id_reshaped.expand(*batch_shape, -1, 3), new_x
            )
            new_center = new_sum / atom_count.clamp_min(1)
            new_center = new_center[..., atom_chain_id, :]
            new_rg_sum = torch.zeros(
                *batch_shape, n_chain, 1, device=device, dtype=dtype
            ).scatter_add_(
                -2,
                atom_chain_id_reshaped.expand(*batch_shape, -1, 1),
                ((new_x - new_center) ** 2).sum(dim=-1, keepdim=True),
            )
            new_rg = new_rg_sum / atom_count.clamp_min(1)
            new_rg = new_rg[..., atom_chain_id, :].clamp_min(1e-12)
            new_x = (new_x - new_center) * ((rg / new_rg) ** 0.5) + center
            masked_delta_x = new_x - masked_coords
            delta_x[..., mask_atom, :] = masked_delta_x

        return delta_x


@register
class PlanarImproperPotential(Potential):
    """Planarity constraint (typically for sp2 systems) using |dihedral|.

    Expected `feats`:
        - `planar_improper_index`: `[4, M]`
        - `planar_improper_is_carbonyl`: `[M]` float/bool-like.
          If `> 0.5`, prefer |phi| near pi; otherwise prefer |phi| near 0.

    Param:
        - `buffer`: tolerance in radians.
    """

    def __init__(self, default_params: Optional[dict[str, Any]] = None):
        defaults = {"buffer": 0.1309}
        if default_params is not None:
            defaults.update(default_params)
        super().__init__(defaults)

    def _eval(self, coords, feats, params, need_grad: bool):
        """Compute planarity penalty energy (and optionally gradients)."""
        idx = feats["planar_improper_index"]
        is_carbonyl = feats["planar_improper_is_carbonyl"].to(coords.dtype)
        if idx.numel() == 0:
            return (
                _zeros_energy_and_grad(coords) if need_grad else _zeros_energy(coords)
            )

        energy, grad = _planar_improper_value_and_grad(coords, idx, need_grad)
        k = torch.ones_like(is_carbonyl)
        energy = k * energy
        e_sum = _sum_energy(energy)

        if not need_grad:
            return e_sum
        grad_atom = _aggregate_atom_gradients(coords, idx, grad, k)
        return e_sum, grad_atom


@register
class LinearBondPotential(Potential):
    """Linearity constraint (typically for triple bonds) using angles.

    Expected `feats`:
        - `linear_triple_bond_index`: `[3, M]` angle triples.

    Param:
        - `buffer`: tolerance in radians; enforce angle >= pi - buffer.
    """

    def __init__(self, default_params: Optional[dict[str, Any]] = None):
        defaults = {"buffer": 0.08726646259}
        if default_params is not None:
            defaults.update(default_params)
        super().__init__(defaults)

    def _eval(self, coords, feats, params, need_grad: bool):
        """Compute linear-angle penalty energy (and optionally gradients)."""
        idx = feats["linear_triple_bond_index"]
        if idx.numel() == 0:
            return (
                _zeros_energy_and_grad(coords) if need_grad else _zeros_energy(coords)
            )

        value, grad_value = _angle_value_and_grad(coords, idx, need_grad)
        k = torch.ones_like(value)
        lower = torch.full(
            (idx.shape[1],),
            torch.pi - float(params["buffer"]),
            device=coords.device,
            dtype=value.dtype,
        )
        e, dE = _flat_bottom_linear(value, k, lower, None)
        e_sum = _sum_energy(e)
        if not need_grad:
            return e_sum
        grad_atom = _aggregate_atom_gradients(coords, idx, grad_value, dE)
        return e_sum, grad_atom


@register
class ExperimentalTorsionPotential(Potential):
    """Torsion energy from a cosine-series expansion (up to order 6).

    Expected `feats`:
        - `experimental_torsion_index`: `[4, M]`
        - `experimental_torsion_force_constant`: `[..., M, 6]`
        - `experimental_torsion_sign`: `[..., M, 6]`

    Energy:
        E(phi) = Σ_{n=1..6} k_n * (1 + s_n * cos(n*phi))
    """

    def _eval(self, coords, feats, params, need_grad: bool):
        """Compute torsion energy (and optionally gradients)."""
        idx = feats["experimental_torsion_index"]
        if idx.numel() == 0:
            return (
                _zeros_energy_and_grad(coords) if need_grad else _zeros_energy(coords)
            )

        force_constants = feats["experimental_torsion_force_constant"]
        signs = feats["experimental_torsion_sign"]

        phi, grad_phi = _dihedral_value_and_grad(coords, idx, need_grad)

        cosphi = torch.cos(phi)
        cosphi2 = cosphi * cosphi
        cosphi3 = cosphi * cosphi2
        cosphi4 = cosphi * cosphi3
        cosphi5 = cosphi * cosphi4
        cosphi6 = cosphi * cosphi5
        cos2 = 2.0 * cosphi2 - 1.0
        cos3 = 4.0 * cosphi3 - 3.0 * cosphi
        cos4 = 8.0 * cosphi4 - 8.0 * cosphi2 + 1.0
        cos5 = 16.0 * cosphi5 - 20.0 * cosphi3 + 5.0 * cosphi
        cos6 = 32.0 * cosphi6 - 48.0 * cosphi4 + 18.0 * cosphi2 - 1.0

        energy = (
            force_constants[..., 0] * (1.0 + signs[..., 0] * cosphi)
            + force_constants[..., 1] * (1.0 + signs[..., 1] * cos2)
            + force_constants[..., 2] * (1.0 + signs[..., 2] * cos3)
            + force_constants[..., 3] * (1.0 + signs[..., 3] * cos4)
            + force_constants[..., 4] * (1.0 + signs[..., 4] * cos5)
            + force_constants[..., 5] * (1.0 + signs[..., 5] * cos6)
        )
        e_sum = _sum_energy(energy)

        if not need_grad:
            return e_sum

        sinphi = torch.sin(phi)
        dcos2 = 4.0 * cosphi
        dcos3 = 12.0 * cosphi2 - 3.0
        dcos4 = 32.0 * cosphi3 - 16.0 * cosphi
        dcos5 = 80.0 * cosphi4 - 60.0 * cosphi2 + 5.0
        dcos6 = 192.0 * cosphi5 - 192.0 * cosphi3 + 36.0 * cosphi
        dE_dphi = -sinphi * (
            force_constants[..., 0] * signs[..., 0]
            + force_constants[..., 1] * signs[..., 1] * dcos2
            + force_constants[..., 2] * signs[..., 2] * dcos3
            + force_constants[..., 3] * signs[..., 3] * dcos4
            + force_constants[..., 4] * signs[..., 4] * dcos5
            + force_constants[..., 5] * signs[..., 5] * dcos6
        )

        grad_atom = _aggregate_atom_gradients(coords, idx, grad_phi, dE_dphi)
        return e_sum, grad_atom


@register
class VinaStericPotential(Potential):
    """Vina-style steric term for inter-chain (complex) contacts.

    This term builds candidate atom pairs across different chains (excluding
    chain pairs connected by inter-chain covalent bonds), assigns each pair an
    equilibrium distance based on VDW radii, and applies an AutoDock-Vina-like
    empirical potential (two Gaussians + quadratic repulsion).

    Param:
        - `buffer`: active-pair threshold factor. Only pairs with
          `dist < r_eq * (1 - buffer)` are evaluated to reduce O(N^2) cost.
    """

    def __init__(self, default_params: Optional[dict[str, Any]] = None):
        defaults = {"buffer": 0.225}
        if default_params is not None:
            defaults.update(default_params)
        super().__init__(defaults)

        # Cache collision candidates since they depend only on `feats` (static
        # per sample across diffusion steps) and are expensive to rebuild.
        self._cache_key = None
        self._cache_val = None

    @staticmethod
    def _cache_key_from_feats(feats: dict[str, Any]) -> tuple:
        """Build a stable cache key from feature tensor identities/shapes."""
        asym_id = feats["asym_id"]
        atom_to_token_idx = feats["atom_to_token_idx"]
        ref_element = feats["ref_element"]
        b_idx = feats.get("interchain_bond_index")

        def _tensor_id(x: torch.Tensor) -> tuple:
            return (x.data_ptr(), tuple(x.shape), str(x.dtype))

        return (
            str(asym_id.device),
            _tensor_id(asym_id),
            _tensor_id(atom_to_token_idx),
            _tensor_id(ref_element),
            None if b_idx is None else _tensor_id(b_idx),
        )

    def _get_collision_candidates(self, feats: dict[str, Any]):
        """Select atom pairs for the Vina-style steric term and equilibrium distances.

        Returns:
            sel_idx: `[2, M]` atom pair indices.
            r_eq: `[M]` equilibrium distance proxy (sum of VDW radii).
        """
        key = self._cache_key_from_feats(feats)
        if key == self._cache_key and self._cache_val is not None:
            return self._cache_val

        def _cache_and_return(out):
            self._cache_key = key
            self._cache_val = out
            return out

        device = feats["asym_id"].device
        c_map = feats["asym_id"][..., feats["atom_to_token_idx"]]
        n_atoms = c_map.size(-1)

        # This term is inter-chain only. If there is a single chain, early-exit
        # to avoid constructing O(N^2) candidate pairs.
        if n_atoms == 0:
            return _cache_and_return(
                torch.empty((2, 0), device=device, dtype=torch.long),
                torch.empty((0,), device=device, dtype=torch.float32),
            )
        # If all atoms map to the same chain id, there is no inter-chain pair.
        if torch.all(c_map == c_map[..., :1]):
            return _cache_and_return(
                torch.empty((2, 0), device=device, dtype=torch.long),
                torch.empty((0,), device=device, dtype=torch.float32),
            )

        n_chains = int(c_map.max().item()) + 1

        vdw_lib = _get_vdw_radii_128(device)
        r_atom = vdw_lib[feats["ref_element"].argmax(dim=-1)]

        prohibited = torch.eye(n_chains, dtype=torch.bool, device=device)

        b_idx = feats.get("interchain_bond_index")
        if b_idx is not None and b_idx.numel() > 0:
            ca, cb = c_map[b_idx[0]], c_map[b_idx[1]]
            prohibited[ca, cb] = True
            prohibited[cb, ca] = True

        # Build candidate pairs **only across chains** to avoid generating
        # intra-chain O(N^2) pairs and then masking them out.
        order = torch.argsort(c_map)  # [N]
        c_sorted = c_map.index_select(0, order)
        chain_ids, counts = torch.unique_consecutive(c_sorted, return_counts=True)
        if chain_ids.numel() <= 1:
            out = (
                torch.empty((2, 0), device=device, dtype=torch.long),
                torch.empty((0,), device=device, dtype=torch.float32),
            )
            return _cache_and_return(out)

        # Convert small metadata to CPU once; avoid device->host sync per pair.
        chain_ids_cpu = chain_ids.to("cpu").tolist()
        counts_cpu = counts.to("cpu").tolist()
        prohibited_cpu = prohibited.to("cpu")

        # Legacy behavior (previous `is_complex[i] & is_complex[j]`): only keep
        # pairs where BOTH atoms come from chains that have more than 1 atom.
        # This is often used to avoid over-constraining coordination-like cases
        # represented as singleton chains (e.g., metal ions).
        valid_chain = {
            int(cid) for cid, cnt in zip(chain_ids_cpu, counts_cpu) if int(cnt) > 1
        }
        if len(valid_chain) <= 1:
            out = (
                torch.empty((2, 0), device=device, dtype=torch.long),
                torch.empty((0,), device=device, dtype=torch.float32),
            )
            return _cache_and_return(out)

        # Build per-chain atom index lists from the sorted permutation.
        chain_atoms: dict[int, torch.Tensor] = {}
        start = 0
        for cid, cnt in zip(chain_ids_cpu, counts_cpu):
            end = start + int(cnt)
            if int(cid) in valid_chain:
                chain_atoms[int(cid)] = order[start:end]
            start = end

        i_list: list[torch.Tensor] = []
        j_list: list[torch.Tensor] = []
        for ii, ca in enumerate(chain_ids_cpu):
            if int(ca) not in valid_chain:
                continue
            a_atoms = chain_atoms[int(ca)]
            na = int(a_atoms.numel())
            if na == 0:
                continue
            for cb in chain_ids_cpu[ii + 1 :]:
                if int(cb) not in valid_chain:
                    continue
                if bool(prohibited_cpu[int(ca), int(cb)]):
                    continue
                b_atoms = chain_atoms[int(cb)]
                nb = int(b_atoms.numel())
                if nb == 0:
                    continue
                # Cartesian product of atom indices across chains.
                i_list.append(a_atoms.repeat_interleave(nb))
                j_list.append(b_atoms.repeat(na))

        if len(i_list) == 0:
            out = (
                torch.empty((2, 0), device=device, dtype=torch.long),
                torch.empty((0,), device=device, dtype=torch.float32),
            )
            return _cache_and_return(out)

        i_all = torch.cat(i_list, dim=0)
        j_all = torch.cat(j_list, dim=0)
        sel_idx = torch.stack((i_all, j_all), dim=0)
        r_eq = r_atom[sel_idx[0]] + r_atom[sel_idx[1]]

        out = (sel_idx, r_eq)
        return _cache_and_return(out)

    def _eval(self, coords, feats, params, need_grad: bool):
        """Compute Vina-style steric energy (and optionally gradients)."""
        idx, eq = self._get_collision_candidates(feats)
        if idx.numel() == 0:
            return (
                _zeros_energy_and_grad(coords) if need_grad else _zeros_energy(coords)
            )
        # Only compute energy/grad on active pairs:
        # - First compute distances for all candidate pairs to build an active mask
        #   (avoids allocating a huge dv/dx tensor).
        # - Then vectorize energy/grad over active pairs and scatter_add back.
        buf = float(params["buffer"])
        device, dtype = coords.device, coords.dtype

        batch_shape = coords.shape[:-2]
        n_atom = int(coords.shape[-2])
        b = math.prod(batch_shape) if len(batch_shape) > 0 else 1
        coords_b = coords.reshape(b, n_atom, 3)

        # Pass 1: distance only -> active mask
        value0, _ = _distance_value_and_grad(coords_b, idx, False)  # [B, M]
        thr = eq.to(value0.dtype) * (1.0 - buf)  # [M]
        active = value0 < thr.unsqueeze(0)  # [B, M]

        b_idx, m_idx = active.nonzero(as_tuple=True)  # [M_active]
        if b_idx.numel() == 0:
            return (
                _zeros_energy_and_grad(coords) if need_grad else _zeros_energy(coords)
            )

        # Gather active pair atom indices.
        i_atom = idx[0].index_select(0, m_idx)  # [M_active]
        j_atom = idx[1].index_select(0, m_idx)  # [M_active]
        eq_a = eq.index_select(0, m_idx).to(dtype)  # [M_active]

        # Energy on active pairs.
        v = value0[b_idx, m_idx].to(dtype)  # [M_active]
        dist_diff = v - eq_a
        norm_d = dist_diff / 0.5
        g1 = -0.0356 * torch.exp(-(norm_d**2))
        g2 = -0.00516 * torch.exp(-(((dist_diff - 3.0) / 2.0) ** 2))
        rep = 0.840 * torch.where(dist_diff < 0.0, dist_diff**2, 0.0)
        e_pair = g1 + g2 + rep  # [M_active]

        e_out = torch.zeros((b,), device=device, dtype=dtype)  # [B]
        e_out.scatter_add_(0, b_idx, e_pair)

        if not need_grad:
            if len(batch_shape) == 0:
                return e_out[0]
            return e_out.reshape(*batch_shape)

        # Pass 2: grad only on active pairs.
        ri = coords_b[b_idx, i_atom]  # [M_active, 3]
        rj = coords_b[b_idx, j_atom]
        r = ri - rj
        norm = torch.linalg.norm(r, dim=-1).clamp_min(1e-8)  # [M_active]
        r_hat = r / norm.unsqueeze(-1)

        # dE/d(dist)
        dg1 = -2.0 * g1 * norm_d * (1.0 / 0.5)
        dg2 = -0.5 * g2 * (dist_diff - 3.0)
        drep = 0.840 * torch.where(dist_diff < 0.0, 2.0 * dist_diff, 0.0)
        d_total = dg1 + dg2 + drep  # [M_active]

        # dE/dx = dE/d(dist) * d(dist)/dx
        gi = d_total.unsqueeze(-1) * r_hat  # [M_active, 3]
        gj = -gi

        g_flat = torch.zeros((b * n_atom, 3), device=device, dtype=dtype)  # [B*N, 3]
        flat_i = b_idx * n_atom + i_atom
        flat_j = b_idx * n_atom + j_atom
        g_flat.scatter_add_(0, flat_i.unsqueeze(-1).expand(-1, 3), gi)
        g_flat.scatter_add_(0, flat_j.unsqueeze(-1).expand(-1, 3), gj)
        g_out = g_flat.reshape(b, n_atom, 3).reshape(*coords.shape)  # [B, N, 3]

        if len(batch_shape) == 0:
            return e_out[0], g_out
        return e_out.reshape(*batch_shape), g_out
