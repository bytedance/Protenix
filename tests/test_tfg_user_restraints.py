import numpy as np
import pytest
import torch


def _user_distance_restraint_potential():
    pytest.importorskip("rdkit")
    from protenix.tfg.potentials import UserDistanceRestraintPotential

    return UserDistanceRestraintPotential()


def _feature_dict(index, lower, upper, weight=None):
    if weight is None:
        weight = [1.0] * len(lower)
    return {
        "user_distance_restraint_index": torch.tensor(index, dtype=torch.long),
        "user_distance_restraint_lower_bound": torch.tensor(
            lower, dtype=torch.float32
        ),
        "user_distance_restraint_upper_bound": torch.tensor(
            upper, dtype=torch.float32
        ),
        "user_distance_restraint_weight": torch.tensor(weight, dtype=torch.float32),
    }


def _empty_user_restraint_features():
    return {
        "user_distance_restraint_index": torch.empty((2, 0), dtype=torch.long),
        "user_distance_restraint_lower_bound": torch.empty((0,), dtype=torch.float32),
        "user_distance_restraint_upper_bound": torch.empty((0,), dtype=torch.float32),
        "user_distance_restraint_weight": torch.empty((0,), dtype=torch.float32),
    }


def test_user_distance_restraint_empty_features_return_zero():
    potential = _user_distance_restraint_potential()
    coords = torch.zeros((3, 3), dtype=torch.float32)
    feats = _empty_user_restraint_features()

    energy, grad = potential.energy_and_grad(coords, feats)

    torch.testing.assert_close(energy, torch.tensor(0.0))
    torch.testing.assert_close(grad, torch.zeros_like(coords))


def test_user_distance_restraint_upper_bound_violation_gradient():
    potential = _user_distance_restraint_potential()
    coords = torch.tensor([[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]])
    feats = _feature_dict([[0], [1]], lower=[0.0], upper=[3.0])

    energy, grad = potential.energy_and_grad(coords, feats)

    torch.testing.assert_close(energy, torch.tensor(0.5))
    assert grad[0, 0] < 0.0
    assert grad[1, 0] > 0.0


def test_user_distance_restraint_lower_bound_violation_gradient():
    potential = _user_distance_restraint_potential()
    coords = torch.tensor([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    feats = _feature_dict([[0], [1]], lower=[2.0], upper=[4.0])

    energy, grad = potential.energy_and_grad(coords, feats)

    torch.testing.assert_close(energy, torch.tensor(0.5))
    assert grad[0, 0] > 0.0
    assert grad[1, 0] < 0.0


def test_user_distance_restraint_zero_inside_bounds_and_weight_scaling():
    potential = _user_distance_restraint_potential()
    coords = torch.tensor([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]])
    inside = _feature_dict([[0], [1]], lower=[2.0], upper=[4.0], weight=[2.0])
    outside = _feature_dict([[0], [1]], lower=[0.0], upper=[2.0], weight=[2.0])

    torch.testing.assert_close(potential.energy(coords, inside), torch.tensor(0.0))
    torch.testing.assert_close(potential.energy(coords, outside), torch.tensor(1.0))


def test_user_distance_restraint_supports_batch_dimensions():
    potential = _user_distance_restraint_potential()
    coords = torch.tensor(
        [
            [[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]],
            [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
        ]
    )
    feats = _feature_dict([[0], [1]], lower=[0.0], upper=[3.0])

    energy = potential.energy(coords, feats)

    torch.testing.assert_close(energy, torch.tensor([0.5, 0.0]))


class DummyAtomArray:
    def __init__(self):
        self.shape = (4,)
        self.label_entity_id = np.array(["1", "1", "2", "2"])
        self.copy_id = np.array([1, 1, 1, 1])
        self.res_id = np.array([10, 10, 20, 20])
        self.atom_name = np.array(["CA", "CB", "N", "C"])


def _token_array():
    pytest.importorskip("biotite")
    from protenix.data.tokenizer import Token, TokenArray

    token1 = Token(1)
    token1.atom_indices = [0, 1]
    token1.atom_names = ["CA", "CB"]
    token1.centre_atom_index = 0
    token2 = Token(1)
    token2.atom_indices = [2, 3]
    token2.atom_names = ["N", "C"]
    token2.centre_atom_index = 2
    return TokenArray([token1, token2])


def test_tfg_contact_features_resolve_atom_contacts_exactly():
    pytest.importorskip("biotite")
    from protenix.data.constraint.constraint_featurizer import ConstraintFeatureGenerator

    constraint = {
        "contact": [
            {
                "entity1": 1,
                "copy1": 1,
                "position1": 10,
                "atom1": "CB",
                "entity2": 2,
                "copy2": 1,
                "position2": 20,
                "atom2": "C",
                "min_distance": 2.0,
                "max_distance": 6.0,
            }
        ]
    }

    feats = ConstraintFeatureGenerator.generate_tfg_contact_restraint_features(
        _token_array(), DummyAtomArray(), [{}, {}], constraint
    )

    torch.testing.assert_close(
        feats["user_distance_restraint_index"], torch.tensor([[1], [3]])
    )
    torch.testing.assert_close(
        feats["user_distance_restraint_lower_bound"], torch.tensor([2.0])
    )
    torch.testing.assert_close(
        feats["user_distance_restraint_upper_bound"], torch.tensor([6.0])
    )


def test_tfg_contact_features_resolve_token_contacts_to_centre_atoms():
    pytest.importorskip("biotite")
    from protenix.data.constraint.constraint_featurizer import ConstraintFeatureGenerator

    constraint = {
        "contact": [
            {
                "entity1": 1,
                "copy1": 1,
                "position1": 10,
                "entity2": 2,
                "copy2": 1,
                "position2": 20,
                "min_distance": 1.5,
                "max_distance": 7.0,
            }
        ]
    }

    feats = ConstraintFeatureGenerator.generate_tfg_contact_restraint_features(
        _token_array(), DummyAtomArray(), [{}, {}], constraint
    )

    torch.testing.assert_close(
        feats["user_distance_restraint_index"], torch.tensor([[0], [2]])
    )
    torch.testing.assert_close(
        feats["user_distance_restraint_lower_bound"], torch.tensor([1.5])
    )
    torch.testing.assert_close(
        feats["user_distance_restraint_upper_bound"], torch.tensor([7.0])
    )


def test_tfg_contact_features_preserve_mixed_atom_residue_contacts():
    pytest.importorskip("biotite")
    from protenix.data.constraint.constraint_featurizer import ConstraintFeatureGenerator

    constraint = {
        "contact": [
            {
                "entity1": 1,
                "copy1": 1,
                "position1": 10,
                "atom1": "CB",
                "entity2": 2,
                "copy2": 1,
                "position2": 20,
                "min_distance": 2.5,
                "max_distance": 6.5,
            }
        ]
    }

    feats = ConstraintFeatureGenerator.generate_tfg_contact_restraint_features(
        _token_array(), DummyAtomArray(), [{}, {}], constraint
    )

    torch.testing.assert_close(
        feats["user_distance_restraint_index"], torch.tensor([[1], [2]])
    )
    torch.testing.assert_close(
        feats["user_distance_restraint_lower_bound"], torch.tensor([2.5])
    )
    torch.testing.assert_close(
        feats["user_distance_restraint_upper_bound"], torch.tensor([6.5])
    )


def test_tfg_contact_features_deduplicate_with_tightest_interval():
    pytest.importorskip("biotite")
    from protenix.data.constraint.constraint_featurizer import ConstraintFeatureGenerator

    constraint = {
        "contact": [
            {
                "entity1": 1,
                "copy1": 1,
                "position1": 10,
                "entity2": 2,
                "copy2": 1,
                "position2": 20,
                "min_distance": 1.0,
                "max_distance": 8.0,
            },
            {
                "entity1": 2,
                "copy1": 1,
                "position1": 20,
                "entity2": 1,
                "copy2": 1,
                "position2": 10,
                "min_distance": 2.0,
                "max_distance": 5.0,
            },
        ]
    }

    feats = ConstraintFeatureGenerator.generate_tfg_contact_restraint_features(
        _token_array(), DummyAtomArray(), [{}, {}], constraint
    )

    torch.testing.assert_close(
        feats["user_distance_restraint_index"], torch.tensor([[0], [2]])
    )
    torch.testing.assert_close(
        feats["user_distance_restraint_lower_bound"], torch.tensor([2.0])
    )
    torch.testing.assert_close(
        feats["user_distance_restraint_upper_bound"], torch.tensor([5.0])
    )


def test_tfg_config_validates_user_restraint_features():
    pytest.importorskip("rdkit")
    from protenix.tfg.config import parse_tfg_config, validate_features

    cfg = parse_tfg_config(
        {
            "enable": True,
            "terms": {
                "UserDistanceRestraintPotential": {
                    "interval": 1,
                    "weight": 1.0,
                }
            },
        }
    )
    feats = _empty_user_restraint_features()

    validate_features(feats, cfg.terms)
