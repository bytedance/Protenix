import math

import numpy as np

from protenix.nanobody import dockq
from protenix.nanobody.metrics import aligned_rmsd, rmsd, summarize_rows


def test_simple_rmsd_calculation_on_synthetic_coordinates():
    true = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    pred = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
    assert math.isclose(rmsd(pred, true), math.sqrt(0.5))


def test_aligned_rmsd_removes_rigid_translation():
    true = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
    pred = true + np.array([10.0, -2.0, 5.0])
    assert aligned_rmsd(pred, true) < 1e-6


def test_dockq_unavailable_is_graceful(monkeypatch, tmp_path):
    monkeypatch.setattr(dockq, "find_dockq_executable", lambda: None)
    result = dockq.compute_dockq(tmp_path / "pred.cif", tmp_path / "native.cif")
    assert result["metric_status"] == "dockq_unavailable"
    assert result["dockq"] is None
    assert "DockQ executable not found" in result["notes"]


def test_summary_aggregates_counts_and_nanobody_metrics():
    rows = [
        {
            "metric_status": "ok",
            "dockq": 0.5,
            "dockq_success_0_23": True,
            "cdrh3_rmsd_ca": 1.0,
            "nanobody_rmsd_ca": 2.0,
        },
        {
            "metric_status": "missing_prediction",
            "dockq": "",
            "cdrh3_rmsd_ca": 3.0,
            "nanobody_rmsd_ca": 4.0,
        },
    ]
    summary = summarize_rows(rows)
    assert summary["n_total"] == 2
    assert summary["n_successful"] == 1
    assert summary["n_missing_prediction"] == 1
    assert summary["dockq_success_rate_0_23"] == 1.0
    assert summary["median_cdrh3_rmsd_ca"] == 2.0
