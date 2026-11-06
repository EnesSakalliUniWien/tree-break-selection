"""End-to-end identity-covariance root diagnostic, including unresolved output."""

import json
from decimal import Decimal

import numpy as np
import pandas as pd

from benchmarks.diagnostics.calibration.selected.hierarchy.diffusion_root_diagnostic import (
    diagnose_diffusion_root,
    write_root_diagnostics,
)


def test_reentry_report_keeps_observations_separate_from_certificates(tmp_path):
    values = np.random.default_rng(20261105).normal(size=(9, 1))
    result = diagnose_diffusion_root(values, anchor_ratios=(1.0, 1.1, 1.2))
    assert [row["status"] for row in result["replay"]] == ["selected", "nonselected", "selected"]
    root = result["root"]
    assert root["feature_dimension"] == 1
    assert root["covariance"] == "known_identity"
    assert root["projection"] == "full_feature"
    assert root["branch_length"] == 0
    assert root["inflation"] == 1
    assert root["contrast_energy_rank_selection"] is False
    assert root["passthrough"] is False
    assert Decimal(root["selection_mass_lower"]) > 0
    assert Decimal(0) <= Decimal(root["p_lower"]) <= Decimal(root["p_upper"]) <= Decimal(1)
    assert root["gate_status"] == "unresolved"
    assert root["clusters"] == 1
    outputs = write_root_diagnostics([result], tmp_path)
    manifest = json.loads(outputs["manifest"].read_text())
    assert manifest["root_counts"]["unresolved"] == 1
    assert manifest["scope"] == "diagnostic_root_only"
    assert "not bitwise" in manifest["arithmetic_contract"]
    assert manifest["source_hashes"]
    assert pd.read_csv(outputs["cells"]).status.eq("selected").any()
    assert pd.read_csv(outputs["roots"]).clusters.tolist() == [1]


def test_unsupported_observation_is_reported_without_opening_root():
    result = diagnose_diffusion_root(np.zeros((9, 1)))
    assert result["root"]["gate_status"] == "unsupported"
    assert result["root"]["clusters"] == 1
    assert result["root"]["passthrough"] is False
    assert result["root"]["p_lower"] == "0"
    assert result["root"]["p_upper"] == "1"
    assert result["root"]["error"]


def test_unsupported_panel_still_writes_readable_empty_tables(tmp_path):
    result = diagnose_diffusion_root(np.zeros((9, 1)))
    outputs = write_root_diagnostics([result], tmp_path)
    assert pd.read_csv(outputs["cells"]).empty
    assert pd.read_csv(outputs["replay"]).empty
    manifest = json.loads(outputs["manifest"].read_text())
    assert manifest["source_hashes"]
    assert len(result["root"]["input_data_sha256"]) == 64


def test_report_exposes_certified_window_expansion_without_dropping_unknown_mass():
    values = np.random.default_rng(20261105).normal(size=(9, 1))
    result = diagnose_diffusion_root(values, observed_half_width=1e-8, expansion_steps=2)
    root = result["root"]
    assert root["observed_half_width"] == 1e-8
    assert root["expansion_steps"] == 2
    assert Decimal(root["selection_mass_lower"]) > Decimal("1e-10")
    assert root["gate_status"] == "unresolved"
    assert root["clusters"] == 1
    assert result["cells"][-1]["status"] == "unresolved"
    assert result["cells"][-1]["hi"] == "inf"
