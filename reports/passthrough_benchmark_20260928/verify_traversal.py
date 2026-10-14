"""Independent three-node BH/traversal reconstruction for the matched rerun."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parent


def main():
    rows = pd.read_csv(BASE / "evaluation/known_branch_traversal_rows.csv")
    old = pd.read_csv(BASE.parent / "known_branch_skew_20260928/evaluation/known_branch_traversal_rows.csv")
    keys = ["pattern", "covariance_profile", "replicate", "passthrough"]
    assert not rows.duplicated(keys).any()
    pd.testing.assert_frame_equal(rows.sort_values(keys).reset_index(drop=True), old.sort_values(keys).reset_index(drop=True))
    expected_null = {
        "global_null": [True, True, True], "root_only": [False, True, True],
        "A_deep_only": [True, False, True], "separate_deep_branches": [True, False, False],
        "strong_root_weak_A": [False, False, True],
    }
    paired_p = rows.groupby(keys[:3]).injected_local_tests.nunique()
    assert paired_p.eq(1).all()
    assert rows.groupby(keys[:3]).size().eq(2).all()
    audited = 0
    for row in rows.itertuples():
        assert row.status == "ok"
        tests = json.loads(row.injected_local_tests)
        p = np.array([tests[name]["p_value"] for name in ("R", "A", "B")])
        root = p[0] <= row.alpha
        child = np.ones(2, dtype=bool) if max(p[1:]) <= row.alpha else p[1:] <= row.alpha/2
        descend = bool(root or (row.passthrough and child.any()))
        split = np.array([root, descend and child[0], descend and child[1]])
        null = np.array(expected_null[row.pattern])
        false = int((split & null).sum())
        assert row.n_reached_splits == split.sum()
        assert row.n_reached_false_splits == false
        assert np.isclose(row.reached_false_split_fdp, false/max(1, split.sum()))
        assert row.n_reached_internal == (3 if descend else 1)
        assert row.num_clusters == (2+int(child.sum()) if descend else 1)
        trace = json.loads(row.traversal_decisions)
        assert trace["R"] == ("split" if root else "pass_through" if descend else "boundary")
        for i, node in enumerate(("A", "B")):
            assert (node in trace) == descend
            if descend:
                assert trace[node] == ("split" if child[i] else "boundary")
        audited += 1
    manifest = json.loads((BASE / "evaluation/manifest.json").read_text())
    frozen = json.loads((BASE / "pre_run_manifest.json").read_text())
    assert frozen["source_sha256"] == manifest["source_sha256"]
    for name, digest in manifest["source_sha256"].items():
        assert hashlib.sha256((BASE.parents[1] / name).read_bytes()).hexdigest() == digest
    result = {"passed": True, "matched_dataset_pairs": len(paired_p), "independently_reconstructed_traversals": audited, "original_rerun_rows_equal": True, "same_injected_tests_within_every_pair": True, "before_after_and_current_source_hashes_verified": len(manifest["source_sha256"]), "scope": "Fixed known tree, known-covariance Gaussian local tests and deliberately open edges; rerun is not additional independent evidence."}
    (BASE / "independent_verification.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
