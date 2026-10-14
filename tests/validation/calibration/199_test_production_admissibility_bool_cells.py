from __future__ import annotations

import pandas as pd
from benchmarks.diagnostics.calibration.traversal.production_admissibility_contract import (
    evaluate_production_admissibility_components,
)


def test_required_for_production_reads_float_typed_flags() -> None:
    # pandas stores a bool column as float once it has held a blank cell.
    components = pd.DataFrame(
        {
            "contract_id": ["c", "c"],
            "component_id": ["a", "b"],
            "component_type": ["gate", "gate"],
            "component_status": ["production_ready", "fail_closed_undefined"],
            "required_for_production": [1.0, 0.0],
        }
    )
    normalized = evaluate_production_admissibility_components(components)
    assert normalized["required_for_production"].tolist() == [True, False]
