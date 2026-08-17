import pandas as pd
import pytest

from perturbseq_pipeline.design_gate import evaluate_design_gate, require_design_gate


def test_design_gate_passes_crossed_batch_design():
    metadata = pd.DataFrame(
        {
            "cell_barcode": ["c1", "c2", "c3", "c4"],
            "batch": ["B1", "B1", "B2", "B2"],
            "donor": ["D1", "D2", "D1", "D2"],
            "assigned_guide": ["NT_g1", "TP53_g1", "NT_g1", "TP53_g1"],
            "target_gene": ["NT", "TP53", "NT", "TP53"],
            "escaped_perturbation": [False, True, True, False],
            "include_primary_analysis": [True] * 4,
        }
    )

    diagnostics, decision = evaluate_design_gate(metadata)

    assert decision["downstream_analysis_allowed"] is True
    assert not diagnostics["completely_confounded"].any()


def test_design_gate_blocks_exact_batch_alias_and_enforces_decision():
    metadata = pd.DataFrame(
        {
            "cell_barcode": ["c1", "c2", "c3", "c4"],
            "batch": ["B1", "B1", "B2", "B2"],
            "donor": ["D1", "D1", "D2", "D2"],
            "assigned_guide": ["NT_g1", "NT_g1", "TP53_g1", "TP53_g1"],
            "target_gene": ["NT", "NT", "TP53", "TP53"],
            "escaped_perturbation": [False, False, True, True],
            "include_primary_analysis": [True] * 4,
        }
    )

    _, decision = evaluate_design_gate(metadata)

    assert decision["status"] == "blocked"
    assert set(decision["blocking_confounders"]) == {
        "donor",
        "assigned_guide",
        "target_gene",
        "escaped_perturbation",
    }
    with pytest.raises(RuntimeError, match="Downstream analysis blocked"):
        require_design_gate(decision)
