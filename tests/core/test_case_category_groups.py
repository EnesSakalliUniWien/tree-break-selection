from benchmarks.shared.cases import get_default_test_cases
from benchmarks.shared.plots.cover_page import category_group


def test_every_default_case_category_has_a_report_group() -> None:
    # The relationship analysis maps each case category to a report section;
    # a category without a group crashes the full benchmark after all cases ran.
    for case in get_default_test_cases():
        assert category_group(case["category"])
