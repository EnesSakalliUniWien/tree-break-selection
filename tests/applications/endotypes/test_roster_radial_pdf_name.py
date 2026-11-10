from applications.endotypes.reports.export_subspace_cluster_rosters import radial_pdf_name


def test_radial_pdf_name_keeps_historical_form() -> None:
    assert radial_pdf_name("My Data/1_") == "my_data_1_radial_tree_clusters.pdf"
    assert radial_pdf_name("_x_") == "x_radial_tree_clusters.pdf"
