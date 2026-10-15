"""PDF-focused utility helpers for benchmark plotting and report assembly."""

from tree_break_selection.plot.backend import configure_matplotlib_backend

configure_matplotlib_backend()

from .layout import PDF_PAGE_SIZE_INCHES, prepare_pdf_figure
from .merge import merge_existing_pdfs
from .session import open_pdf_pages, resolve_pdf_output_path

__all__ = [
    "PDF_PAGE_SIZE_INCHES",
    "prepare_pdf_figure",
    "merge_existing_pdfs",
    "resolve_pdf_output_path",
    "open_pdf_pages",
]
