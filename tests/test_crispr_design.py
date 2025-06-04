import pytest
pytest.importorskip("Bio")
pytest.importorskip("numpy")

from Bio.Seq import Seq
from enhancement_engine.core.crispr import CRISPRDesigner
from enhancement_engine.models.data_classes import CasType


def test_design_guides_cas9_no_keyerror():
    designer = CRISPRDesigner(CasType.CAS9)
    sequence = "ATCG" * 5 + "TGG" + "A" * 5
    guides = designer.design_guides(sequence)
    assert isinstance(guides, list)
