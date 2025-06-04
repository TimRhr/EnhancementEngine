import pytest
pytest.importorskip("Bio")

from enhancement_engine.cli import _parse_gene_list


def test_parse_gene_list_string():
    result = _parse_gene_list('GENE1,GENE2,GENE3')
    assert result == ['GENE1', 'GENE2', 'GENE3']


def test_parse_gene_list_file(tmp_path):
    f = tmp_path / 'genes.txt'
    f.write_text('A\nB\n')
    assert _parse_gene_list(str(f)) == ['A', 'B']
