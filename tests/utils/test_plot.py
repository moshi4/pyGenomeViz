from __future__ import annotations

from pathlib import Path

from pygenomeviz.gui.config import (
    AlignConfig,
    FeatureConfig,
    FigureConfig,
    PgvGuiPlotConfig,
)
from pygenomeviz.gui.plot import plot_by_gui_cfg
from pygenomeviz.parser import Genbank
from tests.marker import skipif_blast_not_installed


@skipif_blast_not_installed
def test_create_genomeviz(gbk_dataset_files: list[Path], tmp_path: Path):
    """Test create genomeviz for GUI"""
    gbk_list = [Genbank(gbk_file) for gbk_file in gbk_dataset_files]

    name2seqid2range = {}
    for gbk in gbk_list:
        seqid2range = {}
        for seqid, size in gbk.get_seqid2size().items():
            seqid2range[seqid] = (0, size)
        name2seqid2range[gbk.name] = seqid2range

    cfg = PgvGuiPlotConfig(
        fig=FigureConfig(),
        feat=FeatureConfig(),
        aln=AlignConfig(method="BLAST (nucleotide)"),
        name2seqid2range=name2seqid2range,
    )

    gv, align_coords = plot_by_gui_cfg(gbk_list, cfg)
    outfile = tmp_path / "test.png"
    gv.savefig(outfile)

    assert len(align_coords) > 0
    assert outfile.exists() is True
