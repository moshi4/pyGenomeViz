import urllib.request
from importlib.util import find_spec

import pytest

from pygenomeviz.align import Blast, MMseqs, MUMmer, ProgressiveMauve

skipif_blast_not_installed = pytest.mark.skipif(
    not Blast.check_installation(exit_on_false=False),
    reason=f"{Blast.get_tool_name()} is not installed.",
)

skipif_mummer_not_installed = pytest.mark.skipif(
    not MUMmer.check_installation(exit_on_false=False),
    reason=f"{MUMmer.get_tool_name()} is not installed.",
)

skipif_mmseqs_not_installed = pytest.mark.skipif(
    not MMseqs.check_installation(exit_on_false=False),
    reason=f"{MMseqs.get_tool_name()} is not installed.",
)

skipif_pmauve_not_installed = pytest.mark.skipif(
    not ProgressiveMauve.check_installation(exit_on_false=False),
    reason=f"{ProgressiveMauve.get_tool_name()} is not installed.",
)


def _check_streamlit_installation() -> bool:
    return True if find_spec("streamlit") else False


skipif_streamlit_not_installed = pytest.mark.skipif(
    not _check_streamlit_installation(),
    reason="streamlit is not installed.",
)


def _check_network_connection() -> bool:
    try:
        urllib.request.urlopen("https://github.com")
        return True
    except Exception:
        return False


skipif_network_connection_failed = pytest.mark.skipif(
    not _check_network_connection(),
    reason="Network connection failed.",
)
