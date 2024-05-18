from pathlib import Path

from streamlit.testing.v1 import AppTest

import pygenomeviz.gui
from tests.marker import skipif_streamlit_not_installed


@skipif_streamlit_not_installed
def test_gui_app():
    """Test streamlit from app file"""
    app_file = Path(pygenomeviz.gui.__file__).parent / "app.py"
    app_test = AppTest.from_file(str(app_file))
    app_test.run()
    assert not app_test.exception
