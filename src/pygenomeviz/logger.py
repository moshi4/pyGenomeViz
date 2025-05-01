from __future__ import annotations

import logging
import sys
from pathlib import Path


def init_null_logger():
    """Initialize package root logger with NullHandler

    Configuring package root null logger for a library
    https://docs.python.org/3/howto/logging.html#configuring-logging-for-a-library
    """
    pkg_root_name = __name__.split(".")[0]
    logger = logging.getLogger(pkg_root_name)
    logger.addHandler(logging.NullHandler())


def init_logger(
    *,
    quiet: bool = False,
    verbose: bool = False,
    log_file: str | Path | None = None,
):
    """Initialize package root logger with StreamHandler(& FileHandler)

    Configuring package root default logger for a CLI tool

    Parameters
    ----------
    quiet : bool, optional
        If True, no print info log on screen
    verbose: bool, optional
        If True & quiet=False, print debug log on screen
    log_file : str | Path | None, optional
        Log file
    """
    pkg_root_name = __name__.split(".")[0]
    logger = logging.getLogger(pkg_root_name)

    # Remove existing handler to avoid duplicate logging
    for handler in logger.handlers:
        logger.removeHandler(handler)
        handler.close()

    logger.setLevel(logging.DEBUG)
    log_formatter = logging.Formatter(
        fmt="$asctime | $levelname | $message",
        datefmt="%Y-%m-%d %H:%M:%S",
        style="$",
    )
    # Add stream handler for terminal stderr
    stream_handler = logging.StreamHandler(sys.stderr)
    stream_handler.setFormatter(log_formatter)
    if quiet:
        log_level = logging.WARNING
    else:
        log_level = logging.DEBUG if verbose else logging.INFO
    stream_handler.setLevel(log_level)
    logger.addHandler(stream_handler)

    if log_file:
        # Add file handler for log file
        file_handler = logging.FileHandler(log_file, mode="w", encoding="utf-8")
        file_handler.setFormatter(log_formatter)
        log_level = logging.DEBUG if verbose else logging.INFO
        file_handler.setLevel(log_level)
        logger.addHandler(file_handler)
