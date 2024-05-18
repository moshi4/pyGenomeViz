from __future__ import annotations

import logging
import sys
from pathlib import Path


def get_logger(
    name: str | None = None,
    log_file: str | Path | None = None,
    quiet: bool = False,
) -> logging.Logger:
    """Get logger

    Parameters
    ----------
    name : str | None, optional
        Logger name
    log_file : str | Path | None, optional
        Log file path for log file stream
    quiet : bool, optional
        If True, don't display log message less than WARNING level (DEBUG, INFO).

    Returns
    -------
    logger : logging.Logger
        Logger
    """
    logger = logging.getLogger(name)

    # If logger already exists, only apply quiet
    if logger.hasHandlers():
        for handler in logger.handlers:
            if isinstance(handler, logging.StreamHandler):
                log_level = logging.WARNING if quiet else logging.INFO
                handler.setLevel(log_level)
        return logger

    logger.setLevel(logging.DEBUG)

    # Setup log formatter
    format = "$asctime | $levelname | $message"
    datefmt = "%Y-%m-%d %H:%M:%S"
    formatter = logging.Formatter(fmt=format, datefmt=datefmt, style="$")

    # Setup stream handler
    stream_handler = logging.StreamHandler(sys.stderr)
    stream_handler.setFormatter(formatter)
    log_level = logging.WARNING if quiet else logging.INFO
    stream_handler.setLevel(log_level)
    logger.addHandler(stream_handler)

    if log_file:
        # Setup log file handler
        file_handler = logging.FileHandler(log_file, mode="w")
        file_handler.setFormatter(formatter)
        file_handler.setLevel(logging.DEBUG)
        logger.addHandler(file_handler)

    return logger
