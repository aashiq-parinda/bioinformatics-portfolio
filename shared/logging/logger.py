"""Structured and styled logging module using Rich and standard logging."""

import logging
import os
import sys
from typing import Optional

from rich.console import Console
from rich.logging import RichHandler

_LOGGER_CACHE: dict[str, logging.Logger] = {}


def setup_logger(
    name: str = "bioinformatics",
    level: Optional[str] = None,
    log_file: Optional[str] = None,
) -> logging.Logger:
    """Configure and return a structured logger with rich console output and optional file handler.

    Args:
        name: Name of the logger namespace.
        level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL).
        log_file: Optional path to append structured log records.

    Returns:
        Configured logging.Logger instance.
    """
    if name in _LOGGER_CACHE:
        return _LOGGER_CACHE[name]

    log_level_str = (level or os.getenv("LOG_LEVEL", "INFO")).upper()
    log_level = getattr(logging, log_level_str, logging.INFO)

    logger = logging.getLogger(name)
    logger.setLevel(log_level)
    logger.propagate = False

    # Clear existing handlers to prevent duplicate output
    if logger.hasHandlers():
        logger.handlers.clear()

    # Rich Console Handler
    console = Console(file=sys.stderr, color_system="auto")
    rich_handler = RichHandler(
        console=console,
        show_time=True,
        show_path=False,
        rich_tracebacks=True,
        markup=True,
    )
    rich_handler.setLevel(log_level)
    logger.addHandler(rich_handler)

    # Optional File Handler
    if log_file:
        os.makedirs(os.path.dirname(os.path.abspath(log_file)), exist_ok=True)
        file_handler = logging.FileHandler(log_file, encoding="utf-8")
        file_handler.setLevel(log_level)
        file_formatter = logging.Formatter(
            "%(asctime)s [%(levelname)s] %(name)s: %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)

    _LOGGER_CACHE[name] = logger
    return logger


def get_logger(name: str = "bioinformatics") -> logging.Logger:
    """Retrieve or create an existing logger namespace."""
    if name in _LOGGER_CACHE:
        return _LOGGER_CACHE[name]
    return setup_logger(name=name)
