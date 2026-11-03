"""
Defines requirements of logging configuration
"""

import os
from pathlib import Path

# Define the default log file (hidden file within home directory)
home_dir = os.path.expanduser("~")
DEFAULT_LOG = os.path.join(home_dir, ".PGx_StarFinder.log")

LOG_FILE = DEFAULT_LOG

# Normalise the path to an absolute path to ensure consistency across platforms
LOG_FILE = os.path.abspath(os.path.expanduser(LOG_FILE))

# Check the directory exists before creating logfile
log_dir = os.path.dirname(LOG_FILE)
if log_dir:
    os.makedirs(log_dir, exist_ok=True)

# Logging configuration
LOGGING = {
    "version": 1,
    "disable_existing_loggers": False,

    # --------------------------------------------------
    # FORMATTERS
    # --------------------------------------------------
    "formatters": {
        "standard": {
            "format": (
                "%(asctime)s | %(levelname)s | %(name)s | "
                "%(filename)s:%(lineno)d | %(message)s"
            )
        }
    },

    # --------------------------------------------------
    # HANDLERS
    # --------------------------------------------------
    "handlers": {
        "console": {
            "class": "logging.StreamHandler",
            "level": "INFO",
            "formatter": "standard",
        },

        # --------------------------------------------------
        # ROTATING FILE HANDLER
        # --------------------------------------------------
        #
        # This replaces the basic FileHandler with a
        # RotatingFileHandler to prevent unlimited growth.
        #

        "file": {
            "class": "logging.handlers.RotatingFileHandler",
            "level": "DEBUG",
            "formatter": "standard",
            "filename": LOG_FILE,
            # --------------------------------------------------
            # ROTATION SETTINGS (IMPORTANT TEACHING POINT)
            # --------------------------------------------------
            #
            # maxBytes:
            # Maximum size of the log file BEFORE rotation happens.
            #
            # Example:
            # 1_048_576 bytes = 1 MB
            #
            # Here we set a small size for demonstration.
            #
            "maxBytes": 1024 * 100,
            "backupCount": 3,

            # ensures file opens safely even if reused
            "encoding": "utf-8",
        },
    },

    "loggers": {
        "pgx_starfinder": {
            "level": "DEBUG",
            "handlers": ["console", "file"],
            "propagate": False
        },
    },

    "root": {
        "level": "DEBUG",
        "handlers": ["console", "file"],
    }
}