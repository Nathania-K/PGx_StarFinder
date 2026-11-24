"""Command-line utility for running Django management commands."""

import os 
import sys

def main():
    """
    Django's management command:
    1. Specifies project settings module
    2. Imports Django's managment framework
    3. Delegates execution to Django

    Commands supplied on command line is forwareded unchanged to Django's processing.
    """

    #Tells Django which modules to load.
    os.environ.setdefault(
        "DJANGO_SETTINGS_MODULE",
        "pgx_starfinder.settings",
    )

    #Imports Django's management command framework.
    try:
        from django.core.management import execute_from_command_line

    except ImportError as exc:
        raise ImportError(
            "Couldn't import Django." 
            "Check if installed and if Conda enviroment is activated "
        ) from exc

    #Pass command-line args directly to Django.
    execute_from_command_line(sys.argv)

if __name__ == "__main__":
    main()