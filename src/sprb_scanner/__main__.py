"""Allow the installed interpreter to invoke the same public CLI entry point."""

from sprb_scanner.cli import app

if __name__ == "__main__":
    app()
