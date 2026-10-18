"""Compatibility entry point; prefer the installed sprb-scanner command."""

from sprb_scanner.cli import app

if __name__ == "__main__":
    app()
