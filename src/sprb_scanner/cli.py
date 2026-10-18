"""Thin CLI boundary; scanning commands are added in subsequent batches."""

from pathlib import Path
from typing import Annotated

import typer

from sprb_scanner import __version__
from sprb_scanner.core.config import load_config
from sprb_scanner.core.exceptions import ConfigurationError

app = typer.Typer(no_args_is_help=True, help="SprB-like Scanner V2 — project foundation.")


@app.command()
def version():
    """Show the installed package version."""
    typer.echo(__version__)


@app.command()
def validate(config: Annotated[Path | None, typer.Option(help="YAML configuration path")] = None):
    """Validate configuration only; input/model/result validation arrives in later batches."""
    try:
        settings = load_config(config)
    except ConfigurationError as exc:
        typer.echo(str(exc), err=True)
        raise typer.Exit(2) from exc
    typer.echo(settings.model_dump_json(indent=2))
