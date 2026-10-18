"""Thin CLI adapters for offline preparation and independent NCBI acquisition."""

import json
from pathlib import Path
from typing import Annotated

import typer

from sprb_scanner import __version__
from sprb_scanner.core.config import load_config
from sprb_scanner.core.exceptions import ScannerError
from sprb_scanner.io.exporters import summary, validate_prepared
from sprb_scanner.io.ncbi import download as download_genome
from sprb_scanner.io.prepare import prepare as prepare_genome
from sprb_scanner.io.prepare import prepare_dataset

app = typer.Typer(no_args_is_help=True, help="SprB-like Scanner V2 — genome infrastructure.")
ConfigPath = Annotated[Path | None, typer.Option(help="YAML configuration path")]
InputPath = Annotated[Path | None, typer.Option("--input", help="Input file")]
Format = Annotated[str, typer.Option("--format", help="genbank, gff3 or protein_fasta")]


def fail(exc):
    typer.echo(str(exc), err=True)
    raise typer.Exit(2) from exc


def emit(value):
    typer.echo(json.dumps(value, indent=2, ensure_ascii=False))


@app.command()
def version():
    """Show the installed package version."""
    typer.echo(__version__)


@app.command()
def validate(
    config: ConfigPath = None,
    input: InputPath = None,
    format: Format = "genbank",
    fna: Path | None = None,
    faa: Path | None = None,
    genome_id: str | None = None,
    prepared: Path | None = None,
):
    """Validate configuration, raw inputs or a prepared output directory (offline)."""
    try:
        settings = load_config(config)
        if input is not None and prepared is not None:
            raise ValueError("Choose --input or --prepared, not both")
        if input is None and (fna is not None or faa is not None or genome_id is not None):
            raise ValueError("--fna, --faa and --genome-id require --input")
        if input is not None:
            dataset, _ = prepare_dataset(
                input, format, genome_id, fna, faa, settings.input.default_translation_table
            )
            emit(
                summary(dataset)
                | {"quality_issues": [x.model_dump() for x in dataset.quality_issues]}
            )
            if any(x.severity == "error" for x in dataset.quality_issues):
                raise typer.Exit(2)
        elif prepared is not None:
            result = validate_prepared(prepared)
            emit(result)
            if result["errors"]:
                raise typer.Exit(2)
        else:
            typer.echo(settings.model_dump_json(indent=2))
    except (ScannerError, OSError, ValueError) as exc:
        fail(exc)


@app.command()
def prepare(
    input: Annotated[Path, typer.Option("--input")],
    format: Format,
    output: Annotated[Path, typer.Option(help="New output directory")],
    config: ConfigPath = None,
    genome_id: str | None = None,
    fna: Path | None = None,
    faa: Path | None = None,
):
    """Normalize one annotated assembly or protein FASTA; no scanning is performed."""
    try:
        result = prepare_genome(input, format, output, load_config(config), genome_id, fna, faa)
        emit(result)
    except (ScannerError, OSError, ValueError) as exc:
        fail(exc)


@app.command()
def download(
    accession: Annotated[str, typer.Argument(help="One versioned GCF/GCA accession")],
    output: Annotated[Path, typer.Option(help="New package directory")],
    config: ConfigPath = None,
):
    """Download one NCBI Datasets assembly package with provenance and checksums."""
    try:
        manifest = download_genome(accession, output, load_config(config).download)
        emit(
            {
                "accession": manifest["accession"],
                "output": str(output.resolve()),
                "files": len(manifest["files"]),
            }
        )
    except (ScannerError, OSError, ValueError) as exc:
        fail(exc)
