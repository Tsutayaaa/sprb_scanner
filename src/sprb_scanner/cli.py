"""Thin CLI adapters for offline preparation and independent NCBI acquisition."""

import json
from pathlib import Path
from typing import Annotated

import typer

from sprb_scanner import __version__
from sprb_scanner.core.config import load_config, load_profile
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


@app.command("hmm-search")
def hmm_search(
    prepared: Annotated[Path, typer.Option(help="Validated prepared dataset directory")],
    profile: Annotated[Path, typer.Option(help="ProfileConfig JSON with provenance and HMM path")],
    output: Annotated[Path, typer.Option(help="New scan output directory")],
    threads: int = 1,
    sequence_evalue: float = 10.0,
    domain_evalue: float = 10.0,
    max_sensitivity: bool = False,
):
    """Search one protein HMM; export all returned domains without classification."""
    from sprb_scanner.core.models import GenomeDataset
    from sprb_scanner.scanners.hmm import HMMSearchConfig, run_hmm

    try:
        validation = validate_prepared(prepared)
        if validation["errors"]:
            raise ValueError("Prepared dataset has quality errors")
        dataset = GenomeDataset.model_validate_json((prepared / "dataset.json").read_text())
        model = load_profile(profile)
        settings = HMMSearchConfig(
            threads=threads,
            sequence_report_evalue=sequence_evalue,
            domain_report_evalue=domain_evalue,
            max_sensitivity=max_sensitivity,
        )
        emit(run_hmm(dataset.proteins, model, output, settings))
    except (ScannerError, OSError, ValueError) as exc:
        fail(exc)


@app.command("sprf-search")
def sprf_search(
    prepared: Annotated[Path, typer.Option(help="Validated prepared dataset directory")],
    tigr_profile: Annotated[Path, typer.Option(help="TIGR03519 ProfileConfig JSON")],
    pfam_profile: Annotated[Path, typer.Option(help="PF11751 ProfileConfig JSON")],
    output: Annotated[Path, typer.Option(help="New PorP/SprF output directory")],
    threads: int = 1,
    sequence_evalue: float = 10.0,
    domain_evalue: float = 10.0,
    max_sensitivity: bool = False,
):
    """Scan PorP/SprF family with TIGR03519 and PF11751 independently."""
    from sprb_scanner.core.models import GenomeDataset
    from sprb_scanner.scanners.hmm import HMMSearchConfig
    from sprb_scanner.scanners.sprf import run_sprf

    try:
        validation = validate_prepared(prepared)
        if validation["errors"]:
            raise ValueError("Prepared dataset has quality errors")
        dataset = GenomeDataset.model_validate_json((prepared / "dataset.json").read_text())
        settings = HMMSearchConfig(
            threads=threads,
            sequence_report_evalue=sequence_evalue,
            domain_report_evalue=domain_evalue,
            max_sensitivity=max_sensitivity,
        )
        report = run_sprf(
            dataset.proteins,
            (load_profile(tigr_profile), load_profile(pfam_profile)),
            output,
            settings,
        )
        emit(report)
        if report["state"] != "succeeded":
            raise typer.Exit(2)
    except (ScannerError, OSError, ValueError) as exc:
        fail(exc)


@app.command("sprb-repeat-search")
def sprb_repeat_search(
    prepared: Annotated[Path, typer.Option(help="Validated prepared dataset directory")],
    profile: Annotated[Path, typer.Option(help="PF13573 SprB repeat profile JSON")],
    output: Annotated[Path, typer.Option(help="New output directory")],
    threads: int = 1,
    sequence_evalue: float = 10.0,
    domain_evalue: float = 10.0,
    max_sensitivity: bool = False,
):
    """Scan PF13573 and retain all internal repeat domain alignments."""
    from sprb_scanner.core.models import GenomeDataset
    from sprb_scanner.scanners.hmm import HMMSearchConfig
    from sprb_scanner.scanners.repeats import run_repeats

    try:
        if validate_prepared(prepared)["errors"]:
            raise ValueError("Prepared dataset has quality errors")
        dataset = GenomeDataset.model_validate_json((prepared / "dataset.json").read_text())
        settings = HMMSearchConfig(
            threads=threads,
            sequence_report_evalue=sequence_evalue,
            domain_report_evalue=domain_evalue,
            max_sensitivity=max_sensitivity,
        )
        emit(run_repeats(dataset.proteins, load_profile(profile), output, settings))
    except (ScannerError, OSError, ValueError) as exc:
        fail(exc)


@app.command("bctd-hmm-search")
def bctd_hmm_search(
    prepared: Annotated[Path, typer.Option(help="Validated prepared dataset directory")],
    tigr_profile: Annotated[Path, typer.Option(help="TIGR04131 profile JSON")],
    pfam_profile: Annotated[Path, typer.Option(help="PF13585 CHU_C profile JSON")],
    output: Annotated[Path, typer.Option(help="New output directory")],
    threads: int = 1,
    sequence_evalue: float = 10.0,
    domain_evalue: float = 10.0,
    max_sensitivity: bool = False,
):
    """Scan Type B CTD with TIGR04131 and PF13585 independently."""
    from sprb_scanner.core.models import GenomeDataset
    from sprb_scanner.scanners.bctd import run_bctd
    from sprb_scanner.scanners.hmm import HMMSearchConfig

    try:
        if validate_prepared(prepared)["errors"]:
            raise ValueError("Prepared dataset has quality errors")
        dataset = GenomeDataset.model_validate_json((prepared / "dataset.json").read_text())
        settings = HMMSearchConfig(
            threads=threads,
            sequence_report_evalue=sequence_evalue,
            domain_report_evalue=domain_evalue,
            max_sensitivity=max_sensitivity,
        )
        report = run_bctd(
            dataset.proteins,
            (load_profile(tigr_profile), load_profile(pfam_profile)),
            output,
            settings,
        )
        emit(report)
        if report["state"] != "succeeded":
            raise typer.Exit(2)
    except (ScannerError, OSError, ValueError) as exc:
        fail(exc)


@app.command("neighborhood")
def neighborhood(
    prepared: Annotated[Path, typer.Option(help="Prepared annotated dataset directory")],
    gene: Annotated[str, typer.Option(help="Internal CDS gene_id")],
    output: Annotated[Path, typer.Option(help="New neighborhood output directory")],
    window: Annotated[int, typer.Option(min=0, help="CDS ranks on each side")] = 5,
):
    """Query objective CDS neighbors in coordinate and target transcription directions."""
    from sprb_scanner.neighborhood.service import run_neighborhood

    try:
        result = run_neighborhood(prepared, gene, output, window)
        typer.echo(
            f"{gene} | {result.contig.source_id} | strand {result.target.strand} | "
            f"{len(result.neighbors)} neighbors | 0-based half-open"
        )
        for neighbor in result.neighbors:
            typer.echo(
                f"{neighbor.gene.gene_id} | {neighbor.relative_position or 'unknown'} "
                f"(genomic {neighbor.genomic_position}) | "
                f"{neighbor.intervening_cds} intervening CDS | {neighbor.distance_bp} bp | "
                f"overlap={neighbor.overlap_bp} bp signed_gap={neighbor.signed_gap_bp} bp | "
                f"adjacent={neighbor.directly_adjacent} same_strand={neighbor.same_strand}"
            )
        typer.echo(f"Output: {output.resolve()}")
    except (ScannerError, OSError, ValueError) as exc:
        fail(exc)


@app.command("typeb-sprf-associate")
def typeb_sprf_associate(
    prepared: Annotated[Path, typer.Option(help="Prepared annotated genome")],
    bctd_scan: Annotated[Path, typer.Option(help="Complete bctd-hmm-search output")],
    sprf_scan: Annotated[Path, typer.Option(help="Complete sprf-search output")],
    output: Annotated[Path, typer.Option(help="New association output directory")],
    window: Annotated[int, typer.Option(min=0, help="CDS ranks on each side")] = 5,
):
    """Associate existing Type B / PorP-SprF HMM domains by local genomic context."""
    from sprb_scanner.association.service import run_association

    try:
        emit(run_association(prepared, bctd_scan, sprf_scan, output, window))
    except (ScannerError, OSError, ValueError) as exc:
        fail(exc)
