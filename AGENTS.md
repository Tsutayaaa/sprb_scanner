# Development rules

- Use src layout, Python 3.11+, uv and the committed lockfile.
- Core must not import scanners. Scanners share core records and remain independent.
- Keep scanning logic out of CLI and visualization. Configure significant parameters.
- Internal genomic and protein coordinates are 0-based, half-open.
- Preserve raw hits, provenance and failures; never silently discard low-confidence evidence.
- Never compare or combine HMM and sequence-search scores as interchangeable statistics.
- No hard full-protein length or full-length coverage filters.
- Each new feature needs meaningful tests. Public model changes require checking all consumers.
- Document Purpose, Input, Output, Dependencies, Parameters, Error Handling and Test Cases.
- Use a separate codex/ feature branch for each independently reviewable batch; never commit to main.
- Submit changes through PRs. Never bypass CI. Merge commits only; no squash or rebase merges.
- PRs require explicit confirmation by Tsutayaaa. Do not merge without that confirmation.
- Do not add a license until the owner requests one.
- Do not commit large genome data, downloaded models or generated results; tiny fixtures are allowed.
- Future modules stay documented interfaces until their development phase.

Checks: uv sync --locked; uv run pytest; uv run ruff check .; uv run ruff format --check .;
uv run sprb-scanner --help; uv build.
