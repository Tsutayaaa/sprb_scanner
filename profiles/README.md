# Model provenance

Biological HMM files remain outside Git. `TIGR04131.profile.json` is a checked provenance
manifest, not the model itself. Copy it next to the downloaded `TIGR04131.1.HMM`, or change
its `path` to the local absolute path. Relative paths resolve against the manifest directory.

`hmm-search --profile FILE` accepts one ProfileConfig JSON containing method, path and
provenance (model_id, version, source, sha256). The engine validates checksum, HMM syntax,
amino-acid alphabet, single-model count, identity and versioned accession before searching.
The generic engine accepts other amino-acid models with their own manifests; do not invent
checksums or substitute models under the TIGR04131 identity.

TIGR04131.1 was acquired from NCBI PGAP HMM release 20.0. The `current` URL is mutable:
the recorded SHA-256 pins the actual bytes. No models are redistributed by this repository.
See [HMM engine usage and validation](../docs/hmm_engine.md).

The older `bctd.profiles` YAML field remains a configuration interface for future orchestration;
this single-model command uses its explicit JSON manifest, not that YAML field.

PorP/SprF family manifests: TIGR03519.profile.json and PF11751.profile.json. Models are
not bundled; copy manifests next to the corresponding downloaded files or set absolute paths.
The PF11751 download is gzip; its manifest pins the decompressed HMM bytes. See
[sprf-search](../docs/sprf_scanner.md) for dual-model evidence and failure semantics.

SprB/Pfam manifests: PF13573.profile.json (SprB repeat) and PF13585.profile.json (CHU_C).
Both source downloads are gzip and checksums refer to decompressed HMM files. See
[repeat and Type B CTD HMM usage](../docs/sprb_pfam_models.md). Models stay outside Git.
