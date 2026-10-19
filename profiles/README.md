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
