# Model provenance

No biological profiles are bundled in the foundation batch. Configure model records under
`bctd.profiles` in the YAML config; each requires provenance (model_id, version, source, sha256),
method (`hmm` or `mmseqs2`) and a local path. Relative paths resolve against the YAML directory.
Only structural configuration validation currently runs; file existence/checksum/profile identity
validation is implemented with the scanner in batch 3. Never invent a checksum or model version.
