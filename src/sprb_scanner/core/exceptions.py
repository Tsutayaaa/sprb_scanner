"""Errors suitable for reporting at application boundaries."""


class ScannerError(Exception):
    """Base application error."""


class ConfigurationError(ScannerError):
    """Invalid or unreadable configuration."""


class InputError(ScannerError):
    """Input is invalid or outside the documented contract."""


class DownloadError(ScannerError):
    """Download or archive validation failed."""


class HMMError(ScannerError):
    """Profile loading, provenance or HMM search failed."""


class MMseqsError(ScannerError):
    """MMseqs2 database construction, execution or alignment validation failed."""


class ReferenceRegistryError(ScannerError):
    """Invalid resource identity, dependency graph, checksum or artifact format."""


class ReferenceUnavailableError(ReferenceRegistryError):
    """A registered resource has no local artifact; it has not been downloaded."""
