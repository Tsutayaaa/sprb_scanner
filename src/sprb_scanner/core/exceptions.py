"""Errors suitable for reporting at application boundaries."""


class ScannerError(Exception):
    """Base application error."""


class ConfigurationError(ScannerError):
    """Invalid or unreadable configuration."""


class InputError(ScannerError):
    """Input is invalid or outside the documented contract."""


class DownloadError(ScannerError):
    """Download or archive validation failed."""
