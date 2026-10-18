"""Errors suitable for reporting at application boundaries."""


class ScannerError(Exception):
    """Base application error."""


class ConfigurationError(ScannerError):
    """Invalid or unreadable configuration."""
