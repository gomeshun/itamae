"""Composable variance adapters, integration, and explicit safe caching.

Existing SASHIMI implementations can use the callable adapter during migration.
WDM and FDM variants can compose their model-owned spectra with ITAMAE windows,
integration, and cache-key machinery without moving physical defaults here.
"""

from .cache import (
    VARIANCE_CACHE_SCHEMA_VERSION,
    load_variance_cache,
    save_variance_cache,
    variance_cache_key,
)
from .callable import CallableVarianceModel
from .integrated import IntegratedVarianceModel

__all__ = [
    "CallableVarianceModel",
    "IntegratedVarianceModel",
    "VARIANCE_CACHE_SCHEMA_VERSION",
    "load_variance_cache",
    "save_variance_cache",
    "variance_cache_key",
]
