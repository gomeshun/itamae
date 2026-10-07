"""Explicit, content-addressed variance-table cache helpers."""

from __future__ import annotations

from hashlib import sha256
import json
import os
from pathlib import Path
import tempfile
from typing import Any, Mapping

import numpy as np

VARIANCE_CACHE_SCHEMA_VERSION = "2.0"


def _payload_digest(key: str, mass: np.ndarray, variance: np.ndarray) -> str:
    """Bind the request key to both payload arrays in a portable representation."""
    digest = sha256(f"{VARIANCE_CACHE_SCHEMA_VERSION}:{key}".encode("ascii"))
    for name, value in (("mass", mass), ("variance", variance)):
        digest.update(name.encode("ascii"))
        digest.update(np.asarray(value.shape, dtype="<i8").tobytes())
        digest.update(np.ascontiguousarray(value, dtype="<f8").tobytes())
    return digest.hexdigest()


def variance_cache_key(
    model_identifier: str,
    mass: Any,
    *,
    backend_identifier: str,
    settings: Mapping[str, Any] | None = None,
) -> str:
    """Return a stable SHA-256 key covering model, backend, grid, and settings."""
    for name, value in {
        "model_identifier": model_identifier,
        "backend_identifier": backend_identifier,
    }.items():
        if not isinstance(value, str) or not value.strip():
            raise ValueError(f"{name} must be a non-empty string.")
    mass_array = np.asarray(mass, dtype="<f8")
    if (
        mass_array.ndim != 1
        or mass_array.size < 2
        or not np.all(np.isfinite(mass_array))
        or np.any(mass_array <= 0.0)
        or np.any(np.diff(mass_array) <= 0.0)
    ):
        raise ValueError("Cache mass grid must be finite, positive, and strictly increasing.")
    payload = {
        "schema_version": VARIANCE_CACHE_SCHEMA_VERSION,
        "model_identifier": model_identifier,
        "backend_identifier": backend_identifier,
        "settings": {} if settings is None else dict(settings),
        "mass_shape": mass_array.shape,
        "mass_sha256": sha256(np.ascontiguousarray(mass_array).tobytes()).hexdigest(),
    }
    try:
        encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)
    except (TypeError, ValueError) as exc:
        raise ValueError("Cache settings must be finite JSON-compatible values.") from exc
    return sha256(encoded.encode("utf-8")).hexdigest()


def save_variance_cache(
    path: str | os.PathLike[str],
    *,
    key: str,
    mass: Any,
    variance: Any,
) -> None:
    """Atomically save a variance table with its request key and payload digest.

    The digest detects accidental corruption; it is not an authentication
    mechanism for files from untrusted authors.
    """
    if not isinstance(key, str) or len(key) != 64:
        raise ValueError("key must be a SHA-256 hexadecimal digest.")
    try:
        int(key, 16)
    except ValueError as exc:
        raise ValueError("key must be a SHA-256 hexadecimal digest.") from exc
    mass_array = np.asarray(mass, dtype=float)
    variance_array = np.asarray(variance, dtype=float)
    if mass_array.ndim != 1 or variance_array.shape != mass_array.shape:
        raise ValueError("Mass and variance cache arrays must be aligned one-dimensional data.")
    if mass_array.size < 2:
        raise ValueError("Variance cache grids must contain at least two masses.")
    if (
        not np.all(np.isfinite(mass_array))
        or np.any(mass_array <= 0.0)
        or np.any(np.diff(mass_array) <= 0.0)
        or not np.all(np.isfinite(variance_array))
        or np.any(variance_array < 0.0)
    ):
        raise ValueError("Variance cache arrays contain invalid values.")

    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            dir=destination.parent,
            prefix=f".{destination.name}.",
            suffix=".npz",
            delete=False,
        ) as temporary:
            temporary_name = temporary.name
            np.savez_compressed(
                temporary,
                schema_version=np.asarray(VARIANCE_CACHE_SCHEMA_VERSION),
                key=np.asarray(key),
                payload_sha256=np.asarray(_payload_digest(key, mass_array, variance_array)),
                mass=mass_array,
                variance=variance_array,
            )
        os.replace(temporary_name, destination)
    finally:
        if temporary_name is not None and os.path.exists(temporary_name):
            os.unlink(temporary_name)


def load_variance_cache(
    path: str | os.PathLike[str],
    *,
    expected_key: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Load only when schema, request key, and payload integrity all match.

    Schema 1 caches lack payload verification and must be recomputed.
    """
    with np.load(Path(path), allow_pickle=False) as cache:
        required = {"schema_version", "key", "payload_sha256", "mass", "variance"}
        if set(cache.files) != required:
            raise ValueError("Variance cache has an unsupported field set.")
        if str(cache["schema_version"]) != VARIANCE_CACHE_SCHEMA_VERSION:
            raise ValueError("Variance cache schema version does not match.")
        if str(cache["key"]) != expected_key:
            raise ValueError("Variance cache key does not match the requested model.")
        mass = np.asarray(cache["mass"], dtype=float)
        variance = np.asarray(cache["variance"], dtype=float)
        payload_digest = str(cache["payload_sha256"])
    if (
        mass.ndim != 1
        or mass.size < 2
        or variance.shape != mass.shape
        or not np.all(np.isfinite(mass))
        or np.any(mass <= 0.0)
        or np.any(np.diff(mass) <= 0.0)
        or not np.all(np.isfinite(variance))
        or np.any(variance < 0.0)
    ):
        raise ValueError("Variance cache contains invalid arrays.")
    if _payload_digest(expected_key, mass, variance) != payload_digest:
        raise ValueError("Variance cache payload digest does not match; the cache is corrupt.")
    return mass, variance


__all__ = [
    "VARIANCE_CACHE_SCHEMA_VERSION",
    "load_variance_cache",
    "save_variance_cache",
    "variance_cache_key",
]
