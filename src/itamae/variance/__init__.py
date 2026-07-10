"""Composable variance adapters and future power-spectrum implementations.

The initial module provides a callable adapter so existing SASHIMI variance
implementations can satisfy the common protocol before their physical formulae
are moved into dedicated ITAMAE components.
"""

from .callable import CallableVarianceModel

__all__ = ["CallableVarianceModel"]
