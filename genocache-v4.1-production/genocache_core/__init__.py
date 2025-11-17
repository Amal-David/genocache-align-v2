"""
GenoCache V4.1 Production - Core Module

Components:
- encoder: Neural DNA encoder (128D embeddings)
- adaptive_seeding: Seeding with top-k candidate support
- extend_phase: EXTEND phase (THE FIX for chromosome accuracy)
- fast_alignment: Parasail-based alignment
"""

from .encoder import GenoCacheEncoder
from .adaptive_seeding import AdaptiveSeeder
from .extend_phase import ExtendPhase
from .fast_alignment import FastAligner

__version__ = "4.1.0"
__all__ = [
    "GenoCacheEncoder",
    "AdaptiveSeeder",
    "ExtendPhase",
    "FastAligner"
]
