"""GenoCache production: verified native alignment and exact batch caching."""

__version__ = "0.1.0"
__all__ = ["Engine", "EngineError", "ValidationError", "IntegrityError", "ExecutionError", "__version__"]


def __getattr__(name: str):
    # Keep package import light so the deadline-controlled conversion subprocess
    # can execute ``python -m genocache.io`` without importing that module twice.
    if name == "Engine":
        from .engine import Engine
        return Engine
    if name in {"EngineError", "ValidationError", "IntegrityError", "ExecutionError"}:
        from . import io
        return getattr(io, name)
    raise AttributeError(name)
