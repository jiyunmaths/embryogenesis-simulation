"""Opt-in native backend retaining the fertilization assay's cue and clamp hooks."""
from .fertilization_cue import CueSimulation
from .native_mechanics import NativeSimulation


class NativeCueSimulation(CueSimulation, NativeSimulation):
    """Cue hooks wrap native steps/mechanics; original checkpoints stay readable."""
