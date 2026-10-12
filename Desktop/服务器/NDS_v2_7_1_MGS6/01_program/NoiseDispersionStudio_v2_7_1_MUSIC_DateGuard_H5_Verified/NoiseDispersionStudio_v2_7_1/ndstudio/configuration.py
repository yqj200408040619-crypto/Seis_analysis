"""Preserve current on-disk settings while applying explicit GUI edits."""
from copy import deepcopy


def apply_editor_changes(base,current,edited):
    """Three-way merge: untouched editor values cannot undo worker updates."""
    result=deepcopy(current)
    for key,value in edited.items():
        if isinstance(value,dict) and isinstance(base.get(key),dict) and isinstance(current.get(key),dict):
            result[key]=apply_editor_changes(base[key],current[key],value)
        elif key not in base or value!=base[key]:result[key]=deepcopy(value)
    return result
