"""Strict articulated-input file admission, independent of the GUI."""

from .articulated import state_vector, tree_links
from .document import MAX_PROJECT_BYTES, Project
from .model import finite_number, read_json


def load_articulated_input(path):
    data = read_json(path, MAX_PROJECT_BYTES)
    if (
        not isinstance(data, dict)
        or set(data) != {"format", "schema", "project", "state"}
        or data["format"] != "vinkulum-articulated-input"
        or type(data["schema"]) is not int
        or data["schema"] != 1
    ):
        raise ValueError("Unknown articulated input format.")
    project = Project.from_dict(data["project"])
    links = tree_links(project)
    state = data["state"]
    if not isinstance(state, dict) or set(state) != {
        "q",
        "velocity",
        "acceleration",
        "effort",
        "time_s",
    }:
        raise ValueError("Incomplete articulated state.")
    for key in ("q", "velocity", "acceleration", "effort"):
        state_vector(state[key], links, key)
    if (
        not finite_number(state["time_s"])
        or not 0 <= state["time_s"] <= project.duration
    ):
        raise ValueError("Invalid load time.")
    return project, state
