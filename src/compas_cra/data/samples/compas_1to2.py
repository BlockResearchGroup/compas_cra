import json
import pathlib

old_new = {
    "value": "data",
    "dna": "default_node_attributes",
    "dea": "default_edge_attributes",
    "dva": "default_vertex_attributes",
    "dfa": "default_face_attributes",
    "dca": "default_cell_attributes",
}


def replace_dict_keys(olddata, newdata):
    for key, value in olddata.items():
        key = old_new.get(key, key)
        if isinstance(value, dict):
            newdata[key] = {}
            replace_dict_keys(value, newdata[key])
        else:
            newdata[key] = value


def replace_untyped_graph_values(olddata, newdata):
    for key, value in olddata.items():
        if key == "graph":
            if "dtype" not in value:
                newdata[key] = {
                    "dtype": "compas.datastructures/Graph",
                    "data": value,
                }
            else:
                newdata[key] = value
        else:
            if isinstance(value, dict):
                newdata[key] = {}
                replace_untyped_graph_values(value, newdata[key])
            else:
                newdata[key] = value


def clean_up_interface_data(olddata):
    if isinstance(olddata, dict):
        for key, value in olddata.items():
            if isinstance(value, dict):
                if "type" in value and "interaction" in value and "viewmesh" in value:
                    print("here")
                    del value["type"]
                    del value["interaction"]
                    del value["viewmesh"]
            clean_up_interface_data(value)
    elif isinstance(olddata, list):
        for item in olddata:
            clean_up_interface_data(item)


here = pathlib.Path(__file__).parent

for filepath in here.iterdir():
    if filepath.suffix != ".json":
        continue

    with open(filepath) as f:
        olddata = json.load(f)
        newdata = {}
        replace_dict_keys(olddata, newdata)

        olddata = newdata
        newdata = {}
        replace_untyped_graph_values(olddata, newdata)

        clean_up_interface_data(newdata)

    with open(here / f"{filepath.stem}.json", "w+") as f:
        json.dump(newdata, f)
