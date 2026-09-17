"""Extract native CPEG8T fields without reconstructing stresses or reactions."""
import json
import os
import sys

from abaqusConstants import DOUBLE_PRECISION
from odbAccess import openOdb


def data(value):
    result = value.dataDouble if value.precision == DOUBLE_PRECISION else value.data
    try:
        return [float(component) for component in result]
    except TypeError:
        return float(result)


path = os.path.abspath(sys.argv[1])
odb = openOdb(path, readOnly=True)
result = {"steps": {}, "history": {}}
for name, step in odb.steps.items():
    frames = []
    for frame in step.frames:
        fields = {}
        for field_name, field in frame.fieldOutputs.items():
            values = []
            for value in field.values:
                row = {"data": data(value)}
                for key in ("nodeLabel", "elementLabel", "integrationPoint"):
                    identifier = getattr(value, key, None)
                    if identifier is not None:
                        row[key] = identifier
                values.append(row)
            fields[field_name] = {"components": list(field.componentLabels), "values": values}
        frames.append({"time": frame.frameValue, "fields": fields})
    result["steps"][name] = frames
    result["history"][name] = {
        region_name: {key: list(output.data) for key, output in region.historyOutputs.items()}
        for region_name, region in step.historyRegions.items()
    }
odb.close()
with open(os.path.splitext(path)[0] + "_fields.json", "w") as stream:
    json.dump(result, stream, indent=2, sort_keys=True)
