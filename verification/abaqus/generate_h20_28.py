from __future__ import print_function

import io
import sys


if len(sys.argv) not in (2, 3):
    raise RuntimeError(
        "usage: generate_h20_28.py <h20_28_output.inp> [default|linear|quadratic]"
    )


smoothing = sys.argv[2].lower() if len(sys.argv) == 3 else "default"
if smoothing not in ("default", "linear", "quadratic"):
    raise RuntimeError("unknown sliding-transition smoothing: %s" % smoothing)
sliding_transition = ""
if smoothing != "default":
    sliding_transition = ", sliding transition=%s smoothing" % smoothing.upper()


primary_divisions = 8
nodes = []
primary_node_labels = {}
elements = []


def primary_node(key):
    if key not in primary_node_labels:
        primary_node_labels[key] = len(nodes) + 1
        nodes.append(
            (
                0.5 * key[0],
                float(key[1]) / (2 * primary_divisions),
                float(key[2]) / (2 * primary_divisions),
            )
        )
    return primary_node_labels[key]


def hex20_keys(x0, x1, y0, y1, z0, z1):
    xm = (x0 + x1) // 2
    ym = (y0 + y1) // 2
    zm = (z0 + z1) // 2
    return [
        (x0, y0, z0),
        (x1, y0, z0),
        (x1, y1, z0),
        (x0, y1, z0),
        (x0, y0, z1),
        (x1, y0, z1),
        (x1, y1, z1),
        (x0, y1, z1),
        (xm, y0, z0),
        (x1, ym, z0),
        (xm, y1, z0),
        (x0, ym, z0),
        (xm, y0, z1),
        (x1, ym, z1),
        (xm, y1, z1),
        (x0, ym, z1),
        (x0, y0, zm),
        (x1, y0, zm),
        (x1, y1, zm),
        (x0, y1, zm),
    ]


for y_cell in range(primary_divisions):
    for z_cell in range(primary_divisions):
        keys = hex20_keys(
            0,
            2,
            2 * y_cell,
            2 * y_cell + 2,
            2 * z_cell,
            2 * z_cell + 2,
        )
        elements.append([primary_node(key) for key in keys])

primary_all = list(range(1, len(nodes) + 1))
primary_face = sorted(
    label for key, label in primary_node_labels.items() if key[0] == 2
)

secondary_keys = hex20_keys(2, 4, 0, 2 * primary_divisions, 0, 2 * primary_divisions)
secondary_labels = []
for key in secondary_keys:
    secondary_labels.append(len(nodes) + 1)
    nodes.append(
        (
            0.5 * key[0],
            float(key[1]) / (2 * primary_divisions),
            float(key[2]) / (2 * primary_divisions),
        )
    )
secondary_face_local = [0, 3, 7, 4, 11, 19, 15, 16]
secondary_face = [secondary_labels[index] for index in secondary_face_local]


def write_labels(output, labels):
    for begin in range(0, len(labels), 16):
        output.write(", ".join(str(value) for value in labels[begin : begin + 16]) + "\n")


output = io.open(sys.argv[1], "w", newline="\n")
output.write(
    "*Heading\n"
    "** H20.28 Abaqus/Standard small-sliding C3D20 averaging-region refinement probe.\n"
    "** One secondary face opposes an 8-by-8 primary face mesh.\n"
    "*Preprint, echo=NO, model=NO, history=NO, contact=YES\n"
    "*Node\n"
)
for label, point in enumerate(nodes, 1):
    output.write("%d, %.16g, %.16g, %.16g\n" % (label, point[0], point[1], point[2]))
output.write("*Element, type=C3D20, elset=PRIMARY\n")
for label, connectivity in enumerate(elements, 1):
    output.write("%d, %s\n" % (label, ", ".join(str(value) for value in connectivity)))
secondary_element = len(elements) + 1
output.write("*Element, type=C3D20, elset=SECONDARY\n")
output.write("%d, %s\n" % (secondary_element, ", ".join(str(value) for value in secondary_labels)))
output.write("*Nset, nset=PRIMARY_ALL\n")
write_labels(output, primary_all)
output.write("*Nset, nset=PRIMARY_FACE\n")
write_labels(output, primary_face)
output.write("*Nset, nset=SECONDARY_ALL\n")
write_labels(output, secondary_labels)
output.write("*Nset, nset=SECONDARY_FACE\n")
write_labels(output, secondary_face)
for local_node, label in enumerate(secondary_face, 1):
    output.write("*Nset, nset=S%d\n%d\n" % (local_node, label))
output.write(
    (
        "*Surface, type=ELEMENT, name=PRIMARY_CONTACT\n"
        "PRIMARY, S4\n"
        "*Surface, type=ELEMENT, name=SECONDARY_CONTACT\n"
        "SECONDARY, S6\n"
        "*Material, name=ELASTIC\n"
        "*Elastic\n"
        "1.e9, 0.\n"
        "*Solid Section, elset=PRIMARY, material=ELASTIC\n"
        ",\n"
        "*Solid Section, elset=SECONDARY, material=ELASTIC\n"
        ",\n"
        "*Surface Interaction, name=LINEAR_PENALTY\n"
        "*Surface Behavior, penalty=LINEAR\n"
        "1.e8,\n"
        "*Contact Pair, interaction=LINEAR_PENALTY, type=SURFACE TO SURFACE, small sliding, adjust=0.%s\n"
        "SECONDARY_CONTACT, PRIMARY_CONTACT\n"
        "*Step, name=BASE, nlgeom=NO, inc=40\n"
        "*Static\n"
        "0.1, 1., 1.e-8, 0.1\n"
        "*Boundary\n"
        "PRIMARY_ALL, 1, 3, 0.\n"
        "SECONDARY_ALL, 2, 3, 0.\n"
    )
    % sliding_transition
)
for local_node in range(1, 9):
    output.write("S%d, 1, 1, -1.e-4\n" % local_node)
output.write(
    "*Output, field, frequency=1\n"
    "*Node Output\n"
    "COORD, RF, U\n"
    "*Contact Output\n"
    "CSTRESS, CDISP, CFORCE\n"
    "*End Step\n"
)
for input_node in range(1, 9):
    for suffix, value in (("PLUS", "-1.01e-4"), ("MINUS", "-9.9e-5")):
        output.write(
            "*Step, name=S%d_%s, nlgeom=NO, inc=40\n"
            "*Static\n"
            "0.1, 1., 1.e-8, 0.1\n"
            "*Boundary, op=MOD\n" % (input_node, suffix)
        )
        for local_node in range(1, 9):
            output.write("S%d, 1, 1, %s\n" % (local_node, value if local_node == input_node else "-1.e-4"))
        output.write("*End Step\n")
output.close()
