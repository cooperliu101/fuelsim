from __future__ import print_function

import sys


if len(sys.argv) != 3:
    raise RuntimeError(
        "usage: generate_h20_27_transfer.py <h20_24_source.inp> <transfer_output.inp>"
    )


source_lines = open(sys.argv[1], "rb").read().decode("ascii").splitlines()
coordinates = {}
reading_nodes = False
for line in source_lines:
    if line == "*Node":
        reading_nodes = True
        continue
    if reading_nodes and line.startswith("*"):
        break
    if not reading_nodes:
        continue
    values = [value.strip() for value in line.split(",")]
    label = int(values[0])
    if label > 76:
        continue
    point = [float(values[1]), float(values[2]), float(values[3])]
    if label >= 57:
        point[1] += 0.35
    coordinates[label] = point
if len(coordinates) != 76:
    raise RuntimeError("H20.24 source did not provide the expected first 76 nodes")


output = open(sys.argv[2], "wb")


def write(text):
    output.write(text.encode("ascii"))


write(
    "*Heading\n"
    "** H20.27 Abaqus/Standard small-sliding transfer, release, and recontact identification.\n"
    "** Four primary C3D20 faces oppose one translated secondary C3D20 face.\n"
    "*Preprint, echo=NO, model=NO, history=NO, contact=YES\n"
    "*Node\n"
)
for label in range(1, 77):
    point = coordinates[label]
    write("%d, %.16g, %.16g, %.16g\n" % (label, point[0], point[1], point[2]))
write(
    "*Element, type=C3D20, elset=PRIMARY\n"
    "1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 17, 18, 19, 20, 13, 14, 15, 16\n"
    "2, 4, 3, 21, 22, 8, 7, 23, 24, 11, 25, 26, 27, 19, 30, 31, 32, 16, 15, 28, 29\n"
    "3, 22, 21, 33, 34, 24, 23, 35, 36, 26, 37, 38, 39, 31, 42, 43, 44, 29, 28, 40, 41\n"
    "4, 34, 33, 45, 46, 36, 35, 47, 48, 38, 49, 50, 51, 43, 54, 55, 56, 41, 40, 52, 53\n"
    "*Element, type=C3D20, elset=SECONDARY\n"
    "5, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 73, 74, 75, 76, 69, 70, 71, 72\n"
    "*Elset, elset=ALL\n"
    "PRIMARY, SECONDARY\n"
    "*Nset, nset=PRIMARY_ALL, generate\n"
    "1, 56, 1\n"
    "*Nset, nset=SECONDARY_ALL, generate\n"
    "57, 76, 1\n"
    "*Surface, type=ELEMENT, name=PRIMARY_CONTACT\n"
    "PRIMARY, S4\n"
    "*Surface, type=ELEMENT, name=SECONDARY_CONTACT\n"
    "SECONDARY, S6\n"
    "*Material, name=ELASTIC\n"
    "*Elastic\n"
    "1.e9, 0.\n"
    "*Solid Section, elset=ALL, material=ELASTIC\n"
    ",\n"
    "*Surface Interaction, name=LINEAR_PENALTY\n"
    "*Surface Behavior, penalty=LINEAR\n"
    "1.e8,\n"
    "*Contact Pair, interaction=LINEAR_PENALTY, type=SURFACE TO SURFACE, small sliding, adjust=0.\n"
    "SECONDARY_CONTACT, PRIMARY_CONTACT\n"
    "*Step, name=BASE, nlgeom=NO, inc=40\n"
    "*Static\n"
    "0.1, 1., 1.e-8, 0.1\n"
    "*Boundary\n"
    "PRIMARY_ALL, 1, 3, 0.\n"
    "SECONDARY_ALL, 1, 1, -1.e-4\n"
    "SECONDARY_ALL, 2, 3, 0.\n"
    "*Output, field, frequency=1\n"
    "*Node Output\n"
    "COORD, RF, U\n"
    "*Contact Output\n"
    "CSTRESS, CDISP, CFORCE\n"
    "*Output, history, frequency=1\n"
    "*Contact Output\n"
    "CFN, CMN, CAREA, XN\n"
    "*End Step\n"
    "*Step, name=SHIFT, nlgeom=NO, inc=40\n"
    "*Static\n"
    "0.1, 1., 1.e-8, 0.1\n"
    "*Boundary, op=MOD\n"
    "SECONDARY_ALL, 2, 2, 0.1\n"
    "*End Step\n"
    "*Step, name=OPEN, nlgeom=NO, inc=40\n"
    "*Static\n"
    "0.1, 1., 1.e-8, 0.1\n"
    "*Boundary, op=MOD\n"
    "SECONDARY_ALL, 1, 1, 1.e-4\n"
    "*End Step\n"
    "*Step, name=RECLOSE, nlgeom=NO, inc=40\n"
    "*Static\n"
    "0.1, 1., 1.e-8, 0.1\n"
    "*Boundary, op=MOD\n"
    "SECONDARY_ALL, 1, 1, -1.e-4\n"
    "*End Step\n"
    "*Step, name=RETURN, nlgeom=NO, inc=40\n"
    "*Static\n"
    "0.1, 1., 1.e-8, 0.1\n"
    "*Boundary, op=MOD\n"
    "SECONDARY_ALL, 2, 2, 0.\n"
    "*End Step\n"
)
output.close()
