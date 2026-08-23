from __future__ import print_function

import sys


if len(sys.argv) != 3:
    raise RuntimeError(
        "usage: generate_h20_27.py <h20_24_source.inp> <h20_27_output.inp>"
    )


secondary_labels = [57, 60, 64, 61, 68, 72, 76, 69, 78, 80, 83, 85, 88]


source = open(sys.argv[1], "rb").read()
marker = b"*Nset, nset=PRIMARY_X0"
position = source.find(marker)
if position < 0:
    raise RuntimeError("H20.24 source input does not contain the expected set marker")

output = open(sys.argv[2], "wb")
def write(text):
    output.write(text.encode("ascii"))


write(
    "*Heading\n"
    "** H20.27 Abaqus/Standard nonmatching C3D20 surface-to-surface operator identification.\n"
    "** Geometry and connectivity are copied from the tracked H20.24 input.\n"
    "*Preprint, echo=NO, model=NO, history=NO, contact=YES\n"
)
source_body = source[source.find(b"*Node"):position]
output.write(source_body)
write("*Nset, nset=PRIMARY_ALL, generate\n1, 56, 1\n")
write("*Nset, nset=SECONDARY_ALL, generate\n57, 88, 1\n")
for local_node, label in enumerate(secondary_labels, 1):
    write("*Nset, nset=S%02d\n%d\n" % (local_node, label))
write(
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
    "SECONDARY_ALL, 2, 3, 0.\n"
)
for local_node in range(1, len(secondary_labels) + 1):
    write("S%02d, 1, 1, -1.e-4\n" % local_node)
write(
    "*Output, field, frequency=1\n"
    "*Node Output\n"
    "COORD, RF, U\n"
    "*Contact Output\n"
    "CSTRESS, CDISP, CFORCE\n"
    "*Output, history, frequency=1\n"
    "*Contact Output\n"
    "CFN, CMN, CAREA, XN\n"
    "*End Step\n"
)

for input_node in range(1, len(secondary_labels) + 1):
    for suffix, perturbed_value in (("PLUS", "-1.01e-4"), ("MINUS", "-9.9e-5")):
        write(
            "*Step, name=S%02d_%s, nlgeom=NO, inc=40\n"
            "*Static\n"
            "0.1, 1., 1.e-8, 0.1\n"
            "*Boundary, op=MOD\n" % (input_node, suffix)
        )
        for local_node in range(1, len(secondary_labels) + 1):
            value = perturbed_value if local_node == input_node else "-1.e-4"
            write("S%02d, 1, 1, %s\n" % (local_node, value))
        write("*End Step\n")

for step_name, value in (("OPEN", "1.e-4"), ("RECLOSE", "-1.e-4")):
    write(
        "*Step, name=%s, nlgeom=NO, inc=40\n"
        "*Static\n"
        "0.1, 1., 1.e-8, 0.1\n"
        "*Boundary, op=MOD\n" % step_name
    )
    for local_node in range(1, len(secondary_labels) + 1):
        write("S%02d, 1, 1, %s\n" % (local_node, value))
    write("*End Step\n")

output.close()
