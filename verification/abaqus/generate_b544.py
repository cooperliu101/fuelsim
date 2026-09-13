from pathlib import Path


X_NODES = 5
Y_NODES = 2
Z_NODES = 3


def label(x, y, z):
    return z * X_NODES * Y_NODES + y * X_NODES + x + 1


def coordinate(x, y, z):
    x_coordinate = float(x)
    y_coordinate = float(y)
    z_coordinate = 0.5 * float(z)
    if x == 0 or x + 1 == X_NODES:
        return x_coordinate, y_coordinate, z_coordinate
    alternating = -1.0 if x % 2 == 0 else 1.0
    return (
        x_coordinate + 0.06 * (2.0 * float(y) - 1.0) * (float(z) - 1.0),
        y_coordinate + 0.04 * alternating * (float(z) - 1.0),
        z_coordinate + 0.05 * alternating * (float(y) - 0.5),
    )


def generate():
    lines = [
        "*Heading",
        "** B5.44 connected distorted C3D8RT cantilever with temperature-dependent properties.",
        "** The physical coefficients are prescribed independently of the comparison results.",
        "*Preprint, echo=NO, model=NO, history=NO, contact=NO",
        "*Node",
    ]
    for z in range(Z_NODES):
        for y in range(Y_NODES):
            for x in range(X_NODES):
                lines.append("%d, %.16g, %.16g, %.16g" % ((label(x, y, z),) + coordinate(x, y, z)))
    lines.append("*Element, type=C3D8RT, elset=VOLUME")
    element = 0
    for z in range(Z_NODES - 1):
        for x in range(X_NODES - 1):
            element += 1
            connectivity = (
                label(x, 0, z),
                label(x + 1, 0, z),
                label(x + 1, 1, z),
                label(x, 1, z),
                label(x, 0, z + 1),
                label(x + 1, 0, z + 1),
                label(x + 1, 1, z + 1),
                label(x, 1, z + 1),
            )
            lines.append("%d, %s" % (element, ", ".join(str(value) for value in connectivity)))
    left = [label(0, y, z) for z in range(Z_NODES) for y in range(Y_NODES)]
    right = [label(X_NODES - 1, y, z) for z in range(Z_NODES) for y in range(Y_NODES)]
    lines.extend(
        [
            "*Nset, nset=ALL_NODES, generate",
            "1, 30, 1",
            "*Nset, nset=LEFT",
            ", ".join(str(value) for value in left),
            "*Nset, nset=RIGHT",
            ", ".join(str(value) for value in right),
            "*Surface, type=ELEMENT, name=RIGHT_SURFACE",
            "4, S4",
            "8, S4",
            "*Surface, type=ELEMENT, name=TOP_SURFACE",
            "5, S2",
            "6, S2",
            "7, S2",
            "8, S2",
            "*Surface, type=ELEMENT, name=Y_HIGH_SURFACE",
        ]
    )
    lines.extend("%d, S5" % value for value in range(1, 9))
    lines.extend(
        [
            "*Material, name=TEMPERATURE_DEPENDENT",
            "*Elastic",
            "2.04e8, 0.25, 280.0",
            "2.0e8, 0.25, 300.0",
            "1.2e8, 0.25, 700.0",
            "1.0e8, 0.25, 800.0",
            "*Expansion, zero=300.0",
            "1.0e-5",
            "*Conductivity",
            "9.6, 280.0",
            "10.0, 300.0",
            "18.0, 700.0",
            "20.0, 800.0",
            "*Density",
            "2000.0",
            "*Specific Heat",
            "475.0, 280.0",
            "500.0, 300.0",
            "1000.0, 700.0",
            "1125.0, 800.0",
            "*Solid Section, elset=VOLUME, material=TEMPERATURE_DEPENDENT",
            ",",
            "*Hourglass Stiffness",
            "400000.0",
            "*Initial Conditions, type=TEMPERATURE",
            "ALL_NODES, 300.0",
            "*Amplitude, name=CONSTANT_LOAD, definition=TABULAR, time=TOTAL TIME",
            "0.0, 1.0, 1.0e6, 1.0",
            "*Amplitude, name=RIGHT_TEMPERATURE, definition=TABULAR, time=TOTAL TIME",
            "0.0, 0.4285714285714286, 1.0e6, 1.0",
            "*Step, name=PATH, nlgeom=YES, inc=10",
            "*Coupled Temperature-Displacement, deltmx=1000.0",
            "1.0e5, 1.0e6, 1.0e5, 1.0e5",
            "*Boundary",
            "LEFT, 11, 11, 300.0",
            "*Boundary, amplitude=RIGHT_TEMPERATURE",
            "RIGHT, 11, 11, 700.0",
            "*Boundary",
            "LEFT, 1, 3, 0.0",
            "*Dflux, amplitude=CONSTANT_LOAD",
            "VOLUME, BF, 0.5",
            "5, S2, 50.0",
            "6, S2, 50.0",
            "7, S2, 50.0",
            "8, S2, 50.0",
            "*Film, amplitude=CONSTANT_LOAD",
            "1, F5, 280.0, 5.0",
            "2, F5, 280.0, 5.0",
            "3, F5, 280.0, 5.0",
            "4, F5, 280.0, 5.0",
            "5, F5, 280.0, 5.0",
            "6, F5, 280.0, 5.0",
            "7, F5, 280.0, 5.0",
            "8, F5, 280.0, 5.0",
            "*Dsload, follower=NO, amplitude=CONSTANT_LOAD",
            "RIGHT_SURFACE, TRVEC, 2.0e4, 0.0, 0.0, -1.0",
            "*Output, field, frequency=1",
            "*Node Output",
            "COORD, NT, RFL, RF, U",
            "*Element Output, directions=YES",
            "CE, CEEQ, COORD, EE, HFL, IVOL, LE, PE, PEEQ, S, TEMP",
            "*Output, history, frequency=1",
            "*Energy Output",
            "ALLAE, ALLCD, ALLFD, ALLIE, ALLPD, ALLSE, ALLWK",
            "*End Step",
        ]
    )
    output = Path(__file__).with_name("b544_hex8_c3d8rt_distorted_bending.inp")
    output.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    generate()
