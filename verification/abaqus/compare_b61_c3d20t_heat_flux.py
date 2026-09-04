import csv
import json
import math
import sys


if len(sys.argv) != 5:
    raise RuntimeError(
        "usage: compare_b61_c3d20t_heat_flux.py <mesh.json> <previous_nodal.csv> "
        "<current_nodal.csv> <diagnostic_integration.csv>"
    )


GAUSS3 = (
    (-math.sqrt(3.0 / 5.0), 5.0 / 9.0),
    (0.0, 8.0 / 9.0),
    (math.sqrt(3.0 / 5.0), 5.0 / 9.0),
)
SIGNS = (
    (-1.0, -1.0, -1.0),
    (1.0, -1.0, -1.0),
    (1.0, 1.0, -1.0),
    (-1.0, 1.0, -1.0),
    (-1.0, -1.0, 1.0),
    (1.0, -1.0, 1.0),
    (1.0, 1.0, 1.0),
    (-1.0, 1.0, 1.0),
)


def read_nodal(path):
    with open(path, newline="") as source:
        return {int(row["id"]): row for row in csv.DictReader(source)}


def determinant(matrix):
    return (
        matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
    )


def inverse(matrix):
    value = determinant(matrix)
    return (
        (
            (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) / value,
            (matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2]) / value,
            (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]) / value,
        ),
        (
            (matrix[1][2] * matrix[2][0] - matrix[1][0] * matrix[2][2]) / value,
            (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0]) / value,
            (matrix[0][2] * matrix[1][0] - matrix[0][0] * matrix[1][2]) / value,
        ),
        (
            (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]) / value,
            (matrix[0][1] * matrix[2][0] - matrix[0][0] * matrix[2][1]) / value,
            (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]) / value,
        ),
    )


def displacement_shape_derivatives(xi, eta, zeta):
    derivatives = [[0.0] * 3 for _ in range(20)]
    for node, (sx, sy, sz) in enumerate(SIGNS):
        ax, ay, az = 1.0 + sx * xi, 1.0 + sy * eta, 1.0 + sz * zeta
        total = sx * xi + sy * eta + sz * zeta - 2.0
        derivatives[node] = [
            0.125 * sx * ay * az * (total + ax),
            0.125 * sy * ax * az * (total + ay),
            0.125 * sz * ax * ay * (total + az),
        ]

    def xi_edge(node, sy, sz):
        derivatives[node] = [
            -0.5 * xi * (1.0 + sy * eta) * (1.0 + sz * zeta),
            0.25 * sy * (1.0 - xi * xi) * (1.0 + sz * zeta),
            0.25 * sz * (1.0 - xi * xi) * (1.0 + sy * eta),
        ]

    def eta_edge(node, sx, sz):
        derivatives[node] = [
            0.25 * sx * (1.0 - eta * eta) * (1.0 + sz * zeta),
            -0.5 * eta * (1.0 + sx * xi) * (1.0 + sz * zeta),
            0.25 * sz * (1.0 - eta * eta) * (1.0 + sx * xi),
        ]

    def zeta_edge(node, sx, sy):
        derivatives[node] = [
            0.25 * sx * (1.0 - zeta * zeta) * (1.0 + sy * eta),
            0.25 * sy * (1.0 - zeta * zeta) * (1.0 + sx * xi),
            -0.5 * zeta * (1.0 + sx * xi) * (1.0 + sy * eta),
        ]

    xi_edge(8, -1.0, -1.0)
    eta_edge(9, 1.0, -1.0)
    xi_edge(10, 1.0, -1.0)
    eta_edge(11, -1.0, -1.0)
    zeta_edge(12, -1.0, -1.0)
    zeta_edge(13, 1.0, -1.0)
    zeta_edge(14, 1.0, 1.0)
    zeta_edge(15, -1.0, 1.0)
    xi_edge(16, -1.0, 1.0)
    eta_edge(17, 1.0, 1.0)
    xi_edge(18, 1.0, 1.0)
    eta_edge(19, -1.0, 1.0)
    return derivatives


def temperature_shape_derivatives(xi, eta, zeta):
    return tuple(
        (
            0.125 * sx * (1.0 + sy * eta) * (1.0 + sz * zeta),
            0.125 * sy * (1.0 + sx * xi) * (1.0 + sz * zeta),
            0.125 * sz * (1.0 + sx * xi) * (1.0 + sy * eta),
        )
        for sx, sy, sz in SIGNS
    )


def temperature_shapes(xi, eta, zeta):
    return tuple(
        0.125 * (1.0 + sx * xi) * (1.0 + sy * eta) * (1.0 + sz * zeta)
        for sx, sy, sz in SIGNS
    )


def geometry(coordinates, xi, eta, zeta, quadratic=True):
    derivatives = (
        displacement_shape_derivatives(xi, eta, zeta)
        if quadratic
        else temperature_shape_derivatives(xi, eta, zeta)
    )
    jacobian = [[0.0] * 3 for _ in range(3)]
    for coordinate, derivative in zip(coordinates, derivatives):
        for physical in range(3):
            for natural in range(3):
                jacobian[physical][natural] += coordinate[physical] * derivative[natural]
    return determinant(jacobian), inverse(jacobian)


def metrics(calculated, reference):
    differences = []
    magnitudes = []
    for key in sorted(reference):
        differences.append(
            math.sqrt(sum((left - right) ** 2 for left, right in zip(calculated[key], reference[key])))
        )
        magnitudes.append(math.sqrt(sum(value * value for value in reference[key])))
    return (
        100.0 * math.sqrt(sum(value * value for value in differences) / sum(value * value for value in magnitudes)),
        100.0 * max(differences) / max(magnitudes),
        100.0 * max(difference / magnitude for difference, magnitude in zip(differences, magnitudes) if magnitude),
        max(differences),
    )


with open(sys.argv[1]) as source:
    mesh = json.load(source)
previous = read_nodal(sys.argv[2])
current = read_nodal(sys.argv[3])
with open(sys.argv[4], newline="") as source:
    rows = list(csv.DictReader(source))
required = {"hfl1_w_m2", "hfl2_w_m2", "hfl3_w_m2"}
if not rows or not required.issubset(rows[0]):
    raise RuntimeError("integration file was not written with the diagnostics extraction mode")
reference_flux = {
    (int(row["element"]), int(row["integration_point"])): tuple(
        float(row["hfl%d_w_m2" % component]) for component in range(1, 4)
    )
    for row in rows
}
reference_coordinates = {node["label"]: tuple(node["coordinates"]) for node in mesh["nodes"]}
calculated_flux = {name: {} for name in ("reference", "increment_midpoint", "current")}
corner_nodes = set(node for element in mesh["elements"] for node in element["nodes"][:8])
fixed_temperature_nodes = set()
for node_set in mesh["node_sets"]:
    if node_set["name"].lower() in ("plate_back", "plate_front"):
        fixed_temperature_nodes.update(node_set["nodes"])
thermal_residual = {node: 0.0 for node in corner_nodes}

for element in mesh["elements"]:
    labels = element["nodes"]
    configurations = {
        "reference": [reference_coordinates[label] for label in labels],
        "increment_midpoint": [
            tuple(
                reference_coordinates[label][component]
                + 0.5
                * (
                    float(previous[label]["displacement_" + "xyz"[component]])
                    + float(current[label]["displacement_" + "xyz"[component]])
                )
                for component in range(3)
            )
            for label in labels
        ],
        "current": [
            tuple(
                reference_coordinates[label][component]
                + float(current[label]["displacement_" + "xyz"[component]])
                for component in range(3)
            )
            for label in labels
        ],
    }
    temperatures = [float(current[label]["temperature"]) for label in labels[:8]]
    previous_temperatures = [float(previous[label]["temperature"]) for label in labels[:8]]
    if element["block"] == "MEAT":
        conductivity, capacity, source = 3.0, 10970.0 * 300.0, 2.0e8
    else:
        conductivity, capacity, source = 16.0, 6500.0 * 330.0, 0.0
    point = 0
    for zeta, weight_z in GAUSS3:
        for eta, weight_y in GAUSS3:
            for xi, weight_x in GAUSS3:
                point += 1
                weight = weight_x * weight_y * weight_z
                shapes = temperature_shapes(xi, eta, zeta)
                natural_derivatives = temperature_shape_derivatives(xi, eta, zeta)
                for name, coordinates in configurations.items():
                    mapping = geometry(coordinates, xi, eta, zeta)[1]
                    gradient = tuple(
                        sum(
                            natural_derivatives[node][natural] * mapping[natural][physical] * temperatures[node]
                            for node in range(8)
                            for natural in range(3)
                        )
                        for physical in range(3)
                    )
                    calculated_flux[name][(element["label"], point)] = tuple(
                        -conductivity * component for component in gradient
                    )
                current_determinant = geometry(configurations["current"], xi, eta, zeta)[0]
                midpoint_mapping = geometry(configurations["increment_midpoint"], xi, eta, zeta)[1]
                linear_current_determinant = geometry(
                    configurations["current"][:8], xi, eta, zeta, quadratic=False
                )[0]
                midpoint_gradients = [
                    tuple(
                        sum(
                            natural_derivatives[node][natural] * midpoint_mapping[natural][physical]
                            for natural in range(3)
                        )
                        for physical in range(3)
                    )
                    for node in range(8)
                ]
                midpoint_temperature_gradient = tuple(
                    sum(midpoint_gradients[node][physical] * temperatures[node] for node in range(8))
                    for physical in range(3)
                )
                temperature_rate = sum(
                    shapes[node] * (temperatures[node] - previous_temperatures[node]) / 2.0
                    for node in range(8)
                )
                for node, label in enumerate(labels[:8]):
                    conduction = sum(
                        midpoint_gradients[node][physical] * midpoint_temperature_gradient[physical]
                        for physical in range(3)
                    )
                    thermal_residual[label] += weight * (
                        conductivity * current_determinant * conduction
                        + capacity * current_determinant * shapes[node] * temperature_rate
                        - source * linear_current_determinant * shapes[node]
                    )

if set(reference_flux) != set(calculated_flux["reference"]):
    raise RuntimeError("mesh and integration files do not contain the same material points")
for name in ("reference", "increment_midpoint", "current"):
    values = metrics(calculated_flux[name], reference_flux)
    print(
        "%s_relative_l2_percent=%.12g\n"
        "%s_relative_peak_percent=%.12g\n"
        "%s_maximum_pointwise_percent=%.12g\n"
        "%s_maximum_absolute_difference_w_m2=%.12g"
        % ((name, values[0], name, values[1], name, values[2], name, values[3]))
    )
free_residual = [
    value for node, value in thermal_residual.items() if node not in fixed_temperature_nodes
]
print(
    "abaqus_operator_free_thermal_residual_l2_w=%.12g\n"
    "abaqus_operator_maximum_free_thermal_residual_w=%.12g"
    % (
        math.sqrt(sum(value * value for value in free_residual)),
        max(abs(value) for value in free_residual),
    )
)
