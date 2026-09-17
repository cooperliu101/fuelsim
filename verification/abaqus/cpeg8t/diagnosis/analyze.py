"""Identify native integration and contact operators without solving a FE model.

This is a diagnostic, not a replacement acceptance test. Original failed
comparisons and their 0.5 percent gates remain in force.
"""
import hashlib
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT.parent


def verify(directory, manifest):
    for line in (directory / manifest).read_text().splitlines():
        expected, name = line.split(maxsplit=1)
        if hashlib.sha256((directory / name).read_bytes()).hexdigest() != expected:
            raise AssertionError("Reference checksum mismatch: " + name)


def frames(directory, name):
    return json.loads((directory / (name + "_fields.json")).read_text())["steps"]


def material(frame, field, element):
    values = [v for v in frame["fields"][field]["values"] if v["elementLabel"] == element]
    values.sort(key=lambda v: v["integrationPoint"])
    assert [v["integrationPoint"] for v in values] == list(range(1, 10))
    return np.asarray([v["data"] for v in values])


def creep_identity(name, finite):
    history = frames(SOURCE, name)["HISTORY"]
    records = []
    for element in (2, 3):
        old_stress = np.zeros((9, 4))
        old_creep = np.zeros((9, 4))
        old_plastic = np.zeros(9)
        old_time = 0.0
        for frame in history[1:]:
            time = frame["time"]
            stress = material(frame, "S", element)
            creep = material(frame, "CE", element)
            # Abaqus strain output uses engineering shear; constitutive tensors do not.
            creep[:, 3] *= 0.5
            plastic = material(frame, "PEEQ", element) if element == 3 else np.zeros(9)
            dev = stress.copy()
            dev[:, :3] -= np.mean(stress[:, :3], axis=1)[:, None]
            old_dev = old_stress.copy()
            old_dev[:, :3] -= np.mean(old_stress[:, :3], axis=1)[:, None]
            coefficient = 1.5e-7 * (time - old_time)
            forward = coefficient * old_dev
            backward = coefficient * dev
            increment = creep - old_creep
            # These phases follow the prescribed load/hold/reversal history.
            active_plastic = element == 3 and (time <= 1.0 or time >= 2.5)
            implicit = finite or active_plastic
            np.testing.assert_allclose(increment, backward if implicit else forward, rtol=0, atol=1e-12)
            if active_plastic:
                assert np.all(plastic > old_plastic)
            else:
                np.testing.assert_allclose(plastic, old_plastic, rtol=0, atol=1e-12)
            records.append(dict(element=element, time=time,
                                identified_method="backward_euler" if implicit else "forward_euler",
                                forward_maximum_residual=float(np.max(abs(increment - forward))),
                                backward_maximum_residual=float(np.max(abs(increment - backward)))))
            old_time, old_stress, old_creep, old_plastic = time, stress, creep, plastic
    return records


def contact_vectors(frame):
    field = next(v for k, v in frame["fields"].items() if k.startswith("CNORMF"))
    labels = [v["nodeLabel"] for v in field["values"]]
    assert len(set(labels)) == len(labels)
    return {v["nodeLabel"]: np.asarray(v["data"]) for v in field["values"]}


def contact_operator(name, perturbation):
    history = frames(ROOT, name)
    columns = []
    for node in (9, 10, 13):
        forces = []
        for sign in ("POS", "NEG"):
            values = contact_vectors(history[f"N{node}_{sign}"][-1])
            forces.append(np.array([values[n][1] for n in (9, 10, 13)]))
        # Normalize by penalty * initial thickness * edge half-length.
        columns.append(-(forces[0] - forces[1]) / (2 * perturbation * 1e6))
    measured = np.array(columns).T
    consistent = np.array([[4., -1., 2.], [-1., 4., 2.], [2., 2., 16.]]) / 15
    np.testing.assert_allclose(measured.sum(axis=1), [1/3, 1/3, 4/3], rtol=0, atol=1e-6)
    return dict(normalization_N_per_m=1e6, perturbation_m=perturbation,
                native=measured.tolist(), pointwise_integral=consistent.tolist(),
                relative_frobenius_difference=float(np.linalg.norm(measured-consistent)/np.linalg.norm(consistent)),
                maximum_entry_relative_difference=float(np.max(abs((measured-consistent)/consistent))))


def averaged_constraint_reconstruction(operator):
    history = frames(ROOT, "contact_operator_shallow")
    columns = []
    for node in (9, 10, 13):
        gaps = []
        for sign in ("POS", "NEG"):
            frame = history[f"N{node}_{sign}"][-1]
            field = next(v for k, v in frame["fields"].items() if k.startswith("COPEN"))
            values = {v["nodeLabel"]: v["data"] for v in field["values"]}
            gaps.append(np.array([values[n] for n in (9, 10, 13)]))
        columns.append((gaps[0] - gaps[1]) / 2e-10)
    averaging = np.array(columns).T
    # Identify weights from the uniform-force mode, without using the measured
    # nonuniform force derivative matrix. Then predict that matrix independently.
    uniform = contact_vectors(frames(ROOT, "contact_modes")["UNIFORM"][-1])
    uniform_weights = np.array([uniform[n][1] for n in (9, 10, 13)]) / 100.
    weights = np.linalg.solve(averaging.T, uniform_weights)
    assert np.all(weights > 0.)
    predicted = averaging.T @ np.diag(weights) @ averaging
    measured = np.array(operator["native"])
    np.testing.assert_allclose(predicted, measured, rtol=0, atol=1e-7)
    return dict(gap_averaging_matrix=averaging.tolist(), constraint_weights=weights.tolist(),
                reconstructed_force_matrix=predicted.tolist(),
                maximum_absolute_matrix_residual=float(np.max(abs(predicted-measured))),
                scope="Flat single-edge linearization only; these measured coefficients are not production constants.")


def normal_direction():
    secondary_tilt = contact_vectors(frames(ROOT, "contact_modes")["LINEAR"][-1])
    primary_tilt = contact_vectors(frames(ROOT, "contact_primary_tilt")["PRIMARY_TILT"][-1])
    secondary_ratio = np.array([secondary_tilt[n][0] / secondary_tilt[n][1] for n in (9, 10, 13)])
    primary_ratio = np.array([primary_tilt[n][0] / primary_tilt[n][1] for n in (9, 10, 13)])
    np.testing.assert_allclose(secondary_ratio, -.005, rtol=0, atol=1e-12)
    np.testing.assert_allclose(primary_ratio, 0., rtol=0, atol=1e-12)
    return dict(secondary_tilt_native_Fx_over_Fy=secondary_ratio.tolist(),
                secondary_tilt_primary_normal_prediction=0.,
                primary_tilt_native_Fx_over_Fy=primary_ratio.tolist(),
                primary_tilt_primary_normal_prediction=-.005)


def crossing():
    frame = next(f for f in frames(SOURCE, "contact_cycle")["CYCLE"] if f["time"] == 1.)
    reaction = {v["nodeLabel"]: np.asarray(v["data"]) for v in frame["fields"]["RF"]["values"]}
    contact = contact_vectors(frame)
    # The lower body has zero displacement, fixed controls and T=300 K.
    # Zero bulk stress makes its RF a direct contact-force reference.
    for element in (1, 2):
        np.testing.assert_allclose(material(frame, "S", element), 0., rtol=0, atol=1e-7)
    nodes = (4, 7, 3, 13, 10)
    x = np.array([-.02, -.01, 0., .01, .02])
    exact = np.array([-25/6, 125/6, 200/3, 125/6, -25/6])
    native = np.array([reaction[n][1] for n in nodes])
    np.testing.assert_allclose(native.sum(), 100., rtol=0, atol=1e-8)
    np.testing.assert_allclose(native @ x, 0., rtol=0, atol=1e-8)
    return dict(time=1., primary_nodes=nodes, exact_nodal_force_N=exact.tolist(),
                native_RF_N=native.tolist(), native_minus_exact_N=(native-exact).tolist(),
                native_CNORMF_y_N=[float(contact[n][1]) for n in nodes],
                maximum_absolute_difference_N=float(np.max(abs(native-exact))),
                note="Primary CNORMF is not equal to -RF at the end nodes; retain RF and its negative entries.")


if __name__ == "__main__":
    verify(SOURCE, "extended_reference.sha256")
    verify(ROOT, "reference.sha256")
    error = (ROOT / "creep_implicit_parser_error.txt").read_text()
    assert 'Illegal value "IMPLICIT"' in error
    shallow = contact_operator("contact_operator_shallow", 1e-10)
    wide = contact_operator("contact_operator_wide", 1e-10)
    np.testing.assert_allclose(wide["native"], shallow["native"], rtol=0, atol=1e-7)
    report = dict(small_strain_creep=creep_identity("inelastic_history", False),
                  finite_strain_creep=creep_identity("inelastic_finite", True),
                  contact_operator=contact_operator("contact_operator", 1e-7),
                  shallow_contact_operator=shallow,
                  wide_primary_contact_operator=wide,
                  averaged_constraint_reconstruction=averaged_constraint_reconstruction(shallow),
                  contact_normal=normal_direction(), crossing=crossing(),
                  creep_implicit_keyword="rejected by native input processor")
    (ROOT / "diagnosis.json").write_text(json.dumps(report, indent=2) + "\n")
    print("Native creep integration identities, contact operators, normal directions and crossing reactions verified")
