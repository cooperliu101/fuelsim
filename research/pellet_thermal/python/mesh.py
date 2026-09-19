"""Create only fixed research meshes, never production input cards."""

from pathlib import Path
import numpy as np
from scipy.io import netcdf_file

ROOT = Path(__file__).resolve().parents[1] / "cases"


def create(coupled):
    radius, height = 0.004, 0.008
    plane = []
    for y in (-1, 0, 1):
        for x in (-1, 0, 1):
            norm = np.hypot(x, y)
            plane.append((radius * x / norm, radius * y / norm) if norm else (0.0, 0.0))
    points = [(x, y, z) for z in (0.0, height / 2, height) for x, y in plane]
    pellets = []
    for z in range(2):
        for y in range(2):
            for x in range(2):
                a = 9 * z + 3 * y + x
                pellets.append([a, a + 1, a + 4, a + 3, a + 9, a + 10, a + 13, a + 12])
    # Counterclockwise circumference starting at +x; all pellet side facets.
    ring = [5, 8, 7, 6, 3, 0, 1, 2]
    sides = []
    face_edges = [(0, 1), (1, 2), (2, 3), (3, 0)]
    for e, conn in enumerate(pellets):
        for face, (i, j) in enumerate(face_edges):
            if conn[i] % 9 in ring and conn[j] % 9 in ring:
                sides.append((e + 1, face + 1))
    blocks = [("pellet", pellets)]
    surfaces = [
        ("pellet_side", sides),
        ("pellet_boundary", sides + [(i + 1, 5) for i in range(4)] + [(i + 5, 6) for i in range(4)]),
    ]
    if coupled:
        for z in (0.0, height / 2, height):
            for r in (0.0041, 0.0047):
                for i in ring:
                    x, y = plane[i]
                    points.append((x * r / radius, y * r / radius, z))
        clad = []
        for z in range(2):
            for i in range(8):
                a, b = 27 + 16 * z + i, 27 + 16 * z + (i + 1) % 8
                # radial outward then counterclockwise yields positive Jacobian
                clad.append([a, a + 8, b + 8, b, a + 16, a + 24, b + 24, b + 16])
        blocks.append(("cladding", clad))
        surfaces += [
            ("clad_inner", [(9 + i, 4) for i in range(16)]),
            ("clad_outer", [(9 + i, 2) for i in range(16)]),
        ]
    boundary = [i for i in range(27) if i != 13]
    name = "coupled.e" if coupled else "pellet.e"
    with netcdf_file(str(ROOT / name), "w") as f:
        f.api_version = np.float32(7.22)
        f.version = np.float32(7.22)
        f.floating_point_word_size = np.int32(8)
        f.file_size = np.int32(1)
        f.title = "Fixed octagonal pellet thermal research mesh"
        dims = {
            "time_step": None,
            "len_name": 33,
            "num_dim": 3,
            "num_nodes": len(points),
            "num_elem": sum(len(b[1]) for b in blocks),
            "num_el_blk": len(blocks),
            "num_side_sets": len(surfaces),
            "num_node_sets": len(boundary),
        }
        for key, value in dims.items():
            f.createDimension(key, value)

        def var(name, dims, value, typ="i"):
            v = f.createVariable(name, typ, dims)
            v[:] = value
            return v

        def names(name, dim, strings):
            data = np.zeros((len(strings), 33), dtype="S1")
            for i, s in enumerate(strings):
                data[i, : len(s)] = np.frombuffer(s.encode(), dtype="S1")
            var(name, (dim, "len_name"), data, "c")

        for d in range(3):
            var("coord" + "xyz"[d], ("num_nodes",), np.asarray(points)[:, d], "d")
        var("eb_prop1", ("num_el_blk",), np.arange(1, len(blocks) + 1)).name = "ID"
        var("eb_status", ("num_el_blk",), np.ones(len(blocks), dtype=int))
        names("eb_names", "num_el_blk", [b[0] for b in blocks])
        for i, (_, conn) in enumerate(blocks, 1):
            f.createDimension(f"num_el_in_blk{i}", len(conn))
            f.createDimension(f"num_nod_per_el{i}", 8)
            var(
                f"connect{i}", (f"num_el_in_blk{i}", f"num_nod_per_el{i}"), np.asarray(conn) + 1
            ).elem_type = "HEX8"
        var("ss_prop1", ("num_side_sets",), np.arange(1, len(surfaces) + 1)).name = "ID"
        var("ss_status", ("num_side_sets",), np.ones(len(surfaces), dtype=int))
        names("ss_names", "num_side_sets", [s[0] for s in surfaces])
        for i, (_, sides) in enumerate(surfaces, 1):
            dim = f"num_side_ss{i}"
            f.createDimension(dim, len(sides))
            var(f"elem_ss{i}", (dim,), [s[0] for s in sides])
            var(f"side_ss{i}", (dim,), [s[1] for s in sides])
        var("ns_prop1", ("num_node_sets",), np.arange(1, len(boundary) + 1)).name = "ID"
        var("ns_status", ("num_node_sets",), np.ones(len(boundary), dtype=int))
        names("ns_names", "num_node_sets", [f"surface_{i}" for i in boundary])
        for j, i in enumerate(boundary, 1):
            dim = f"num_nod_ns{j}"
            f.createDimension(dim, 1)
            var(f"node_ns{j}", (dim,), [i + 1])
    np.savetxt(
        ROOT / "surface_nodes.txt",
        np.column_stack([boundary, np.asarray(points)[boundary]]),
        fmt=["%d", "%.17g", "%.17g", "%.17g"],
        header="source_internal_id_zero_based x_m y_m z_m",
    )


if __name__ == "__main__":
    ROOT.mkdir(parents=True, exist_ok=True)
    create(False)
    create(True)
