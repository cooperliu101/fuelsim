"""Geometry only: preserve the three original PCMI meshes and enrich their edges."""
from scipy.io import netcdf_file
from generate_b10_cax8t_meshes import ROOT, MESH_ROOT, enrich_quad4_mesh, write_exodus, append_nset


def names(variable):
    return [b"".join(row).decode("ascii").rstrip("\0 ") for row in variable.data]


def write_mesh(group, source):
    with netcdf_file(str(ROOT.parent / "moose" / source), "r", mmap=False) as file:
        coordinates = list(zip(file.variables["coordx"].data.tolist(), file.variables["coordy"].data.tolist()))
        blocks = [(name, file.variables["connect%d" % i].data.tolist())
                  for i, name in enumerate(names(file.variables["eb_names"]), 1)]
        boundaries = [(name, file.variables["elem_ss%d" % i].data.tolist(),
                       file.variables["side_ss%d" % i].data.tolist())
                      for i, name in enumerate(names(file.variables["ss_names"]), 1)]
    for element in ("cax4t", "cax4rt", "cax8t", "cax8rt"):
        xy, topology = enrich_quad4_mesh(coordinates, blocks) if "8" in element else (coordinates, blocks)
        flat = [e for _, elements in topology for e in elements]
        node_sets = [("all", list(range(1, len(xy) + 1)))]
        for name, elements, sides in boundaries:
            nodes = set()
            for label, side in zip(elements, sides):
                e, s = flat[label - 1], side - 1
                nodes.update((e[s], e[(s + 1) % 4]))
                if len(e) == 8:
                    nodes.add(e[4 + s])
            node_sets.append((name, sorted(nodes)))
        prefix = "b13_%s_%s" % (group, element)
        if "8" in element:
            write_exodus(MESH_ROOT / (prefix + ".e"), prefix, xy, topology, boundaries, node_sets)
        lines = ["** Original PCMI geometry: " + source, "*Node"]
        lines += ["%d, %.17g, %.17g" % (n, r, z) for n, (r, z) in enumerate(xy, 1)]
        label = 0
        for name, elements in topology:
            lines.append("*Element, type=%s, elset=%s" % (element.upper(), name.upper()))
            for e in elements:
                label += 1
                lines.append(", ".join(map(str, [label] + e)))
        for name, nodes in node_sets:
            append_nset(lines, name, nodes)
        for name, elements, sides in boundaries:
            lines.append("*Surface, type=ELEMENT, name=S_" + name.upper())
            lines += ["%d, S%d" % pair for pair in zip(elements, sides)]
        (ROOT / (prefix + "_mesh.inc")).write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    write_mesh("small", "m23_pcmi_coupled_cladding_rz_mesh.e")
    write_mesh("finite", "m41_finite_strain_pcmi_rz_mesh.e")
    write_mesh("integrated", "m57_integrated_fuel_cladding_rz_mesh.e")
