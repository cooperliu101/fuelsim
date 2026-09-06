"""Write only the geometry for the CAX8T contact recovery experiment."""
from generate_b10_cax8t_meshes import ROOT, MESH_ROOT, enrich_quad4_mesh, write_exodus, append_nset


coordinates, blocks = enrich_quad4_mesh(
    [(r, z) for z in (0.001, 0.002, 0.0035) for r in (0.001, 0.002)] +
    [(0.002001, 0), (0.003001, 0), (0.003001, 0.004), (0.002001, 0.004)],
    [("inner", [[1, 2, 4, 3], [3, 4, 6, 5]]), ("outer", [[7, 8, 9, 10]])])
edges = [blocks[0][1][0], blocks[0][1][1]]
contact_nodes = [2, 4, 6, edges[0][5], edges[1][5]]
node_sets = [("all", list(range(1, len(coordinates) + 1))),
             ("fixed", [n for n in range(1, len(coordinates) + 1) if n not in contact_nodes])]
node_sets += [("c%d" % i, [n]) for i, n in enumerate(contact_nodes)]
boundaries = [("inner_right", [1, 2], [2, 2]), ("outer_left", [3], [4])]
write_exodus(MESH_ROOT / "b114_cax8t_recovery.e", "CAX8T two-edge contact recovery",
             coordinates, blocks, boundaries, node_sets)
lines = ["** Geometry for the two-edge CAX8T recovery experiment", "*Node"]
lines += ["%d, %.17g, %.17g" % (i, r, z) for i, (r, z) in enumerate(coordinates, 1)]
label = 1
for name, elements in blocks:
    lines.append("*Element, type=CAX8T, elset=" + name.upper())
    for element in elements:
        lines.append(", ".join(map(str, [label] + element)))
        label += 1
for name, nodes in node_sets:
    append_nset(lines, name, nodes)
lines += ["*Surface, type=ELEMENT, name=S_INNER", "INNER, S2",
          "*Surface, type=ELEMENT, name=S_OUTER", "OUTER, S4"]
(ROOT / "b10_cax8t_recovery_mesh.inc").write_text("\n".join(lines) + "\n")
