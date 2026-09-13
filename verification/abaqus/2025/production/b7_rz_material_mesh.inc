** 1 mm radius, 1 mm height. Same nodes and connectivity as b7_rz_material.e.
*Node
1, 0, 0
2, 0.001, 0
3, 0.001, 0.001
4, 0, 0.001
*Element, type=CAX4T, elset=SOLID
1, 1, 2, 3, 4
*Nset, nset=ALL
1, 2, 3, 4
*Nset, nset=AXIS
1, 4
*Nset, nset=BOTTOM
1, 2
*Nset, nset=TOP
3, 4
