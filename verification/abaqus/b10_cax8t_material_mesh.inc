** CAX8T material mesh matching b10_cax8t_material.e
*Node
1, 0, 0
2, 0.001, 0
3, 0.001, 0.001
4, 0, 0.001
5, 0.00050000000000000001, 0
6, 0.001, 0.00050000000000000001
7, 0.00050000000000000001, 0.001
8, 0, 0.00050000000000000001
*Element, type=CAX8T, elset=SOLID
1, 1, 2, 3, 4, 5, 6, 7, 8
*Nset, nset=ALL
1, 2, 3, 4, 5, 6, 7, 8
*Nset, nset=LEFT
1, 4, 8
*Nset, nset=BOTTOM
1, 2, 5
*Nset, nset=TOP
3, 4, 7
*Nset, nset=AXIS
1, 4, 8
