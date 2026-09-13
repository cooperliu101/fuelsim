** Same geometry and labels as b7_rz_contact.e
*Node
1, 0.001, 0
2, 0.002, 0
3, 0.002, 0.001
4, 0.001, 0.001
5, 0.001, 0.0010009999999999999
6, 0.002, 0.0010009999999999999
7, 0.002, 0.0020010000000000002
8, 0.001, 0.0020010000000000002
*Element, type=CAX4T, elset=LOWER
1, 1, 2, 3, 4
*Element, type=CAX4T, elset=UPPER
2, 5, 6, 7, 8
*Nset, nset=BOTTOM
1, 2
*Nset, nset=LOWER_CONTACT
3, 4
*Nset, nset=UPPER_CONTACT
5, 6
*Nset, nset=TOP
7, 8
*Nset, nset=LOWER_RADIAL
1, 2, 3, 4
*Nset, nset=UPPER_RADIAL
5, 6, 7, 8
*Nset, nset=RADIAL_FIXED, generate
1, 8, 1
*Nset, nset=ALL, generate
1, 8, 1
*Surface, type=ELEMENT, name=S_LOWER
LOWER, S3
*Surface, type=ELEMENT, name=S_UPPER
UPPER, S1
