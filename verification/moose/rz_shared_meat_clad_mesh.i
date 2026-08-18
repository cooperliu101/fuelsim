[Mesh]
  [base]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 2
    ny = 1
    xmin = 1
    xmax = 3
    ymin = 0
    ymax = 1
    boundary_name_prefix = plate
    subdomain_name = meat
  []
  [clad]
    type = SubdomainBoundingBoxGenerator
    input = base
    block_id = 1
    block_name = clad
    bottom_left = '2 0 0'
    top_right = '3 1 0'
  []
  coord_type = RZ
[]

[Problem]
  solve = false
[]

[Executioner]
  type = Steady
[]

[Outputs]
  exodus = true
[]
