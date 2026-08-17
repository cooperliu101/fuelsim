[Mesh]
  [base]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 2
    ny = 1
    nz = 1
    xmin = 0
    xmax = 2
    ymin = 0
    ymax = 1
    zmin = 0
    zmax = 1
    boundary_name_prefix = plate
    subdomain_name = meat
  []
  [clad]
    type = SubdomainBoundingBoxGenerator
    input = base
    block_id = 1
    block_name = clad
    bottom_left = '1 0 0'
    top_right = '2 1 1'
  []
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
