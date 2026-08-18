[Mesh]
  [base]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 4
    ny = 4
    nz = 4
    xmin = 0
    xmax = 0.004
    ymin = 0
    ymax = 0.004
    zmin = 0
    zmax = 0.004
    boundary_name_prefix = plate
    subdomain_name = clad
  []
  [meat]
    type = SubdomainBoundingBoxGenerator
    input = base
    block_id = 1
    block_name = meat
    bottom_left = '0.001 0.001 0.001'
    top_right = '0.003 0.003 0.003'
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
