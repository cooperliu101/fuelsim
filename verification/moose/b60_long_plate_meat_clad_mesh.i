[Mesh]
  [base]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 50
    ny = 6
    nz = 4
    xmin = 0
    xmax = 0.1
    ymin = 0
    ymax = 0.006
    zmin = 0
    zmax = 0.001
    boundary_name_prefix = plate
    subdomain_name = clad
  []
  [meat]
    type = SubdomainBoundingBoxGenerator
    input = base
    block_id = 1
    block_name = meat
    bottom_left = '0.04 0.001 0.00025'
    top_right = '0.06 0.005 0.00075'
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
