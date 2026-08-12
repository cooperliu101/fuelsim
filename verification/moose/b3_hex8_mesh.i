[Mesh]
  [mesh]
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
    boundary_name_prefix = solid
    subdomain_name = solid
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
