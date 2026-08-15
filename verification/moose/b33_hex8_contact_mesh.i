[Mesh]
  [primary]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 1
    ny = 1
    nz = 1
    xmin = 0
    xmax = 0.01
    ymin = -0.0001
    ymax = 0.0101
    zmin = -0.0001
    zmax = 0.0101
    boundary_name_prefix = primary
  []
  [secondary]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 1
    ny = 1
    nz = 1
    xmin = 0.01
    xmax = 0.02
    ymin = 0
    ymax = 0.01
    zmin = 0
    zmax = 0.01
    boundary_name_prefix = secondary
    boundary_id_offset = 10
  []
  [secondary_id]
    type = RenameBlockGenerator
    input = secondary
    old_block = 0
    new_block = 1
  []
  [combined]
    type = MeshCollectionGenerator
    inputs = 'primary secondary_id'
  []
  [blocks]
    type = RenameBlockGenerator
    input = combined
    old_block = '0 1'
    new_block = 'primary secondary'
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
