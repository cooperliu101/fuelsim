[Mesh]
  [primary_a_mesh]
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
    boundary_name_prefix = primary_a
  []
  [secondary_a_mesh]
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
    boundary_name_prefix = secondary_a
    boundary_id_offset = 10
  []
  [primary_b_mesh]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 1
    ny = 1
    nz = 1
    xmin = 0
    xmax = 0.01
    ymin = 0.0199
    ymax = 0.0301
    zmin = -0.0001
    zmax = 0.0101
    boundary_name_prefix = primary_b
    boundary_id_offset = 20
  []
  [secondary_b_mesh]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 1
    ny = 1
    nz = 1
    xmin = 0.01
    xmax = 0.02
    ymin = 0.02
    ymax = 0.03
    zmin = 0
    zmax = 0.01
    boundary_name_prefix = secondary_b
    boundary_id_offset = 30
  []
  [primary_a_block]
    type = RenameBlockGenerator
    input = primary_a_mesh
    old_block = 0
    new_block = 1
  []
  [secondary_a_block]
    type = RenameBlockGenerator
    input = secondary_a_mesh
    old_block = 0
    new_block = 2
  []
  [primary_b_block]
    type = RenameBlockGenerator
    input = primary_b_mesh
    old_block = 0
    new_block = 3
  []
  [secondary_b_block]
    type = RenameBlockGenerator
    input = secondary_b_mesh
    old_block = 0
    new_block = 4
  []
  [combined]
    type = MeshCollectionGenerator
    inputs = 'primary_a_block secondary_a_block primary_b_block secondary_b_block'
  []
  [blocks]
    type = RenameBlockGenerator
    input = combined
    old_block = '1 2 3 4'
    new_block = 'primary_a secondary_a primary_b secondary_b'
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
