[Mesh]
  [fuel]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 2
    ny = 4
    nz = 1
    xmin = 0
    xmax = 0.00412
    ymin = 0
    ymax = 0.005
    zmin = 0
    zmax = 0.001
    boundary_name_prefix = fuel
  []
  [clad]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 1
    ny = 16
    nz = 1
    xmin = 0.00412
    xmax = 0.004692
    ymin = -0.001
    ymax = 0.006
    zmin = -0.0002
    zmax = 0.0012
    boundary_name_prefix = clad
    boundary_id_offset = 10
  []
  [clad_id]
    type = RenameBlockGenerator
    input = clad
    old_block = 0
    new_block = 1
  []
  [combined]
    type = MeshCollectionGenerator
    inputs = 'fuel clad_id'
  []
  [blocks]
    type = RenameBlockGenerator
    input = combined
    old_block = '0 1'
    new_block = 'fuel clad'
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
