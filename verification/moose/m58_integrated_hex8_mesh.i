[Mesh]
  [fuel_cross_section]
    type = ConcentricCircleMeshGenerator
    num_sectors = 4
    radii = '0.00412'
    rings = '4'
    has_outer_square = false
    preserve_volumes = false
  []
  [fuel_extruded]
    type = AdvancedExtruderGenerator
    input = fuel_cross_section
    heights = '0.005'
    num_layers = '8'
    direction = '0 0 1'
    bottom_boundary = fuel_bottom
    top_boundary = fuel_top
  []
  [fuel_named]
    type = RenameBlockGenerator
    input = fuel_extruded
    old_block = 1
    new_block = fuel
  []
  [fuel_outer]
    type = RenameBoundaryGenerator
    input = fuel_named
    old_boundary = outer
    new_boundary = fuel_outer
  []
  [clad_cross_section]
    type = AnnularMeshGenerator
    nr = 2
    nt = 16
    rmin = 0.0041235
    rmax = 0.004692
    quad_subdomain_id = 10
    tri_subdomain_id = 11
    boundary_name_prefix = clad
    boundary_id_offset = 20
  []
  [clad_extruded]
    type = AdvancedExtruderGenerator
    input = clad_cross_section
    heights = '0.007'
    num_layers = '16'
    direction = '0 0 1'
    bottom_boundary = clad_bottom
    top_boundary = clad_top
  []
  [clad_rotated]
    type = ParsedNodeTransformGenerator
    input = clad_extruded
    x_function = '0.999998476913*x - 0.001745328366*y'
    y_function = '0.001745328366*x + 0.999998476913*y'
    z_function = 'z'
  []
  [clad_shifted]
    type = TransformGenerator
    input = clad_rotated
    transform = translate
    vector_value = '0 0 -0.001'
  []
  [clad_named]
    type = RenameBlockGenerator
    input = clad_shifted
    old_block = 10
    new_block = clad
  []

  [combined]
    type = MeshCollectionGenerator
    inputs = 'fuel_outer clad_named'
  []
  [blocks]
    type = RenameBlockGenerator
    input = combined
    old_block = '1 10'
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
