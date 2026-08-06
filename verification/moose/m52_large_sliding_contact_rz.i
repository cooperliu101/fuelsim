# M5.2 dynamic-search verification. The lower secondary face stretches from
# 4 mm to 4.23 mm while the upper primary face remains 8 mm wide and uses fine
# 0.0625 mm segments. Contact nodes therefore cross multiple primary segments
# without leaving the full chain.

[Mesh]
  [lower_mesh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 4
    ny = 2
    xmin = 0
    xmax = 0.004
    ymin = 0
    ymax = 0.001
    boundary_name_prefix = lower
  []
  [upper_mesh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 128
    ny = 2
    xmin = 0
    xmax = 0.008
    ymin = 0.001002
    ymax = 0.002002
    boundary_name_prefix = upper
    boundary_id_offset = 10
  []
  [upper_id]
    type = RenameBlockGenerator
    input = upper_mesh
    old_block = 0
    new_block = 1
  []
  [combined]
    type = MeshCollectionGenerator
    inputs = 'lower_mesh upper_id'
  []
  [blocks]
    type = RenameBlockGenerator
    input = combined
    old_block = '0 1'
    new_block = 'lower_pellet upper_pellet'
  []
  coord_type = RZ
  patch_update_strategy = iteration
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Variables]
  [T]
    initial_condition = 800
  []
[]

[Physics/SolidMechanics/QuasiStatic]
  [lower]
    block = lower_pellet
    strain = SMALL
    incremental = true
    add_variables = true
  []
  [upper]
    block = upper_pellet
    strain = SMALL
    incremental = true
    add_variables = true
  []
[]

[Kernels]
  [heat]
    type = HeatConduction
    variable = T
    block = 'lower_pellet upper_pellet'
    use_displaced_mesh = false
  []
[]

[Functions]
  [radial_ramp]
    type = PiecewiseLinear
    x = '0 1'
    y = '0 0.00023'
  []
  [upper_radial_ramp]
    type = PiecewiseLinear
    x = '0 1'
    y = '0 0.00004'
  []
  [axial_ramp]
    type = PiecewiseLinear
    x = '0 1'
    y = '0 3e-6'
  []
  [upper_axial_ramp]
    type = PiecewiseLinear
    x = '0 1'
    y = '0 -1e-6'
  []
[]

[Contact]
  [pellet_contact]
    primary = upper_bottom
    secondary = lower_top
    formulation = penalty
    model = frictionless
    penalty = 1e14
    normalize_penalty = true
  []
[]

[BCs]
  [lower_axis]
    type = DirichletBC
    variable = disp_x
    boundary = lower_left
    value = 0
  []
  [lower_outer]
    type = ADFunctionDirichletBC
    variable = disp_x
    boundary = lower_right
    function = radial_ramp
  []
  [upper_axis]
    type = DirichletBC
    variable = disp_x
    boundary = upper_left
    value = 0
  []
  [upper_outer]
    type = ADFunctionDirichletBC
    variable = disp_x
    boundary = upper_right
    function = upper_radial_ramp
  []
  [lower_bottom]
    type = DirichletBC
    variable = disp_y
    boundary = lower_bottom
    value = 0
  []
  [lower_contact_displacement]
    type = ADFunctionDirichletBC
    variable = disp_y
    boundary = lower_top
    function = axial_ramp
  []
  [upper_bottom]
    type = DirichletBC
    variable = disp_y
    boundary = upper_bottom
    value = 0
  []
  [upper_top]
    type = ADFunctionDirichletBC
    variable = disp_y
    boundary = upper_top
    function = upper_axial_ramp
  []
  [lower_temperature]
    type = DirichletBC
    variable = T
    boundary = lower_right
    value = 800
  []
  [upper_temperature]
    type = DirichletBC
    variable = T
    boundary = upper_right
    value = 800
  []
[]

[Materials]
  [conductivity]
    type = GenericConstantMaterial
    block = 'lower_pellet upper_pellet'
    prop_names = thermal_conductivity
    prop_values = 10
  []
  [elasticity]
    type = ComputeIsotropicElasticityTensor
    block = 'lower_pellet upper_pellet'
    youngs_modulus = 2e11
    poissons_ratio = 0.3
  []
  [lower_stress]
    type = ComputeMultipleInelasticStress
    block = lower_pellet
    inelastic_models = ''
  []
  [upper_stress]
    type = ComputeMultipleInelasticStress
    block = upper_pellet
    inelastic_models = ''
  []
[]

[Executioner]
  type = Transient
  end_time = 1
  dt = 0.05
  solve_type = NEWTON
  line_search = basic
  nl_max_its = 80
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-8
  automatic_scaling = true
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[VectorPostprocessors]
  [all_nodes]
    type = NodalValueSampler
    variable = 'T disp_x disp_y'
    sort_by = id
    use_displaced_mesh = false
  []
  [lower_surface]
    type = NodalValueSampler
    boundary = lower_top
    variable = 'T disp_x disp_y contact_pressure nodal_area penetration'
    sort_by = x
    use_displaced_mesh = false
  []
[]

[Outputs]
  csv = true
[]
