# M3.3 verification: two independently meshed RZ pellets with axial thermal
# and mechanical contact.  The lower/upper contact faces deliberately use
# four and six Line2 segments.

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
    nx = 6
    ny = 2
    xmin = 0
    xmax = 0.0041
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

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [lower]
        block = lower_pellet
        strain = SMALL
        incremental = true
        add_variables = true
        temperature = T
        eigenstrain_names = lower_thermal_strain
      []
      [upper]
        block = upper_pellet
        strain = SMALL
        incremental = true
        add_variables = true
        temperature = T
        eigenstrain_names = upper_thermal_strain
      []
    []
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

[ThermalContact]
  [pellet_gap]
    type = GapHeatTransfer
    variable = T
    primary = upper_bottom
    secondary = lower_top
    gap_conductivity = 0.2
    quadrature = true
    gap_geometry_type = PLATE
    min_gap = 1e-6
    min_gap_order = 0
    max_gap = 1e6
    emissivity_primary = 0
    emissivity_secondary = 0
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
  [upper_axis]
    type = DirichletBC
    variable = disp_x
    boundary = upper_left
    value = 0
  []
  [lower_bottom]
    type = DirichletBC
    variable = disp_y
    boundary = lower_bottom
    value = 0
  []
  [upper_top]
    type = DirichletBC
    variable = disp_y
    boundary = upper_top
    value = 0
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
  [lower_expansion]
    type = ComputeThermalExpansionEigenstrain
    block = lower_pellet
    temperature = T
    stress_free_temperature = 600
    thermal_expansion_coeff = 1e-5
    eigenstrain_name = lower_thermal_strain
  []
  [upper_expansion]
    type = ComputeThermalExpansionEigenstrain
    block = upper_pellet
    temperature = T
    stress_free_temperature = 600
    thermal_expansion_coeff = 1e-5
    eigenstrain_name = upper_thermal_strain
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
  dt = 1
  solve_type = NEWTON
  line_search = basic
  nl_max_its = 50
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
