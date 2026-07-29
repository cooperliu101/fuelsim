# M1 verification: RZ thermoelastic fuel and cladding, gap heat transfer,
# and area-normalized node-face penalty contact.

[Mesh]
  [fuel_mesh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 40
    ny = 10
    xmin = 0.0
    xmax = 0.00412
    ymin = 0.0
    ymax = 0.010
    boundary_name_prefix = fuel
  []
  [clad_mesh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 6
    ny = 10
    xmin = 0.004122
    xmax = 0.004692
    ymin = 0.0
    ymax = 0.010020
    boundary_name_prefix = clad
    boundary_id_offset = 10
  []
  [clad_id]
    type = RenameBlockGenerator
    input = clad_mesh
    old_block = '0'
    new_block = '1'
  []
  [combined]
    type = MeshCollectionGenerator
    inputs = 'fuel_mesh clad_id'
  []
  [rename_block]
    type = RenameBlockGenerator
    input = combined
    old_block = '0 1'
    new_block = 'fuel clad'
  []
  coord_type = RZ
  patch_update_strategy = iteration
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Variables]
  [T]
    initial_condition = 600.0
  []
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [fuel]
        block = fuel
        strain = SMALL
        incremental = true
        add_variables = true
        temperature = T
        eigenstrain_names = 'fuel_thermal_strain'
        generate_output = 'vonmises_stress'
      []
      [clad]
        block = clad
        strain = SMALL
        incremental = true
        add_variables = true
        temperature = T
        eigenstrain_names = 'clad_thermal_strain'
        generate_output = 'vonmises_stress'
      []
    []
  []
[]

[Kernels]
  [heat_conduction]
    type = HeatConduction
    variable = T
    block = 'fuel clad'
    use_displaced_mesh = false
  []
  [heat_source_fuel]
    type = BodyForce
    variable = T
    block = fuel
    value = 2.0e8
    function = t
    use_displaced_mesh = false
  []
[]

[ThermalContact]
  [gap]
    type = GapHeatTransfer
    variable = T
    primary = clad_left
    secondary = fuel_right
    gap_conductivity = 0.4
    quadrature = true
    gap_geometry_type = CYLINDER
    min_gap = 1e-6
    min_gap_order = 0
    max_gap = 1e6
    emissivity_primary = 0
    emissivity_secondary = 0
  []
[]

[Contact]
  [mechanical]
    primary = clad_left
    secondary = fuel_right
    formulation = penalty
    model = frictionless
    penalty = 1e14
    normalize_penalty = true
  []
[]

[BCs]
  [u_center]
    type = DirichletBC
    variable = disp_x
    boundary = fuel_left
    value = 0.0
  []
  [u_bottom_f]
    type = DirichletBC
    variable = disp_y
    boundary = fuel_bottom
    value = 0.0
  []
  [u_bottom_c]
    type = DirichletBC
    variable = disp_y
    boundary = clad_bottom
    value = 0.0
  []
  [T_outer]
    type = DirichletBC
    variable = T
    boundary = clad_right
    value = 600.0
  []
[]

[Materials]
  [fuel_tk]
    type = ParsedMaterial
    block = fuel
    property_name = thermal_conductivity
    expression = '3824.0 / T + 0.61'
    coupled_variables = T
  []
  [fuel_elasticity]
    type = ComputeIsotropicElasticityTensor
    block = fuel
    youngs_modulus = 2.0e11
    poissons_ratio = 0.316
  []
  [fuel_expansion]
    type = ComputeThermalExpansionEigenstrain
    block = fuel
    temperature = T
    stress_free_temperature = 600.0
    thermal_expansion_coeff = 10.0e-6
    eigenstrain_name = fuel_thermal_strain
  []
  [fuel_stress]
    type = ComputeMultipleInelasticStress
    block = fuel
    inelastic_models = ''
  []

  [clad_tk]
    type = GenericConstantMaterial
    block = clad
    prop_names = 'thermal_conductivity'
    prop_values = '16.0'
  []
  [clad_elasticity]
    type = ComputeIsotropicElasticityTensor
    block = clad
    youngs_modulus = 75.0e9
    poissons_ratio = 0.3
  []
  [clad_expansion]
    type = ComputeThermalExpansionEigenstrain
    block = clad
    temperature = T
    stress_free_temperature = 600.0
    thermal_expansion_coeff = 5.0e-6
    eigenstrain_name = clad_thermal_strain
  []
  [clad_stress]
    type = ComputeMultipleInelasticStress
    block = clad
    inelastic_models = ''
  []
[]

[Executioner]
  type = Transient
  end_time = 1
  dt = 0.05
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  line_search = 'basic'
  nl_max_its = 50
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-6
  abort_on_solve_fail = true
  automatic_scaling = true
  compute_scaling_once = false
  off_diagonals_in_auto_scaling = true
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Postprocessors]
  [T_center]
    type = PointValue
    point = '0 0.005 0'
    variable = T
    use_displaced_mesh = false
  []
  [T_fuel_surface]
    type = PointValue
    point = '0.00412 0.005 0'
    variable = T
    use_displaced_mesh = false
  []
  [T_clad_inner]
    type = PointValue
    point = '0.004122 0.005 0'
    variable = T
    use_displaced_mesh = false
  []
  [u_fuel_outer]
    type = PointValue
    point = '0.00412 0.005 0'
    variable = disp_x
    use_displaced_mesh = false
  []
  [u_clad_inner]
    type = PointValue
    point = '0.004122 0.005 0'
    variable = disp_x
    use_displaced_mesh = false
  []
  [u_z_fuel_top]
    type = PointValue
    point = '0 0.010 0'
    variable = disp_y
    use_displaced_mesh = false
  []
  [T_mid_ref]
    type = PointValue
    point = '0.00206 0.005 0'
    variable = T
    use_displaced_mesh = false
  []
[]

[VectorPostprocessors]
  [fuel_surface]
    type = NodalValueSampler
    boundary = fuel_right
    variable = 'T disp_x disp_y contact_pressure nodal_area penetration'
    sort_by = y
    use_displaced_mesh = false
  []
  [clad_inner]
    type = NodalValueSampler
    boundary = clad_left
    variable = 'T disp_x disp_y'
    sort_by = y
    use_displaced_mesh = false
  []
[]

[Outputs]
  csv = true
[]
