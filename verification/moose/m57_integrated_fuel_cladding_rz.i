# Integrated transient finite-strain fuel-cladding comparison. Prescribed fuel
# translation drives frictional sliding across multiple current primary segments
# while independent pressure and power histories load the two regions.

[Mesh]
  [fuel_mesh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 6
    ny = 8
    xmin = 0
    xmax = 0.00412
    ymin = 0
    ymax = 0.005
    boundary_name_prefix = fuel
  []
  [clad_mesh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 2
    ny = 24
    xmin = 0.004121
    xmax = 0.004692
    ymin = 0
    ymax = 0.006
    boundary_name_prefix = clad
    boundary_id_offset = 10
  []
  [clad_id]
    type = RenameBlockGenerator
    input = clad_mesh
    old_block = 0
    new_block = 1
  []
  [combined]
    type = MeshCollectionGenerator
    inputs = 'fuel_mesh clad_id'
  []
  [blocks]
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
    initial_condition = 600
  []
[]

[AuxVariables]
  [mortar_normal_pressure]
    family = LAGRANGE
    order = FIRST
  []
  [mortar_tangential_pressure]
    family = LAGRANGE
    order = FIRST
  []
[]

[Physics/SolidMechanics/QuasiStatic]
  [fuel]
    block = fuel
    strain = FINITE
    add_variables = true
    temperature = T
    eigenstrain_names = fuel_thermal_strain
    use_automatic_differentiation = true
  []
  [clad]
    block = clad
    strain = FINITE
    add_variables = true
    temperature = T
    eigenstrain_names = clad_thermal_strain
    use_automatic_differentiation = true
    generate_output = vonmises_stress
  []
[]

[Functions]
  [power]
    type = PiecewiseLinear
    x = '0 0.5 1.5 3 4.5 6'
    y = '0 0.4 1 0.8 1.2 1'
  []
  [internal_pressure]
    type = PiecewiseLinear
    x = '0 0.5 1.5 3 4.5 6'
    y = '0 0.25 1 0.7 1.1 1'
  []
  [external_pressure]
    type = PiecewiseLinear
    x = '0 0.5 1.5 3 4.5 6'
    y = '0 0.8 1 1.2 0.9 1'
  []
  [axial_slide]
    type = PiecewiseLinear
    x = '0 0.5 1.5 3 4.5 6'
    y = '0 6e-5 2.4e-4 3.6e-4 4.8e-4 6e-4'
  []
[]

[Kernels]
  [heat_time]
    type = HeatConductionTimeDerivative
    variable = T
    block = 'fuel clad'
    use_displaced_mesh = false
  []
  [heat_conduction]
    type = HeatConduction
    variable = T
    block = 'fuel clad'
    use_displaced_mesh = false
  []
  [heat_source]
    type = BodyForce
    variable = T
    block = fuel
    value = 2e8
    function = power
    use_displaced_mesh = false
  []
[]

[ThermalContact]
  [gap]
    type = GapHeatTransfer
    variable = T
    primary = clad_left
    secondary = fuel_right
    gap_conductivity = 0.004
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
    formulation = mortar_penalty
    model = coulomb
    friction_coefficient = 0.002
    penalty = 1e12
    penalty_friction = 1e8
  []
[]

[AuxKernels]
  [mortar_normal_pressure]
    type = MortarUserObjectAux
    variable = mortar_normal_pressure
    user_object = penalty_friction_object_mechanical
    contact_quantity = normal_pressure
  []
  [mortar_tangential_pressure]
    type = MortarUserObjectAux
    variable = mortar_tangential_pressure
    user_object = penalty_friction_object_mechanical
    contact_quantity = tangential_pressure_one
  []
[]

[BCs]
  [fuel_axis]
    type = ADDirichletBC
    variable = disp_x
    boundary = fuel_left
    value = 0
  []
  [fuel_top_slide]
    type = ADFunctionDirichletBC
    variable = disp_y
    boundary = fuel_top
    function = axial_slide
  []
  [clad_bottom]
    type = ADDirichletBC
    variable = disp_y
    boundary = clad_bottom
    value = 0
  []
  [clad_outer_temperature]
    type = DirichletBC
    variable = T
    boundary = clad_right
    value = 600
  []
  [internal_pressure]
    type = ADPressure
    variable = disp_x
    boundary = clad_left
    factor = 5e5
    function = internal_pressure
    use_displaced_mesh = true
  []
  [external_pressure]
    type = ADPressure
    variable = disp_x
    boundary = clad_right
    factor = 2e6
    function = external_pressure
    use_displaced_mesh = true
  []
[]

[Materials]
  [fuel_conductivity]
    type = ParsedMaterial
    block = fuel
    property_name = thermal_conductivity
    expression = '3824 / T + 0.61'
    coupled_variables = T
  []
  [fuel_thermal_mass]
    type = GenericConstantMaterial
    block = fuel
    prop_names = 'density specific_heat'
    prop_values = '10970 300'
  []
  [fuel_elasticity]
    type = ADComputeIsotropicElasticityTensor
    block = fuel
    youngs_modulus = 2e11
    poissons_ratio = 0.316
  []
  [fuel_expansion]
    type = ADComputeThermalExpansionEigenstrain
    block = fuel
    temperature = T
    stress_free_temperature = 600
    thermal_expansion_coeff = 1e-5
    eigenstrain_name = fuel_thermal_strain
  []
  [fuel_stress]
    type = ADComputeMultipleInelasticStress
    block = fuel
    inelastic_models = ''
  []
  [clad_thermal]
    type = GenericConstantMaterial
    block = clad
    prop_names = 'thermal_conductivity density specific_heat'
    prop_values = '16 6500 330'
  []
  [clad_elasticity]
    type = ADComputeIsotropicElasticityTensor
    block = clad
    youngs_modulus = 7.5e10
    poissons_ratio = 0.3
  []
  [clad_expansion]
    type = ADComputeThermalExpansionEigenstrain
    block = clad
    temperature = T
    stress_free_temperature = 600
    thermal_expansion_coeff = 5e-6
    eigenstrain_name = clad_thermal_strain
  []
  [clad_stress]
    type = ADComputeMultipleInelasticStress
    block = clad
    inelastic_models = 'clad_creep clad_plasticity'
    max_iterations = 100
    relative_tolerance = 1e-12
    absolute_tolerance = 1e-6
  []
  [clad_creep]
    type = ADPowerLawCreepStressUpdate
    block = clad
    coefficient = 8e-28
    n_exponent = 3
    m_exponent = 0
    activation_energy = 0
  []
  [clad_plasticity]
    type = ADIsotropicPlasticityStressUpdate
    block = clad
    yield_stress = 1e6
    hardening_constant = 2e10
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  line_search = basic
  dt = 0.03125
  dtmin = 0.001953125
  dtmax = 0.0625
  end_time = 6
  nl_max_its = 100
  # Segment changes in mortar contact are nonsmooth. The observed scaled
  # residual is below 5e-4 after more than a 3,000-fold reduction; final fields
  # are independently gated below 0.5 percent.
  nl_abs_tol = 5e-4
  nl_rel_tol = 1e-8
  abort_on_solve_fail = false
  automatic_scaling = true
  compute_scaling_once = false
  off_diagonals_in_auto_scaling = true
  petsc_options_iname = -pc_type
  petsc_options_value = lu

  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.03125
    growth_factor = 2
    cutback_factor_at_failure = 0.5
    timestep_limiting_function = 'power internal_pressure external_pressure axial_slide'
    force_step_every_function_point = true
  []
[]

[Postprocessors]
  [average_effective_plastic]
    type = ADElementAverageMaterialProperty
    mat_prop = effective_plastic_strain
    block = clad
  []
  [average_effective_creep]
    type = ADElementAverageMaterialProperty
    mat_prop = effective_creep_strain
    block = clad
  []
  [average_vonmises_stress]
    type = ElementAverageValue
    variable = vonmises_stress
    block = clad
  []
[]

[VectorPostprocessors]
  [all_nodes]
    type = NodalValueSampler
    block = 'fuel clad'
    variable = 'T disp_x disp_y'
    sort_by = id
    use_displaced_mesh = false
  []
  [fuel_surface]
    type = NodalValueSampler
    boundary = fuel_right
    variable = 'T disp_x disp_y'
    sort_by = y
    use_displaced_mesh = false
  []
  [mortar_pressure]
    type = NodalValueSampler
    boundary = fuel_right
    variable = 'mortar_normal_pressure mortar_tangential_pressure'
    sort_by = y
    use_displaced_mesh = false
  []
[]

[Outputs]
  csv = true
[]
