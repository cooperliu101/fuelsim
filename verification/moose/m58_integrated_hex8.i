# Three-dimensional counterpart of the M5.7 integrated path. A current-
# thickness has a controlled nonzero extension while the remaining loads
# exercise finite strain, coupled inelasticity, thermal contact, large sliding,
# and vector friction.

[Mesh]
  file = m58_integrated_hex8_mesh.e
  patch_update_strategy = iteration
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Variables]
  [T]
    initial_condition = 600
  []
[]

[AuxVariables]
  [effective_plastic_strain]
    family = MONOMIAL
    order = CONSTANT
    block = clad
  []
  [effective_creep_strain]
    family = MONOMIAL
    order = CONSTANT
    block = clad
  []
[]

[Physics/SolidMechanics/QuasiStatic]
  [fuel]
    block = fuel
    strain = FINITE
    incremental = true
    add_variables = true
    temperature = T
    eigenstrain_names = fuel_thermal_strain
    use_automatic_differentiation = true
  []
  [clad]
    block = clad
    strain = FINITE
    incremental = true
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
    x = '0 0.2 0.5 0.8 1'
    y = '0 0.5 1 1.2 1'
  []
  [internal_pressure]
    type = PiecewiseLinear
    x = '0 0.2 0.5 0.8 1'
    y = '0 0.3 1 0.8 1'
  []
  [external_pressure]
    type = PiecewiseLinear
    x = '0 0.2 0.5 0.8 1'
    y = '0 0.8 1 1.2 1'
  []
  [axial_slide]
    type = PiecewiseLinear
    x = '0 0.2 0.5 0.8 1'
    y = '0 6e-5 3e-4 4.8e-4 6e-4'
  []
  [cladding_axial_slide]
    type = PiecewiseLinear
    x = '0 0.2 0.5 0.8 1'
    y = '0 2e-6 1e-5 1.6e-5 2e-5'
  []
  [thickness_stretch]
    type = PiecewiseLinear
    x = '0 0.2 0.5 0.8 1'
    y = '0 6e-7 3e-6 4.8e-6 6e-6'
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
    value = 2e7
    function = power
    use_displaced_mesh = false
  []
[]

[AuxKernels]
  [effective_plastic_strain]
    type = ADMaterialRealAux
    variable = effective_plastic_strain
    property = effective_plastic_strain
    block = clad
  []
  [effective_creep_strain]
    type = ADMaterialRealAux
    variable = effective_creep_strain
    property = effective_creep_strain
    block = clad
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
    model = coulomb
    friction_coefficient = 0.002
    penalty = 1e12
    normalize_penalty = true
  []
[]

[BCs]
  [fuel_left]
    type = ADDirichletBC
    variable = disp_x
    boundary = fuel_left
    value = 0
  []
  [fuel_back]
    type = ADDirichletBC
    variable = disp_z
    boundary = fuel_back
    value = 0
  []
  [fuel_top_slide]
    type = ADFunctionDirichletBC
    variable = disp_y
    boundary = fuel_top
    function = axial_slide
  []
  [fuel_front_stretch]
    type = ADFunctionDirichletBC
    variable = disp_z
    boundary = fuel_front
    function = thickness_stretch
  []
  [clad_bottom]
    type = ADDirichletBC
    variable = disp_y
    boundary = clad_bottom
    value = 0
  []
  [clad_back]
    type = ADDirichletBC
    variable = disp_z
    boundary = clad_back
    value = 0
  []
  [clad_top_radial_support]
    type = ADDirichletBC
    variable = disp_x
    boundary = clad_top
    value = 0
  []
  [clad_top_slide]
    type = ADFunctionDirichletBC
    variable = disp_y
    boundary = clad_top
    function = cladding_axial_slide
  []
  [clad_front_stretch]
    type = ADFunctionDirichletBC
    variable = disp_z
    boundary = clad_front
    function = thickness_stretch
  []
  [clad_outer_temperature]
    type = DirichletBC
    variable = T
    boundary = clad_right
    value = 600
  []
  [internal_pressure_x]
    type = ADPressure
    variable = disp_x
    boundary = clad_left
    factor = 5e5
    function = internal_pressure
    use_displaced_mesh = true
  []
  [internal_pressure_y]
    type = ADPressure
    variable = disp_y
    boundary = clad_left
    factor = 5e5
    function = internal_pressure
    use_displaced_mesh = true
  []
  [internal_pressure_z]
    type = ADPressure
    variable = disp_z
    boundary = clad_left
    factor = 5e5
    function = internal_pressure
    use_displaced_mesh = true
  []
  [external_pressure_x]
    type = ADPressure
    variable = disp_x
    boundary = clad_right
    factor = 2e6
    function = external_pressure
    use_displaced_mesh = true
  []
  [external_pressure_y]
    type = ADPressure
    variable = disp_y
    boundary = clad_right
    factor = 2e6
    function = external_pressure
    use_displaced_mesh = true
  []
  [external_pressure_z]
    type = ADPressure
    variable = disp_z
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

[Dampers]
  [contact_slip]
    type = ContactSlipDamper
    primary = clad_left
    secondary = fuel_right
    max_iterative_slip = 2e-6
    min_damping = 1e-6
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  line_search = basic
  start_time = 0
  end_time = 1
  dt = 0.05
  nl_max_its = 20
  nl_abs_tol = 1e-7
  nl_rel_tol = 1e-8
  abort_on_solve_fail = true
  automatic_scaling = true
  compute_scaling_once = false
  off_diagonals_in_auto_scaling = true
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
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
    variable = 'T disp_x disp_y disp_z'
    sort_by = id
    use_displaced_mesh = false
    execute_on = timestep_end
  []
  [contact_pressure]
    type = NodalValueSampler
    boundary = fuel_right
    variable = 'contact_pressure nodal_area penetration'
    sort_by = id
    use_displaced_mesh = false
    execute_on = timestep_end
  []
  [clad_state]
    type = ElementValueSampler
    block = clad
    variable = 'vonmises_stress effective_plastic_strain effective_creep_strain'
    sort_by = id
    use_displaced_mesh = false
    execute_on = timestep_end
  []
[]

[Outputs]
  file_base = m58_integrated_hex8
  csv = true
  console = false
[]
