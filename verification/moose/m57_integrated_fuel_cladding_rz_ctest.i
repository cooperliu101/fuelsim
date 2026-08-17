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
    ny = 28
    xmin = 0.004121
    xmax = 0.004692
    ymin = -0.001
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
  [stress_q0]
    order = CONSTANT
    family = MONOMIAL
    block = clad
  []
  [stress_q1]
    order = CONSTANT
    family = MONOMIAL
    block = clad
  []
  [stress_q2]
    order = CONSTANT
    family = MONOMIAL
    block = clad
  []
  [stress_q3]
    order = CONSTANT
    family = MONOMIAL
    block = clad
  []
  [plastic_q0]
    order = CONSTANT
    family = MONOMIAL
    block = clad
  []
  [plastic_q1]
    order = CONSTANT
    family = MONOMIAL
    block = clad
  []
  [plastic_q2]
    order = CONSTANT
    family = MONOMIAL
    block = clad
  []
  [plastic_q3]
    order = CONSTANT
    family = MONOMIAL
    block = clad
  []
  [creep_q0]
    order = CONSTANT
    family = MONOMIAL
    block = clad
  []
  [creep_q1]
    order = CONSTANT
    family = MONOMIAL
    block = clad
  []
  [creep_q2]
    order = CONSTANT
    family = MONOMIAL
    block = clad
  []
  [creep_q3]
    order = CONSTANT
    family = MONOMIAL
    block = clad
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
    x = '0 0.0833333333333333 0.25 0.5 0.75 1'
    y = '0 0.4 1 0.8 1.2 1'
  []
  [internal_pressure]
    type = PiecewiseLinear
    x = '0 0.0833333333333333 0.25 0.5 0.75 1'
    y = '0 0.25 1 0.7 1.1 1'
  []
  [external_pressure]
    type = PiecewiseLinear
    x = '0 0.0833333333333333 0.25 0.5 0.75 1'
    y = '0 0.8 1 1.2 0.9 1'
  []
  [axial_slide]
    type = PiecewiseLinear
    x = '0 0.0833333333333333 0.25 0.5 0.75 1'
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
    value = 2e7
    function = power
    use_displaced_mesh = false
  []
[]

[AuxKernels]
  [stress_q0]
    type = ADMaterialRealAux
    variable = stress_q0
    property = vonmises_stress
    selected_qp = 0
    block = clad
  []
  [stress_q1]
    type = ADMaterialRealAux
    variable = stress_q1
    property = vonmises_stress
    selected_qp = 1
    block = clad
  []
  [stress_q2]
    type = ADMaterialRealAux
    variable = stress_q2
    property = vonmises_stress
    selected_qp = 2
    block = clad
  []
  [stress_q3]
    type = ADMaterialRealAux
    variable = stress_q3
    property = vonmises_stress
    selected_qp = 3
    block = clad
  []
  [plastic_q0]
    type = ADMaterialRealAux
    variable = plastic_q0
    property = effective_plastic_strain
    selected_qp = 0
    block = clad
  []
  [plastic_q1]
    type = ADMaterialRealAux
    variable = plastic_q1
    property = effective_plastic_strain
    selected_qp = 1
    block = clad
  []
  [plastic_q2]
    type = ADMaterialRealAux
    variable = plastic_q2
    property = effective_plastic_strain
    selected_qp = 2
    block = clad
  []
  [plastic_q3]
    type = ADMaterialRealAux
    variable = plastic_q3
    property = effective_plastic_strain
    selected_qp = 3
    block = clad
  []
  [creep_q0]
    type = ADMaterialRealAux
    variable = creep_q0
    property = effective_creep_strain
    selected_qp = 0
    block = clad
  []
  [creep_q1]
    type = ADMaterialRealAux
    variable = creep_q1
    property = effective_creep_strain
    selected_qp = 1
    block = clad
  []
  [creep_q2]
    type = ADMaterialRealAux
    variable = creep_q2
    property = effective_creep_strain
    selected_qp = 2
    block = clad
  []
  [creep_q3]
    type = ADMaterialRealAux
    variable = creep_q3
    property = effective_creep_strain
    selected_qp = 3
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
    model = coulomb
    friction_coefficient = 0.002
    penalty = 1e12
    normalize_penalty = true
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
  [internal_pressure_axial]
    type = ADPressure
    variable = disp_y
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
  [external_pressure_axial]
    type = ADPressure
    variable = disp_y
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
  end_time = 1
  nl_max_its = 100
  nl_abs_tol = 1e-8
  nl_rel_tol = 1e-8
  abort_on_solve_fail = false
  automatic_scaling = true
  compute_scaling_once = false
  off_diagonals_in_auto_scaling = true
  petsc_options_iname = -pc_type
  petsc_options_value = lu

  [TimeStepper]
    type = TimeSequenceStepper
    # Reproduce the two committed half steps of each fuelsim adaptive step.
    time_sequence = '0 0.015625 0.031250 0.0572916666666667 0.0833333333333333 0.1145833333333333
      0.1458333333333333 0.1770833333333333 0.2083333333333333 0.2291666666666667 0.25 0.265625
      0.28125 0.3125 0.34375 0.359375 0.375 0.390625 0.40625 0.421875 0.4375 0.453125
      0.46875 0.4765625 0.484375 0.4921875 0.5 0.515625 0.53125 0.5625 0.59375 0.609375
      0.625 0.65625 0.6875 0.703125 0.71875 0.734375 0.75 0.78125 0.8125 0.828125
      0.84375 0.875 0.90625 0.921875 0.9375 0.9453125 0.953125 0.9609375 0.96875
      0.984375 1.0'
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
  [clad_qp_coordinates]
    type = ElementMaterialSampler
    property = thermal_conductivity
    block = clad
    execute_on = timestep_end
  []
  [clad_qp_values]
    type = ElementValueSampler
    variable = 'stress_q0 stress_q1 stress_q2 stress_q3 plastic_q0 plastic_q1 plastic_q2 plastic_q3 creep_q0 creep_q1 creep_q2 creep_q3'
    block = clad
    sort_by = id
    execute_on = timestep_end
  []
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
  [contact_pressure]
    type = NodalValueSampler
    boundary = fuel_right
    variable = 'contact_pressure nodal_area penetration'
    sort_by = y
    use_displaced_mesh = false
  []
[]

[Outputs]
  csv = true
  file_base = m57_integrated_fuel_cladding_rz_ctest
[]
