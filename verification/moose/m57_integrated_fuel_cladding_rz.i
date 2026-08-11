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
  end_time = 6
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
    # Match the two accepted 0.03125 s fuelsim controller steps through four
    # 0.015625 s material updates, then match its remaining 95 accepted 0.0625 s
    # controller steps through 190 material updates of 0.03125 s.
    time_sequence = '0 0.015625 0.031250 0.046875 0.062500 0.093750 0.125000 0.156250 0.187500 0.218750
      0.250000 0.281250 0.312500 0.343750 0.375000 0.406250 0.437500 0.468750 0.500000 0.531250
      0.562500 0.593750 0.625000 0.656250 0.687500 0.718750 0.750000 0.781250 0.812500 0.843750
      0.875000 0.906250 0.937500 0.968750 1.000000 1.031250 1.062500 1.093750 1.125000 1.156250
      1.187500 1.218750 1.250000 1.281250 1.312500 1.343750 1.375000 1.406250 1.437500 1.468750
      1.500000 1.531250 1.562500 1.593750 1.625000 1.656250 1.687500 1.718750 1.750000 1.781250
      1.812500 1.843750 1.875000 1.906250 1.937500 1.968750 2.000000 2.031250 2.062500 2.093750
      2.125000 2.156250 2.187500 2.218750 2.250000 2.281250 2.312500 2.343750 2.375000 2.406250
      2.437500 2.468750 2.500000 2.531250 2.562500 2.593750 2.625000 2.656250 2.687500 2.718750
      2.750000 2.781250 2.812500 2.843750 2.875000 2.906250 2.937500 2.968750 3.000000 3.031250
      3.062500 3.093750 3.125000 3.156250 3.187500 3.218750 3.250000 3.281250 3.312500 3.343750
      3.375000 3.406250 3.437500 3.468750 3.500000 3.531250 3.562500 3.593750 3.625000 3.656250
      3.687500 3.718750 3.750000 3.781250 3.812500 3.843750 3.875000 3.906250 3.937500 3.968750
      4.000000 4.031250 4.062500 4.093750 4.125000 4.156250 4.187500 4.218750 4.250000 4.281250
      4.312500 4.343750 4.375000 4.406250 4.437500 4.468750 4.500000 4.531250 4.562500 4.593750
      4.625000 4.656250 4.687500 4.718750 4.750000 4.781250 4.812500 4.843750 4.875000 4.906250
      4.937500 4.968750 5.000000 5.031250 5.062500 5.093750 5.125000 5.156250 5.187500 5.218750
      5.250000 5.281250 5.312500 5.343750 5.375000 5.406250 5.437500 5.468750 5.500000 5.531250
      5.562500 5.593750 5.625000 5.656250 5.687500 5.718750 5.750000 5.781250 5.812500 5.843750
      5.875000 5.906250 5.937500 5.968750 6.000000'
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
[]
