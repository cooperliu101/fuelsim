[Mesh]
  file = b36_plate_meat_clad_mesh.e
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
    block = 'meat clad'
  []
  [effective_creep_strain]
    family = MONOMIAL
    order = CONSTANT
    block = 'meat clad'
  []
[]

[Physics/SolidMechanics/QuasiStatic]
  [meat]
    block = meat
    strain = SMALL
    volumetric_locking_correction = true
    incremental = true
    add_variables = true
    temperature = T
    eigenstrain_names = meat_thermal_strain
    use_automatic_differentiation = true
    generate_output = 'stress_xx stress_yy stress_zz stress_xy stress_yz stress_xz'
  []
  [clad]
    block = clad
    strain = SMALL
    volumetric_locking_correction = true
    incremental = true
    add_variables = true
    temperature = T
    eigenstrain_names = clad_thermal_strain
    use_automatic_differentiation = true
    generate_output = 'stress_xx stress_yy stress_zz stress_xy stress_yz stress_xz'
  []
[]

[Kernels]
  [heat_time]
    type = HeatConductionTimeDerivative
    variable = T
    block = 'meat clad'
    use_displaced_mesh = false
  []
  [heat_conduction]
    type = HeatConduction
    variable = T
    block = 'meat clad'
    use_displaced_mesh = false
  []
  [heat_source_meat]
    type = BodyForce
    variable = T
    block = meat
    value = 3.291e5
    use_displaced_mesh = false
  []
  [heat_source_clad]
    type = BodyForce
    variable = T
    block = clad
    value = 2.145e5
    use_displaced_mesh = false
  []
[]

[AuxKernels]
  [effective_plastic_strain]
    type = ADMaterialRealAux
    variable = effective_plastic_strain
    property = effective_plastic_strain
    block = 'meat clad'
  []
  [effective_creep_strain]
    type = ADMaterialRealAux
    variable = effective_creep_strain
    property = effective_creep_strain
    block = 'meat clad'
  []
[]

[Functions]
  [pull]
    type = ParsedFunction
    expression = '4.8e-5*t'
  []
[]

[BCs]
  [fix_x]
    type = ADDirichletBC
    variable = disp_x
    boundary = plate_left
    value = 0
  []
  [fix_y]
    type = ADDirichletBC
    variable = disp_y
    boundary = plate_bottom
    value = 0
  []
  [fix_z]
    type = ADDirichletBC
    variable = disp_z
    boundary = plate_back
    value = 0
  []
  [pull_x]
    type = ADFunctionDirichletBC
    variable = disp_x
    boundary = plate_right
    function = pull
  []
[]

[Materials]
  [meat_thermal]
    type = GenericConstantMaterial
    block = meat
    prop_names = 'thermal_conductivity density specific_heat'
    prop_values = '3 10970 300'
  []
  [meat_elasticity]
    type = ADComputeIsotropicElasticityTensor
    block = meat
    youngs_modulus = 2e11
    poissons_ratio = 0.316
  []
  [meat_expansion]
    type = ADComputeThermalExpansionEigenstrain
    block = meat
    temperature = T
    stress_free_temperature = 600
    thermal_expansion_coeff = 1e-5
    eigenstrain_name = meat_thermal_strain
  []
  [meat_stress]
    type = ADComputeMultipleInelasticStress
    block = meat
    inelastic_models = 'meat_creep meat_plasticity'
    perform_finite_strain_rotations = false
    max_iterations = 100
    relative_tolerance = 1e-12
    absolute_tolerance = 1e-6
  []
  [meat_creep]
    type = ADPowerLawCreepStressUpdate
    block = meat
    coefficient = 1e-31
    n_exponent = 3
    m_exponent = 0
    activation_energy = 0
  []
  [meat_plasticity]
    type = ADIsotropicPlasticityStressUpdate
    block = meat
    yield_stress = 2e8
    hardening_constant = 2e9
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
    youngs_modulus = 2e11
    poissons_ratio = 0.316
  []
  [clad_expansion]
    type = ADComputeThermalExpansionEigenstrain
    block = clad
    temperature = T
    stress_free_temperature = 600
    thermal_expansion_coeff = 1e-5
    eigenstrain_name = clad_thermal_strain
  []
  [clad_stress]
    type = ADComputeMultipleInelasticStress
    block = clad
    inelastic_models = 'clad_creep clad_plasticity'
    perform_finite_strain_rotations = false
    max_iterations = 100
    relative_tolerance = 1e-12
    absolute_tolerance = 1e-6
  []
  [clad_creep]
    type = ADPowerLawCreepStressUpdate
    block = clad
    coefficient = 1e-31
    n_exponent = 3
    m_exponent = 0
    activation_energy = 0
  []
  [clad_plasticity]
    type = ADIsotropicPlasticityStressUpdate
    block = clad
    yield_stress = 2e8
    hardening_constant = 2e9
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
  dt = 0.0125
  end_time = 0.5
  nl_abs_tol = 1e-8
  nl_rel_tol = 1e-10
  nl_max_its = 30
  automatic_scaling = false
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[VectorPostprocessors]
  [all_nodes]
    type = NodalValueSampler
    block = 'meat clad'
    variable = 'T disp_x disp_y disp_z'
    sort_by = id
    use_displaced_mesh = false
    execute_on = timestep_end
  []
  [element_state]
    type = ElementValueSampler
    block = 'meat clad'
    variable = 'stress_xx stress_yy stress_zz stress_xy stress_yz stress_xz effective_plastic_strain effective_creep_strain'
    sort_by = id
    use_displaced_mesh = false
    execute_on = timestep_end
  []
[]

[Outputs]
  file_base = b36_plate_meat_clad
  csv = true
  console = false
[]
