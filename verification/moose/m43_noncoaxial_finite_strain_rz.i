# M4.3 distorted multi-element noncoaxial finite-strain history verification.
# The annulus is stretched, sheared beyond 25 degrees polar rotation, reversed,
# and reloaded while current-configuration pressure and component traction act.
# MOOSE Taylor/Rashid defaults remain unchanged.

[Mesh]
  coord_type = RZ
  [base]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 2
    ny = 2
    xmin = 0.010
    xmax = 0.012
    ymin = 0.0
    ymax = 0.002
  []
  [distort]
    type = ParsedNodeTransformGenerator
    input = base
    x_function = 'x + 1.5e8*(x-0.010)*(0.012-x)*y*(0.002-y)'
    y_function = 'y - 1.0e8*(x-0.010)*(0.012-x)*y*(0.002-y)'
  []
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
  [stress_rr]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_zz]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_hoop]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_rz]
    order = CONSTANT
    family = MONOMIAL
  []
  [elastic_rr]
    order = CONSTANT
    family = MONOMIAL
  []
  [elastic_zz]
    order = CONSTANT
    family = MONOMIAL
  []
  [elastic_hoop]
    order = CONSTANT
    family = MONOMIAL
  []
  [elastic_rz]
    order = CONSTANT
    family = MONOMIAL
  []
  [effective_plastic]
    order = CONSTANT
    family = MONOMIAL
  []
  [effective_creep]
    order = CONSTANT
    family = MONOMIAL
  []
  [combined_rr]
    order = CONSTANT
    family = MONOMIAL
  []
  [combined_zz]
    order = CONSTANT
    family = MONOMIAL
  []
  [combined_hoop]
    order = CONSTANT
    family = MONOMIAL
  []
  [combined_rz]
    order = CONSTANT
    family = MONOMIAL
  []
  [sample_time]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Physics/SolidMechanics/QuasiStatic]
  [solid]
    strain = FINITE
    add_variables = true
    use_automatic_differentiation = true
  []
[]

[Functions]
  [axial_path]
    type = PiecewiseLinear
    x = '0 1 2 3 4 5'
    y = '0 2e-4 2e-4 -1e-4 -1e-4 6e-5'
  []
  [shear_path]
    type = PiecewiseLinear
    x = '0 1 2 3 4 5'
    y = '0 0 2.0e-3 2.0e-3 -1.5e-3 -1.5e-3'
  []
  [load_ramp]
    type = PiecewiseLinear
    x = '0 5'
    y = '0 1'
  []
  [traction_ramp]
    type = PiecewiseLinear
    x = '0 5'
    y = '0 1e6'
  []
  [output_time]
    type = ParsedFunction
    expression = 't'
  []
[]

[Kernels]
  [heat_time]
    type = HeatConductionTimeDerivative
    variable = T
    use_displaced_mesh = false
  []
  [heat_conduction]
    type = HeatConduction
    variable = T
    use_displaced_mesh = false
  []
[]

[AuxKernels]
  # CONSTANT MONOMIAL AuxKernels integrate all four QPs over the reference RZ
  # element volume.  Keep both the measure and execution point explicit.
  [stress_rr]
    type = ADRankTwoAux
    rank_two_tensor = stress
    variable = stress_rr
    index_i = 0
    index_j = 0
    use_displaced_mesh = false
    execute_on = timestep_end
  []
  [stress_zz]
    type = ADRankTwoAux
    rank_two_tensor = stress
    variable = stress_zz
    index_i = 1
    index_j = 1
    use_displaced_mesh = false
    execute_on = timestep_end
  []
  [stress_hoop]
    type = ADRankTwoAux
    rank_two_tensor = stress
    variable = stress_hoop
    index_i = 2
    index_j = 2
    use_displaced_mesh = false
    execute_on = timestep_end
  []
  [stress_rz]
    type = ADRankTwoAux
    rank_two_tensor = stress
    variable = stress_rz
    index_i = 0
    index_j = 1
    use_displaced_mesh = false
    execute_on = timestep_end
  []
  [elastic_rr]
    type = ADRankTwoAux
    rank_two_tensor = elastic_strain
    variable = elastic_rr
    index_i = 0
    index_j = 0
    use_displaced_mesh = false
    execute_on = timestep_end
  []
  [elastic_zz]
    type = ADRankTwoAux
    rank_two_tensor = elastic_strain
    variable = elastic_zz
    index_i = 1
    index_j = 1
    use_displaced_mesh = false
    execute_on = timestep_end
  []
  [elastic_hoop]
    type = ADRankTwoAux
    rank_two_tensor = elastic_strain
    variable = elastic_hoop
    index_i = 2
    index_j = 2
    use_displaced_mesh = false
    execute_on = timestep_end
  []
  [elastic_rz]
    type = ADRankTwoAux
    rank_two_tensor = elastic_strain
    variable = elastic_rz
    index_i = 0
    index_j = 1
    use_displaced_mesh = false
    execute_on = timestep_end
  []
  [effective_plastic]
    type = ADMaterialRealAux
    property = effective_plastic_strain
    variable = effective_plastic
    use_displaced_mesh = false
    execute_on = timestep_end
  []
  [effective_creep]
    type = ADMaterialRealAux
    property = effective_creep_strain
    variable = effective_creep
    use_displaced_mesh = false
    execute_on = timestep_end
  []
  [combined_rr]
    type = ADRankTwoAux
    rank_two_tensor = combined_inelastic_strain
    variable = combined_rr
    index_i = 0
    index_j = 0
    use_displaced_mesh = false
    execute_on = timestep_end
  []
  [combined_zz]
    type = ADRankTwoAux
    rank_two_tensor = combined_inelastic_strain
    variable = combined_zz
    index_i = 1
    index_j = 1
    use_displaced_mesh = false
    execute_on = timestep_end
  []
  [combined_hoop]
    type = ADRankTwoAux
    rank_two_tensor = combined_inelastic_strain
    variable = combined_hoop
    index_i = 2
    index_j = 2
    use_displaced_mesh = false
    execute_on = timestep_end
  []
  [combined_rz]
    type = ADRankTwoAux
    rank_two_tensor = combined_inelastic_strain
    variable = combined_rz
    index_i = 0
    index_j = 1
    use_displaced_mesh = false
    execute_on = timestep_end
  []
  [sample_time]
    type = FunctionAux
    variable = sample_time
    function = output_time
    use_displaced_mesh = false
    execute_on = timestep_end
  []
[]

[BCs]
  [bottom_r]
    type = ADDirichletBC
    variable = disp_x
    boundary = bottom
    value = 0
  []
  [bottom_z]
    type = ADDirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
  []
  [top_r]
    type = ADFunctionDirichletBC
    variable = disp_x
    boundary = top
    function = shear_path
  []
  [top_z]
    type = ADFunctionDirichletBC
    variable = disp_y
    boundary = top
    function = axial_path
  []
  [inner_pressure]
    type = ADPressure
    variable = disp_x
    boundary = left
    factor = 1e6
    function = load_ramp
    use_displaced_mesh = true
  []
  # A pressure boundary condition acts on one displacement equation.  The
  # deformed inner edge is inclined during shear, so retain both components
  # of the current outward normal.
  [inner_pressure_axial]
    type = ADPressure
    variable = disp_y
    boundary = left
    factor = 1e6
    function = load_ramp
    use_displaced_mesh = true
  []
  [outer_axial_traction]
    type = ADFunctionNeumannBC
    variable = disp_y
    boundary = right
    function = traction_ramp
    use_displaced_mesh = true
  []
[]

[Materials]
  [thermal]
    type = GenericConstantMaterial
    prop_names = 'thermal_conductivity density specific_heat'
    prop_values = '1 1 1'
  []
  [elasticity]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 2e11
    poissons_ratio = 0.3
  []
  [stress]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'creep plasticity'
    max_iterations = 100
    relative_tolerance = 1e-12
    absolute_tolerance = 1e-6
  []
  [creep]
    type = ADPowerLawCreepStressUpdate
    coefficient = 1e-30
    n_exponent = 3
    m_exponent = 0
    activation_energy = 0
  []
  [plasticity]
    type = ADIsotropicPlasticityStressUpdate
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
  dt = 0.05
  end_time = 5
  nl_max_its = 40
  nl_abs_tol = 1e-10
  nl_rel_tol = 1e-10
  abort_on_solve_fail = true
  automatic_scaling = true
  compute_scaling_once = false
  off_diagonals_in_auto_scaling = true
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[VectorPostprocessors]
  [element_history]
    type = ElementValueSampler
    variable = 'sample_time stress_rr stress_zz stress_hoop stress_rz elastic_rr elastic_zz elastic_hoop elastic_rz effective_plastic effective_creep combined_rr combined_zz combined_hoop combined_rz'
    sort_by = sample_time
    contains_complete_history = true
    execute_on = timestep_end
    use_displaced_mesh = false
  []
  [all_nodes]
    type = NodalValueSampler
    variable = 'T disp_x disp_y'
    sort_by = id
    use_displaced_mesh = false
  []
[]

[Outputs]
  csv = true
[]
