[Mesh]
  file = b3_hex8_mesh.e
[]
[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]
[AuxVariables]
  [T]
    initial_condition = 600
  []
  [effective_plastic_strain]
    family = MONOMIAL
    order = CONSTANT
  []
  [effective_creep_strain]
    family = MONOMIAL
    order = CONSTANT
  []
[]
[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [solid]
        block = solid
        strain = FINITE
        incremental = true
        add_variables = true
        use_automatic_differentiation = true
        generate_output = 'stress_xx stress_yy stress_zz stress_xy stress_yz stress_xz'
      []
    []
  []
[]
[Functions]
  [pull]
    type = ParsedFunction
    expression = '0.004*t'
  []
[]
[AuxKernels]
  [effective_plastic_strain]
    type = ADMaterialRealAux
    variable = effective_plastic_strain
    property = effective_plastic_strain
    block = solid
  []
  [effective_creep_strain]
    type = ADMaterialRealAux
    variable = effective_creep_strain
    property = effective_creep_strain
    block = solid
  []
[]
[BCs]
  [fix_x]
    type = ADDirichletBC
    variable = disp_x
    boundary = solid_left
    value = 0
  []
  [fix_y]
    type = ADDirichletBC
    variable = disp_y
    boundary = solid_bottom
    value = 0
  []
  [fix_z]
    type = ADDirichletBC
    variable = disp_z
    boundary = solid_back
    value = 0
  []
  [pull_x]
    type = ADFunctionDirichletBC
    variable = disp_x
    boundary = solid_right
    function = pull
  []
[]
[Materials]
  [elasticity]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 2e11
    poissons_ratio = 0.3
    block = solid
  []
  [stress]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'creep plasticity'
    max_iterations = 100
    relative_tolerance = 1e-12
    absolute_tolerance = 1e-6
    block = solid
  []
  [creep]
    type = ADPowerLawCreepStressUpdate
    coefficient = 1e-28
    n_exponent = 3
    m_exponent = 0
    activation_energy = 0
    block = solid
  []
  [plasticity]
    type = ADIsotropicPlasticityStressUpdate
    yield_stress = 2e8
    hardening_constant = 2e9
    block = solid
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
  dt = 0.1
  end_time = 1
  nl_abs_tol = 1e-10
  nl_rel_tol = 1e-10
  nl_max_its = 40
  automatic_scaling = true
  compute_scaling_once = false
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]
[VectorPostprocessors]
  [all_nodes]
    type = NodalValueSampler
    variable = 'T disp_x disp_y disp_z'
    sort_by = id
    use_displaced_mesh = false
  []
  [element_state]
    type = ElementValueSampler
    variable = 'stress_xx stress_yy stress_zz stress_xy stress_yz stress_xz effective_plastic_strain effective_creep_strain'
    sort_by = id
    use_displaced_mesh = false
  []
[]
[Outputs]
  file_base = b32_hex8_finite_coupled
  csv = true
[]
