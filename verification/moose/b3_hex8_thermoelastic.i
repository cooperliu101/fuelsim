[Mesh]
  file = b3_hex8_mesh.e
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Variables]
  [T]
    initial_condition = 300
  []
[]

[AuxVariables]
  [stress_xx]
    family = MONOMIAL
    order = CONSTANT
  []
  [stress_yy]
    family = MONOMIAL
    order = CONSTANT
  []
  [stress_zz]
    family = MONOMIAL
    order = CONSTANT
  []
  [stress_xy]
    family = MONOMIAL
    order = CONSTANT
  []
  [stress_yz]
    family = MONOMIAL
    order = CONSTANT
  []
  [stress_xz]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [solid]
        block = solid
        strain = SMALL
        add_variables = true
        temperature = T
        eigenstrain_names = thermal_strain
        use_automatic_differentiation = true
      []
    []
  []
[]

[Kernels]
  [heat_conduction]
    type = ADMatDiffusion
    variable = T
    block = solid
    diffusivity = thermal_conductivity
  []
[]

[AuxKernels]
  [stress_xx]
    type = ADRankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    index_i = 0
    index_j = 0
    block = solid
    use_displaced_mesh = false
  []
  [stress_yy]
    type = ADRankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 1
    index_j = 1
    block = solid
    use_displaced_mesh = false
  []
  [stress_zz]
    type = ADRankTwoAux
    rank_two_tensor = stress
    variable = stress_zz
    index_i = 2
    index_j = 2
    block = solid
    use_displaced_mesh = false
  []
  [stress_xy]
    type = ADRankTwoAux
    rank_two_tensor = stress
    variable = stress_xy
    index_i = 0
    index_j = 1
    block = solid
    use_displaced_mesh = false
  []
  [stress_yz]
    type = ADRankTwoAux
    rank_two_tensor = stress
    variable = stress_yz
    index_i = 1
    index_j = 2
    block = solid
    use_displaced_mesh = false
  []
  [stress_xz]
    type = ADRankTwoAux
    rank_two_tensor = stress
    variable = stress_xz
    index_i = 0
    index_j = 2
    block = solid
    use_displaced_mesh = false
  []
[]

[BCs]
  [temperature]
    type = DirichletBC
    variable = T
    boundary = solid_left
    value = 300
  []
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
  [traction_x]
    type = ADFunctionNeumannBC
    variable = disp_x
    boundary = solid_right
    function = 1e6
    use_displaced_mesh = false
  []
[]

[Materials]
  [conductivity]
    type = ADGenericConstantMaterial
    prop_names = thermal_conductivity
    prop_values = 10
    block = solid
  []
  [elasticity]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 1e9
    poissons_ratio = 0.25
    block = solid
  []
  [thermal_expansion]
    type = ADComputeThermalExpansionEigenstrain
    temperature = T
    stress_free_temperature = 300
    thermal_expansion_coeff = 1e-5
    eigenstrain_name = thermal_strain
    block = solid
  []
  [stress]
    type = ADComputeLinearElasticStress
    block = solid
  []
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  nl_abs_tol = 1e-10
  nl_rel_tol = 1e-12
  nl_max_its = 20
  automatic_scaling = true
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
  [element_stress]
    type = ElementValueSampler
    variable = 'stress_xx stress_yy stress_zz stress_xy stress_yz stress_xz'
    sort_by = id
    use_displaced_mesh = false
  []
[]

[Outputs]
  file_base = b3_hex8_thermoelastic
  csv = true
[]
