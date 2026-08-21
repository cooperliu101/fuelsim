[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 1
  ny = 1
  nz = 1
  elem_type = HEX20
  xmax = 1
  ymax = 1
  zmax = 1
[]
[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]
[AuxVariables]
  [T]
    initial_condition = 600
  []
[]
[Physics/SolidMechanics/QuasiStatic]
  [solid]
    strain = SMALL
    incremental = true
    add_variables = true
    temperature = T
    eigenstrain_names = thermal_strain
    use_automatic_differentiation = true
  []
[]
[AuxVariables]
  [stress_xx] family = MONOMIAL order = CONSTANT []
  [stress_yy] family = MONOMIAL order = CONSTANT []
  [stress_zz] family = MONOMIAL order = CONSTANT []
  [stress_xy] family = MONOMIAL order = CONSTANT []
  [stress_yz] family = MONOMIAL order = CONSTANT []
  [stress_xz] family = MONOMIAL order = CONSTANT []
  [effective_plastic_strain] family = MONOMIAL order = CONSTANT []
  [effective_creep_strain] family = MONOMIAL order = CONSTANT []
[]
[AuxKernels]
  [sxx] type = ADRankTwoAux rank_two_tensor = stress variable = stress_xx index_i = 0 index_j = 0 []
  [syy] type = ADRankTwoAux rank_two_tensor = stress variable = stress_yy index_i = 1 index_j = 1 []
  [szz] type = ADRankTwoAux rank_two_tensor = stress variable = stress_zz index_i = 2 index_j = 2 []
  [sxy] type = ADRankTwoAux rank_two_tensor = stress variable = stress_xy index_i = 0 index_j = 1 []
  [syz] type = ADRankTwoAux rank_two_tensor = stress variable = stress_yz index_i = 1 index_j = 2 []
  [sxz] type = ADRankTwoAux rank_two_tensor = stress variable = stress_xz index_i = 0 index_j = 2 []
  [ep] type = ADMaterialRealAux variable = effective_plastic_strain property = effective_plastic_strain []
  [ec] type = ADMaterialRealAux variable = effective_creep_strain property = effective_creep_strain []
[]
[BCs]
  [fix_x]
    type = ADDirichletBC
    variable = disp_x
    boundary = left
    value = 0
  []
  [fix_y]
    type = ADDirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
  []
  [fix_z]
    type = ADDirichletBC
    variable = disp_z
    boundary = back
    value = 0
  []
  [pull_x]
    type = ADFunctionDirichletBC
    variable = disp_x
    boundary = right
    function = axial_pull
  []
[]
[Functions]
  [axial_pull]
    type = ParsedFunction
    expression = '0.004*t'
  []
[]
[Materials]
  [thermal]
    type = ADGenericConstantMaterial
    prop_names = 'thermal_conductivity density specific_heat'
    prop_values = '1 1 1'
  []
  [elasticity]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 2e11
    poissons_ratio = 0.3
  []
  [thermal_expansion]
    type = ADComputeThermalExpansionEigenstrain
    temperature = T
    stress_free_temperature = 600
    thermal_expansion_coeff = 0
    eigenstrain_name = thermal_strain
  []
  [stress]
    type = ADComputeMultipleInelasticStress
    inelastic_models = plasticity
    perform_finite_strain_rotations = false
  []
  [plasticity]
    type = ADIsotropicPlasticityStressUpdate
    yield_stress = 2e8
    hardening_constant = 2e9
  []
  [inactive_creep]
    type = ADGenericConstantMaterial
    prop_names = effective_creep_strain
    prop_values = 0
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
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]
[VectorPostprocessors]
  [all_nodes]
    type = NodalValueSampler
    variable = T
    sort_by = id
    use_displaced_mesh = false
  []
  [displacement_nodes]
    type = NodalValueSampler
    variable = 'disp_x disp_y disp_z'
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
  file_base = h20_07_inelastic_plastic
  csv = true
  exodus = true
[]
