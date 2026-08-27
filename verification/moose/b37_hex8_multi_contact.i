# B3.7 three-dimensional multi-contact friction verification.
# Two independent HEX8 contact pairs are solved together. Each pair uses
# penalty mechanical contact with Coulomb friction and a tangential load.

[Mesh]
  file = b37_hex8_multi_contact_mesh.e
  patch_update_strategy = iteration
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
  [pair_a_tangential_force_y]
  []
  [pair_b_tangential_force_y]
  []
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [primary_a]
        block = primary_a
        strain = SMALL
        volumetric_locking_correction = true
        incremental = true
        add_variables = true
        temperature = T
        use_automatic_differentiation = true
      []
      [secondary_a]
        block = secondary_a
        strain = SMALL
        volumetric_locking_correction = true
        incremental = true
        add_variables = true
        temperature = T
        use_automatic_differentiation = true
      []
      [primary_b]
        block = primary_b
        strain = SMALL
        volumetric_locking_correction = true
        incremental = true
        add_variables = true
        temperature = T
        use_automatic_differentiation = true
      []
      [secondary_b]
        block = secondary_b
        strain = SMALL
        volumetric_locking_correction = true
        incremental = true
        add_variables = true
        temperature = T
        use_automatic_differentiation = true
      []
    []
  []
[]

[Kernels]
  [heat]
    type = ADMatDiffusion
    variable = T
    block = 'primary_a secondary_a primary_b secondary_b'
    diffusivity = thermal_conductivity
  []
[]

[Contact]
  [pair_a]
    primary = primary_a_right
    secondary = secondary_a_left
    formulation = penalty
    model = coulomb
    friction_coefficient = 0.001
    penalty = 1e12
    normalize_penalty = true
  []
  [pair_b]
    primary = primary_b_right
    secondary = secondary_b_left
    formulation = penalty
    model = coulomb
    friction_coefficient = 0.001
    penalty = 1e12
    normalize_penalty = true
  []
[]

[BCs]
  [primary_a_temperature]
    type = DirichletBC
    variable = T
    boundary = primary_a_left
    value = 300
  []
  [primary_a_interface_temperature]
    type = DirichletBC
    variable = T
    boundary = primary_a_right
    value = 300
  []
  [primary_a_x]
    type = ADDirichletBC
    variable = disp_x
    boundary = primary_a_left
    value = 0
  []
  [primary_a_y]
    type = ADDirichletBC
    variable = disp_y
    boundary = primary_a_left
    value = 0
  []
  [primary_a_z]
    type = ADDirichletBC
    variable = disp_z
    boundary = primary_a_left
    value = 0
  []
  [secondary_a_temperature]
    type = DirichletBC
    variable = T
    boundary = secondary_a_left
    value = 300
  []
  [secondary_a_right_temperature]
    type = DirichletBC
    variable = T
    boundary = secondary_a_right
    value = 300
  []
  [secondary_a_x]
    type = ADFunctionDirichletBC
    variable = disp_x
    boundary = secondary_a_right
    function = '-2e-6*t'
  []
  [secondary_a_y]
    type = ADFunctionDirichletBC
    variable = disp_y
    boundary = secondary_a_right
    function = '2e-6*t'
  []
  [secondary_a_z]
    type = ADDirichletBC
    variable = disp_z
    boundary = secondary_a_right
    value = 0
  []
  [primary_b_temperature]
    type = DirichletBC
    variable = T
    boundary = primary_b_left
    value = 300
  []
  [primary_b_interface_temperature]
    type = DirichletBC
    variable = T
    boundary = primary_b_right
    value = 300
  []
  [primary_b_x]
    type = ADDirichletBC
    variable = disp_x
    boundary = primary_b_left
    value = 0
  []
  [primary_b_y]
    type = ADDirichletBC
    variable = disp_y
    boundary = primary_b_left
    value = 0
  []
  [primary_b_z]
    type = ADDirichletBC
    variable = disp_z
    boundary = primary_b_left
    value = 0
  []
  [secondary_b_temperature]
    type = DirichletBC
    variable = T
    boundary = secondary_b_left
    value = 300
  []
  [secondary_b_right_temperature]
    type = DirichletBC
    variable = T
    boundary = secondary_b_right
    value = 300
  []
  [secondary_b_x]
    type = ADFunctionDirichletBC
    variable = disp_x
    boundary = secondary_b_right
    function = '-2e-6*t'
  []
  [secondary_b_y]
    type = ADFunctionDirichletBC
    variable = disp_y
    boundary = secondary_b_right
    function = '2e-6*t'
  []
  [secondary_b_z]
    type = ADDirichletBC
    variable = disp_z
    boundary = secondary_b_right
    value = 0
  []
[]

[Materials]
  [conductivity]
    type = ADGenericConstantMaterial
    prop_names = thermal_conductivity
    prop_values = 10
    block = 'primary_a secondary_a primary_b secondary_b'
  []
  [elasticity]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 1e9
    poissons_ratio = 0
    block = 'primary_a secondary_a primary_b secondary_b'
  []
  [primary_a_stress]
    type = ADComputeLinearElasticStress
    block = primary_a
  []
  [secondary_a_stress]
    type = ADComputeLinearElasticStress
    block = secondary_a
  []
  [primary_b_stress]
    type = ADComputeLinearElasticStress
    block = primary_b
  []
  [secondary_b_stress]
    type = ADComputeLinearElasticStress
    block = secondary_b
  []
[]

[AuxKernels]
  [pair_a_tangential_force_y]
    type = PenetrationAux
    variable = pair_a_tangential_force_y
    boundary = secondary_a_left
    paired_boundary = primary_a_right
    quantity = tangential_force_y
    use_displaced_mesh = true
  []
  [pair_b_tangential_force_y]
    type = PenetrationAux
    variable = pair_b_tangential_force_y
    boundary = secondary_b_left
    paired_boundary = primary_b_right
    quantity = tangential_force_y
    use_displaced_mesh = true
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  start_time = 0
  end_time = 1
  dt = 1
  line_search = basic
  nl_abs_tol = 1e-8
  nl_rel_tol = 1e-10
  nl_max_its = 50
  automatic_scaling = true
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[VectorPostprocessors]
  [all_nodes]
    type = NodalValueSampler
    block = 'primary_a secondary_a primary_b secondary_b'
    variable = 'T disp_x disp_y disp_z'
    sort_by = id
    use_displaced_mesh = false
    execute_on = timestep_end
  []
  [pair_a_contact]
    type = NodalValueSampler
    boundary = secondary_a_left
    variable = 'contact_pressure nodal_area penetration pair_a_tangential_force_y'
    sort_by = id
    use_displaced_mesh = false
    execute_on = timestep_end
  []
  [pair_b_contact]
    type = NodalValueSampler
    boundary = secondary_b_left
    variable = 'contact_pressure nodal_area penetration pair_b_tangential_force_y'
    sort_by = id
    use_displaced_mesh = false
    execute_on = timestep_end
  []
[]

[Outputs]
  file_base = b37_hex8_multi_contact
  csv = true
[]
