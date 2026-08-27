[Mesh]
  file = b33_hex8_contact_mesh.e
  patch_update_strategy = iteration
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [primary]
        block = primary
        strain = SMALL
        volumetric_locking_correction = true
        incremental = true
        add_variables = true
        use_automatic_differentiation = true
      []
      [secondary]
        block = secondary
        strain = SMALL
        volumetric_locking_correction = true
        incremental = true
        add_variables = true
        use_automatic_differentiation = true
      []
    []
  []
[]

[Contact]
  [interface]
    primary = primary_right
    secondary = secondary_left
    formulation = penalty
    model = coulomb
    friction_coefficient = 0.2
    penalty = 1e13
    normalize_penalty = true
  []
[]

[BCs]
  [primary_x]
    type = ADDirichletBC
    variable = disp_x
    boundary = primary_left
    value = 0
  []
  [primary_y]
    type = ADDirichletBC
    variable = disp_y
    boundary = primary_left
    value = 0
  []
  [primary_z]
    type = ADDirichletBC
    variable = disp_z
    boundary = primary_left
    value = 0
  []
  [secondary_x]
    type = ADFunctionDirichletBC
    variable = disp_x
    boundary = secondary_right
    function = '-2e-4*t'
  []
  [secondary_y]
    type = ADFunctionDirichletBC
    variable = disp_y
    boundary = secondary_right
    function = '5e-5*t'
  []
  [secondary_z]
    type = ADDirichletBC
    variable = disp_z
    boundary = secondary_right
    value = 0
  []
[]

[Materials]
  [primary_elasticity]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 1e9
    poissons_ratio = 0.25
    block = primary
  []
  [secondary_elasticity]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 1e9
    poissons_ratio = 0.25
    block = secondary
  []
  [primary_stress]
    type = ADComputeLinearElasticStress
    block = primary
  []
  [secondary_stress]
    type = ADComputeLinearElasticStress
    block = secondary
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  start_time = 0
  end_time = 1
  dt = 0.1
  nl_abs_tol = 1e-7
  nl_rel_tol = 1e-10
  nl_max_its = 50
  automatic_scaling = true
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[VectorPostprocessors]
  [all_nodes]
    type = NodalValueSampler
    block = 'primary secondary'
    variable = 'disp_x disp_y disp_z'
    sort_by = id
    use_displaced_mesh = false
    execute_on = timestep_end
  []
  [secondary_contact]
    type = NodalValueSampler
    boundary = secondary_left
    variable = 'contact_pressure nodal_area penetration'
    sort_by = id
    use_displaced_mesh = false
    execute_on = timestep_end
  []
[]

[Outputs]
  file_base = b33_hex8_contact
  csv = true
[]
