[Mesh]
  file = b34_hex8_sliding_contact_mesh.e
  patch_update_strategy = iteration
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Problem]
  extra_tag_vectors = ref
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
        extra_vector_tags = ref
      []
      [secondary]
        block = secondary
        strain = SMALL
        volumetric_locking_correction = true
        incremental = true
        add_variables = true
        use_automatic_differentiation = true
        extra_vector_tags = ref
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
    friction_coefficient = 0.001
    penalty = 1e13
    normalize_penalty = true
  []
[]

[AuxVariables]
  [reaction_x]
  []
  [reaction_y]
  []
  [reaction_z]
  []
[]

[AuxKernels]
  [reaction_x]
    type = ReactionForceAux
    vector_tag = ref
    v = disp_x
    variable = reaction_x
  []
  [reaction_y]
    type = ReactionForceAux
    vector_tag = ref
    v = disp_y
    variable = reaction_y
  []
  [reaction_z]
    type = ReactionForceAux
    vector_tag = ref
    v = disp_z
    variable = reaction_z
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
    function = '-2e-5*t'
  []
  [secondary_y]
    type = ADFunctionDirichletBC
    variable = disp_y
    boundary = secondary_right
    function = '2e-6*t'
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
    poissons_ratio = 0
    block = primary
  []
  [secondary_elasticity]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 1e9
    poissons_ratio = 0
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
  dt = 1
  nl_abs_tol = 1e-8
  nl_rel_tol = 1e-10
  nl_max_its = 50
  automatic_scaling = true
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Postprocessors]
  [reaction_x]
    type = NodalSum
    variable = reaction_x
    boundary = secondary_right
  []
  [reaction_y]
    type = NodalSum
    variable = reaction_y
    boundary = secondary_right
  []
  [reaction_z]
    type = NodalSum
    variable = reaction_z
    boundary = secondary_right
  []
[]

[Outputs]
  file_base = b34_hex8_sliding_contact
  csv = true
[]
