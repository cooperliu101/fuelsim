[Mesh]
  patch_update_strategy = iteration
  [file]
    type = FileMeshGenerator
    file = h20_16_hex20_contact_mesh.e
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Problem]
  extra_tag_vectors = ref
[]

[Variables]
  [disp_x]
    order = SECOND
  []
  [disp_y]
    order = SECOND
  []
  [disp_z]
    order = SECOND
  []
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [primary]
        block = primary
        strain = SMALL
        incremental = true
        use_automatic_differentiation = true
        extra_vector_tags = ref
      []
      [secondary]
        block = secondary
        strain = SMALL
        incremental = true
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
    model = frictionless
    # This is the literal normalize_penalty switch-off diagnostic. Without
    # normalization, this value is a nodal spring stiffness in N/m rather than
    # the surface penalty density in Pa/m used by the normalized comparison.
    penalty = 1e13
    normalize_penalty = false
  []
[]

[AuxVariables]
  [reaction_x]
    order = SECOND
  []
  [reaction_y]
    order = SECOND
  []
  [reaction_z]
    order = SECOND
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
    function = '0'
  []
  [secondary_z]
    type = ADFunctionDirichletBC
    variable = disp_z
    boundary = secondary_right
    function = '0'
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
  nl_abs_tol = 1e-10
  nl_rel_tol = 1e-10
  nl_max_its = 50
  automatic_scaling = true
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  [Quadrature]
    order = FIFTH
  []
[]

[VectorPostprocessors]
  [displacement_nodes]
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
  file_base = h20_22_hex20_unnormalized_mechanical
  csv = true
[]
