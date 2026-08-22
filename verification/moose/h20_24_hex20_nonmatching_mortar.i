# H20.24: quadratic mortar-penalty contact on the tracked nonmatching HEX20 mesh.
# The geometry, elasticity, constraints, normal penalty, and final displacement
# match the Fuelsim and Abaqus H20.24 cases.
[Mesh]
  patch_update_strategy = iteration
  [file]
    type = FileMeshGenerator
    file = h20_24_hex20_nonmatching_contact_mesh.e
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
    primary = primary_contact
    secondary = secondary_contact
    formulation = mortar_penalty
    model = frictionless
    penalty = 1e11
    use_dual = true
  []
[]

[AuxVariables]
  [reaction_x]
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
[]

[BCs]
  [primary_x]
    type = ADDirichletBC
    variable = disp_x
    boundary = primary_x0
    value = 0
  []
  [primary_y]
    type = ADDirichletBC
    variable = disp_y
    boundary = primary_x0
    value = 0
  []
  [primary_z]
    type = ADDirichletBC
    variable = disp_z
    boundary = primary_x0
    value = 0
  []
  [secondary_x]
    type = ADFunctionDirichletBC
    variable = disp_x
    boundary = secondary_x2
    function = '-1e-5*t'
  []
  [secondary_y]
    type = ADDirichletBC
    variable = disp_y
    boundary = secondary_y0
    value = 0
  []
  [secondary_z]
    type = ADDirichletBC
    variable = disp_z
    boundary = secondary_z0
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
  [all_nodes]
    type = NodalValueSampler
    block = 'primary secondary'
    variable = 'disp_x disp_y disp_z contact_pressure'
    sort_by = id
    use_displaced_mesh = false
    execute_on = timestep_end
  []
[]

[Postprocessors]
  [reaction_x]
    type = NodalSum
    variable = reaction_x
    boundary = secondary_x2
  []
[]

[Outputs]
  file_base = h20_24_hex20_nonmatching_mortar
  csv = true
[]
