[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 2
  ny = 2
  xmin = 0.004
  xmax = 0.005
  ymin = 0
  ymax = 0.001
  coord_type = RZ
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[AuxVariables]
  [T]
    initial_condition = 600
  []
[]

[Physics/SolidMechanics/QuasiStatic]
  [solid]
    strain = FINITE
    add_variables = true
    use_automatic_differentiation = true
  []
[]

[BCs]
  [right_radial]
    type = ADDirichletBC
    variable = disp_x
    boundary = right
    value = 0
  []
  [bottom_axial]
    type = ADDirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
  []
  [inner_pressure]
    type = ADPressure
    variable = disp_x
    boundary = left
    factor = 1.0e8
  []
  [top_pressure]
    type = ADPressure
    variable = disp_y
    boundary = top
    factor = 1.0e8
  []
[]

[Materials]
  [elasticity]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 2.0e11
    poissons_ratio = 0.3
  []
  [stress]
    type = ADComputeMultipleInelasticStress
    inelastic_models = ''
  []
[]

[Executioner]
  type = Transient
  dt = 1
  end_time = 1
  solve_type = NEWTON
  line_search = bt
  nl_abs_tol = 1.0e-10
  nl_rel_tol = 1.0e-10
  nl_max_its = 80
  automatic_scaling = true
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[VectorPostprocessors]
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
