[Mesh]
  type = FileMesh
  file = m22_coupled_plastic_creep_traction_rz_mesh.e
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
    generate_output = 'stress_xx stress_yy stress_zz'
  []
[]

[BCs]
  [axis]
    type = ADDirichletBC
    variable = disp_x
    boundary = left
    value = 0
  []
  [bottom]
    type = ADDirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
  []
  [outer_pressure]
    type = ADPressure
    variable = disp_x
    boundary = right
    factor = 1.0e10
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

[Postprocessors]
  [radial_stress]
    type = ElementAverageValue
    variable = stress_xx
  []
  [axial_stress]
    type = ElementAverageValue
    variable = stress_yy
  []
  [hoop_stress]
    type = ElementAverageValue
    variable = stress_zz
  []
[]

[Outputs]
  csv = true
[]
