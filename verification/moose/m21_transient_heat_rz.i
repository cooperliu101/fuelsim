[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 4
  ny = 2
  xmax = 0.004
  ymax = 0.010
  coord_type = RZ
[]

[Variables]
  [T]
    initial_condition = 600
  []
[]

[Kernels]
  [conduction]
    type = ADHeatConduction
    variable = T
    use_displaced_mesh = false
  []
  [capacity]
    type = ADHeatConductionTimeDerivative
    variable = T
    use_displaced_mesh = false
  []
  [source]
    type = ADBodyForce
    variable = T
    value = 3e6
    use_displaced_mesh = false
  []
[]

[Materials]
  [thermal]
    type = ADGenericConstantMaterial
    prop_names = 'thermal_conductivity specific_heat density'
    prop_values = '2 300 10000'
  []
[]

[Postprocessors]
  [T_average]
    type = ElementAverageValue
    variable = T
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [T_error]
    type = NodalL2Error
    variable = T
    function = '600+t'
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[Executioner]
  type = Transient
  scheme = implicit-euler
  solve_type = NEWTON
  dt = 1
  end_time = 10
  nl_abs_tol = 1e-12
  nl_rel_tol = 1e-10
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  csv = true
[]
