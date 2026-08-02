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

[AuxVariables]
  [disp_x]
    initial_condition = 0
  []
  [disp_y]
    initial_condition = 0
  []
[]

[Functions]
  [power]
    type = PiecewiseLinear
    x = '0 2.5 5 10'
    y = '0 1 0.5 1'
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
    function = power
    use_displaced_mesh = false
  []
[]

[BCs]
  [coolant]
    type = ADConvectiveHeatFluxBC
    variable = T
    boundary = right
    T_infinity = 500
    heat_transfer_coefficient = 1000
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

[Executioner]
  type = Transient
  scheme = implicit-euler
  solve_type = NEWTON
  end_time = 10
  nl_abs_tol = 1e-12
  nl_rel_tol = 1e-10
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  [TimeStepper]
    type = TimeSequenceStepper
    time_sequence = '0 2.5 5 8 10'
  []
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
