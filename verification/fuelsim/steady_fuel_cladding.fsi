[Case]
  version = 1
  problem = steady_fuel_cladding
[]

[Mesh]
  type = exodus
  file = ../moose/m1_fuel_cladding_gap_rz_mesh.e

  [fuel]
    block = fuel
    radial_inner = fuel_left
    radial_outer = fuel_right
    bottom = fuel_bottom
    top = fuel_top
  []

  [cladding]
    block = clad
    radial_inner = clad_left
    radial_outer = clad_right
    bottom = clad_bottom
    top = clad_top
  []
[]

[Materials]
  [fuel]
    conductivity_inverse_temperature = 3824
    conductivity_constant = 0.61
    young_modulus = 2e11
    poisson_ratio = 0.316
    thermal_expansion = 1e-5
    reference_temperature = 600
  []

  [cladding]
    conductivity_inverse_temperature = 0
    conductivity_constant = 16
    young_modulus = 7.5e10
    poisson_ratio = 0.3
    thermal_expansion = 5e-6
    reference_temperature = 600
  []
[]

[Physics]
  initial_temperature = 600
  outer_temperature = 600
  final_heat_source = 2e8
  gap_conductivity = 0.4
  minimum_gap = 1e-6
  contact_penalty = 1e14
[]

[Executioner]
  type = steady
  load_steps = 20
[]

[Solver]
  absolute_tolerance = 1e-8
  relative_tolerance = 1e-10
  step_tolerance = 1e-12
  maximum_iterations = 50
[]

[Outputs]
  console = true
[]
