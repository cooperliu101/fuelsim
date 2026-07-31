[Case]
  version = 1
  problem = transient_fuel_cladding
[]

[Mesh]
  type = exodus
  file = ../moose/m23_pcmi_coupled_cladding_rz_mesh.e

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
    density = 10970
    specific_heat = 300
    inelastic_model = elastic
  []

  [cladding]
    conductivity_inverse_temperature = 0
    conductivity_constant = 16
    young_modulus = 7.5e10
    poisson_ratio = 0.3
    thermal_expansion = 0
    reference_temperature = 600
    density = 6500
    specific_heat = 330
    inelastic_model = norton_creep_j2_plasticity
    creep_coefficient = 1e-5
    creep_reference_stress = 5e6
    creep_exponent = 3
    yield_stress = 5e6
    hardening_modulus = 2e9
  []
[]

[Physics]
  initial_temperature = 600
  outer_temperature = 600
  final_heat_source = 2e8
  heat_source_ramp_time = 20
  gap_conductivity = 0.4
  minimum_gap = 1e-6
  contact_penalty = 1e14
[]

[Executioner]
  type = transient
  end_time = 20
  initial_time_step = 1
  minimum_time_step = 0.125
  maximum_time_step = 1
  growth_factor = 1
  cutback_factor = 0.5
  maximum_cutbacks = 3
[]

[Solver]
  absolute_tolerance = 1e-8
  relative_tolerance = 1e-10
  step_tolerance = 1e-12
  maximum_iterations = 80
[]

[Outputs]
  console = true
[]
