[Case]
  version = 1
  problem = transient
[]

[Mesh]
  type = exodus
  file = ../moose/m41_finite_strain_pcmi_rz_mesh.e
[]

[Regions]
  [fuel]
    block = fuel
    strain = finite
    conductivity_inverse_temperature = 3824
    conductivity_constant = 0.61
    young_modulus = 2e11
    poisson_ratio = 0.316
    thermal_expansion = 1e-5
    reference_temperature = 600
    initial_temperature = 600
    volumetric_heat_source = 2e8
    density = 10970
    specific_heat = 300
    inelastic_model = elastic
  []

  [cladding]
    block = clad
    strain = finite
    conductivity_inverse_temperature = 0
    conductivity_constant = 16
    young_modulus = 7.5e10
    poisson_ratio = 0.3
    thermal_expansion = 0
    reference_temperature = 600
    initial_temperature = 600
    volumetric_heat_source = 0
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

[Contact]
  [fuel_cladding]
    primary = clad_left
    secondary = fuel_right

    [thermal]
      gap_conductivity = 0.4
      minimum_gap = 1e-6
    []

    [mechanical]
      formulation = penalty
      penalty = 1e14
    []
  []
[]

[BoundaryConditions]
  [fuel_axis]
    type = dirichlet
    boundary = fuel_left
    field = radial_displacement
    value = 0
  []

  [fuel_bottom]
    type = dirichlet
    boundary = fuel_bottom
    field = axial_displacement
    value = 0
  []

  [cladding_bottom]
    type = dirichlet
    boundary = clad_bottom
    field = axial_displacement
    value = 0
  []

  [cladding_outer_temperature]
    type = dirichlet
    boundary = clad_right
    field = temperature
    value = 600
  []
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
  load_ramp_time = 20
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
