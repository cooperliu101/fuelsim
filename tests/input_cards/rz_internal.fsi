# Shared model for RZ internal contracts; also runnable with fuelsim -i.
[Case]
  version = 3
  problem = transient
  geometry = axisymmetric_rz
[]

[Mesh]
  type = exodus
  file = ../../verification/moose/m23_pcmi_coupled_cladding_rz_mesh.e
[]

[Materials]
  [fuel]
    [thermal]
      function = inverse_temperature_thermophysical
      conductivity_inverse_temperature = 3824
      conductivity_constant = 0.61
      density = 10970
      specific_heat = 300
    []
    [elasticity]
      function = constant_isotropic
      young_modulus = 2e11
      poisson_ratio = 0.316
    []
    [eigenstrains]
      [thermal_expansion]
        function = isotropic_thermal_expansion
        thermal_expansion = 1e-5
        reference_temperature = 600
      []
    []
  []
  [cladding]
    [thermal]
      function = constant_thermophysical
      conductivity = 16
      density = 6500
      specific_heat = 330
    []
    [elasticity]
      function = constant_isotropic
      young_modulus = 7.5e10
      poisson_ratio = 0.3
    []
    [creep]
      function = norton
      coefficient = 1e-5
      reference_stress = 5e6
      stress_exponent = 3
    []
    [plasticity]
      function = linear_isotropic_hardening
      yield_stress = 4e6
      hardening_modulus = 2e9
    []
  []
[]

[Regions]
  [fuel]
    element = cax4t
    block = fuel
    material = fuel
    strain = small
    initial_temperature = 600
    volumetric_heat_source = 2e8
  []
  [cladding]
    element = cax4t
    block = clad
    material = cladding
    strain = small
    initial_temperature = 600
    volumetric_heat_source = 0
  []
[]

[Contact]
  [fuel_cladding]
    primary = clad_left
    secondary = fuel_right
    [thermal]
      gap_conductivity = 0.4
      minimum_gap = 2e-6
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
