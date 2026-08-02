[Case]
  version = 1
  problem = transient
[]

[Mesh]
  type = exodus
  file = ../moose/m22_coupled_plastic_creep_rz_mesh.e
[]

[Regions]
  [material]
    block_id = 0
    conductivity_inverse_temperature = 0
    conductivity_constant = 1
    young_modulus = 2e11
    poisson_ratio = 0.3
    thermal_expansion = 0
    reference_temperature = 600
    initial_temperature = 600
    volumetric_heat_source = 0
    density = 1
    specific_heat = 1
    inelastic_model = norton_creep_j2_plasticity
    creep_coefficient = 1e-4
    creep_reference_stress = 1e8
    creep_exponent = 3
    yield_stress = 2e8
    hardening_modulus = 2e9
  []
[]

[Contact]
[]

[BoundaryConditions]
  [axis]
    type = dirichlet
    boundary = left
    field = radial_displacement
    value = 0
  []
  [bottom]
    type = dirichlet
    boundary = bottom
    field = axial_displacement
    value = 0
  []
  [top_pull]
    type = dirichlet
    boundary = top
    field = axial_displacement
    value = 2e-6
    scale_with_load = true
  []
[]

[Executioner]
  type = transient
  end_time = 1
  initial_time_step = 0.1
  minimum_time_step = 0.1
  maximum_time_step = 0.1
  growth_factor = 1
  cutback_factor = 0.5
  maximum_cutbacks = 0
  load_ramp_time = 1
[]

[Solver]
  absolute_tolerance = 1e-10
  relative_tolerance = 1e-10
  step_tolerance = 1e-12
  maximum_iterations = 40
[]

[Outputs]
  console = false
[]
