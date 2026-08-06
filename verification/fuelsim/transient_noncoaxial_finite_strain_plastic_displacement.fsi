[Case]
  version = 1
  problem = transient
[]

[Mesh]
  type = exodus
  file = ../moose/m43_noncoaxial_finite_strain_rz_mesh.e
[]

[TimeFunctions]
  [axial_path]
    type = piecewise_linear
    times = 0 1 2 3 4 5
    values = 0 0.2 0.2 -0.1 -0.1 0.06
  []
  [shear_path]
    type = piecewise_linear
    times = 0 1 2 3 4 5
    values = 0 0 2.0 2.0 -1.5 -1.5
  []
[]

[Regions]
  [material]
    block_id = 0
    strain = finite
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
    inelastic_model = j2_plasticity
    yield_stress = 2e8
    hardening_modulus = 2e9
  []
[]

[Contact]
[]

[BoundaryConditions]
  [bottom_r]
    type = dirichlet
    boundary = bottom
    field = radial_displacement
    value = 0
  []
  [bottom_z]
    type = dirichlet
    boundary = bottom
    field = axial_displacement
    value = 0
  []
  [top_r]
    type = dirichlet
    boundary = top
    field = radial_displacement
    value = 0.001
    function = shear_path
  []
  [top_z]
    type = dirichlet
    boundary = top
    field = axial_displacement
    value = 0.001
    function = axial_path
  []
[]

[Executioner]
  type = transient
  end_time = 5
  initial_time_step = 0.05
  minimum_time_step = 0.05
  maximum_time_step = 0.05
  growth_factor = 1
  cutback_factor = 0.5
  maximum_cutbacks = 0
  load_ramp_time = 5
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
