[Case]
  version = 3
  problem = transient
  geometry = cartesian_3d
[]

[Mesh]
  type = exodus
  file = ../moose/m58_integrated_hex8_30k_mesh.e
[]

[TimeFunctions]
  [power]
    type = piecewise_linear
    times = 0 0.2 0.5 0.8 1
    values = 0 0.5 1 1.2 1
  []
  [internal_pressure]
    type = piecewise_linear
    times = 0 0.2 0.5 0.8 1
    values = 0 0.3 1 0.8 1
  []
  [external_pressure]
    type = piecewise_linear
    times = 0 0.2 0.5 0.8 1
    values = 0 0.8 1 1.2 1
  []
  [axial_slide]
    type = piecewise_linear
    times = 0 0.2 0.5 0.8 1
    values = 0 6e-5 3e-4 4.8e-4 6e-4
  []
  [cladding_axial_slide]
    type = piecewise_linear
    times = 0 0.2 0.5 0.8 1
    values = 0 2e-6 1e-5 1.6e-5 2e-5
  []
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
    [eigenstrains]
      [thermal_expansion]
        function = isotropic_thermal_expansion
        thermal_expansion = 5e-6
        reference_temperature = 600
      []
    []
    [creep]
      function = norton
      coefficient = 1e-10
      reference_stress = 5e5
      stress_exponent = 3
    []
    [plasticity]
      function = linear_isotropic_hardening
      yield_stress = 1e6
      hardening_modulus = 2e10
    []
  []
[]

[Regions]
  [fuel]
    block = fuel
    material = fuel
    strain = finite
    element = c3d8rt
    initial_temperature = 600
    volumetric_heat_source = 2e7
    heat_source_function = power
  []
  [cladding]
    block = clad
    material = cladding
    strain = finite
    element = c3d8rt
    initial_temperature = 600
    volumetric_heat_source = 0
  []
[]

[Contact]
  [fuel_cladding]
    primary = clad_rmin
    secondary = fuel_outer
    [thermal]
      gap_conductivity = 0.004
      minimum_gap = 1e-6
    []
    [mechanical]
      formulation = penalty
      penalty = 1e10
      mu = 0.002
    []
  []
[]

[BoundaryConditions]
  [fuel_bottom_x]
    type = dirichlet
    boundary = fuel_bottom
    field = displacement_x
    value = 0
  []
  [fuel_bottom_y]
    type = dirichlet
    boundary = fuel_bottom
    field = displacement_y
    value = 0
  []
  [fuel_top_slide]
    type = dirichlet
    boundary = fuel_top
    field = displacement_z
    value = 1
    function = axial_slide
  []
  [cladding_bottom_x]
    type = dirichlet
    boundary = clad_bottom
    field = displacement_x
    value = 0
  []
  [cladding_bottom_y]
    type = dirichlet
    boundary = clad_bottom
    field = displacement_y
    value = 0
  []
  [cladding_bottom_z]
    type = dirichlet
    boundary = clad_bottom
    field = displacement_z
    value = 0
  []
  [cladding_top_slide]
    type = dirichlet
    boundary = clad_top
    field = displacement_z
    value = 1
    function = cladding_axial_slide
  []
  [cladding_outer_temperature]
    type = dirichlet
    boundary = clad_rmax
    field = temperature
    value = 600
  []
  [cladding_internal_pressure]
    type = pressure
    boundary = clad_rmin
    value = 5e5
    function = internal_pressure
    configuration = current
  []
  [cladding_external_pressure]
    type = pressure
    boundary = clad_rmax
    value = 2e6
    function = external_pressure
    configuration = current
  []
[]

[Executioner]
  type = transient
  end_time = 1
  initial_time_step = 0.05
  minimum_time_step = 0.00625
  maximum_time_step = 0.05
  growth_factor = 2
  cutback_factor = 0.5
  maximum_cutbacks = 3
  load_ramp_time = 0
[]

[Solver]
  # Two solver configurations were benchmarked for this 30,148-DOF case.
  # Keep the direct-MUMPS configuration below active as the default. The
  # GMRES field-split configuration is an alternative for dedicated scaling
  # experiments:
  #   linear_solver = gmres
  #   preconditioner = field_split
  # Activate exactly one configuration at a time.
  # One-step measurements use end_time=0.05 in temporary copies of this input,
  # CPUs 0 through 7, one thread per numerical library, and the default MPI
  # transport. Internal time / wall time / aggregate peak memory are:
  #   GMRES + field_split: 1 rank 84.061 / 85.07 s / 1.064 GiB;
  #                        4 ranks 21.552 / 23.34 s / 1.287 GiB;
  #                        8 ranks 16.055 / 19.23 s / 1.879 GiB.
  #   Direct MUMPS:        1 rank 83.041 / 83.97 s / 1.732 GiB;
  #                        4 ranks 42.684 / 44.54 s / 3.311 GiB;
  #                        8 ranks 34.455 / 37.68 s / 4.765 GiB.
  # Internal four-process/eight-process efficiencies are 97.51%/65.45% for
  # GMRES and 48.64%/30.13% for direct MUMPS. Wall-time efficiencies are
  # 91.12%/55.30% and 47.13%/27.86%, respectively. GMRES linear iterations
  # are 482/732/1074; direct MUMPS uses six solves per rank count.
  # Full 1 s (20-step) eight-process rerun with streaming monitoring:
  #   GMRES + field_split: stopped after 5 accepted steps at t=0.25 s,
  #                        maximum_cutbacks / DIVERGED_LINEAR_SOLVE,
  #                        181.77 s internal and 184.98 s wall time.
  #   Direct MUMPS:        completed 20/20 steps,
  #                        794.75 s internal and 797.93 s wall time.
  linear_solver = direct
  preconditioner = lu
  direct_factorization = mumps
  jacobian_lag = 2
  absolute_tolerance = 1e-7
  relative_tolerance = 1e-8
  step_tolerance = 1e-12
  maximum_iterations = 20
  backtracking_fallback = true
  field_residual_scaling = true
  temperature_residual_absolute_tolerance = 1e-8
  mechanical_residual_absolute_tolerance = 1e-4
[]

[Outputs]
  console = true
[]
