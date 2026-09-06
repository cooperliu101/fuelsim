[Case]
  version = 3
  problem = transient
  geometry = axisymmetric_rz
[]
[Mesh]
  type = exodus
  file = ../meshes/b10_rz_probe.e
[]
[TimeFunctions]
  [hot]
    type = piecewise_linear
    times = 0 0.1
    values = 600 620
  []
  [warm]
    type = piecewise_linear
    times = 0 0.1
    values = 600 610
  []
  [cold]
    type = piecewise_linear
    times = 0 0.1
    values = 600 590
  []
[]
[Materials]
  [solid]
    [thermal]
      function = constant_thermophysical
      conductivity = 10
      density = 1000
      specific_heat = 100
    []
    [elasticity]
      function = constant_isotropic
      young_modulus = 1e6
      poisson_ratio = 0.25
    []
    [eigenstrains]
      [expansion]
        function = isotropic_thermal_expansion
        thermal_expansion = 1e-5
        reference_temperature = 600
      []
    []
  []
[]
[Regions]
  [solid]
    block = solid
    material = solid
    element = cax8rt
    strain = finite
    initial_temperature = 600
    volumetric_heat_source = 1e6
  []
[]
[BoundaryConditions]
  [r1]
    type = dirichlet
    boundary = n1
    field = radial_displacement
    value = 0.01
    scale_with_load = true
  []
  [z1]
    type = dirichlet
    boundary = n1
    field = axial_displacement
    value = 0.02
    scale_with_load = true
  []
  [r2]
    type = dirichlet
    boundary = n2
    field = radial_displacement
    value = 0.03
    scale_with_load = true
  []
  [z2]
    type = dirichlet
    boundary = n2
    field = axial_displacement
    value = -0.01
    scale_with_load = true
  []
  [r3]
    type = dirichlet
    boundary = n3
    field = radial_displacement
    value = 0.06
    scale_with_load = true
  []
  [z3]
    type = dirichlet
    boundary = n3
    field = axial_displacement
    value = 0.07
    scale_with_load = true
  []
  [r4]
    type = dirichlet
    boundary = n4
    field = radial_displacement
    value = -0.02
    scale_with_load = true
  []
  [z4]
    type = dirichlet
    boundary = n4
    field = axial_displacement
    value = 0.05
    scale_with_load = true
  []
  [r5]
    type = dirichlet
    boundary = n5
    field = radial_displacement
    value = 0.025
    scale_with_load = true
  []
  [z5]
    type = dirichlet
    boundary = n5
    field = axial_displacement
    value = 0.003
    scale_with_load = true
  []
  [r6]
    type = dirichlet
    boundary = n6
    field = radial_displacement
    value = 0.038
    scale_with_load = true
  []
  [z6]
    type = dirichlet
    boundary = n6
    field = axial_displacement
    value = 0.034
    scale_with_load = true
  []
  [r7]
    type = dirichlet
    boundary = n7
    field = radial_displacement
    value = 0.022
    scale_with_load = true
  []
  [z7]
    type = dirichlet
    boundary = n7
    field = axial_displacement
    value = 0.065
    scale_with_load = true
  []
  [r8]
    type = dirichlet
    boundary = n8
    field = radial_displacement
    value = -0.002
    scale_with_load = true
  []
  [z8]
    type = dirichlet
    boundary = n8
    field = axial_displacement
    value = 0.032
    scale_with_load = true
  []
  [t1]
    type = dirichlet
    boundary = n1
    field = temperature
    value = 600
  []
  [t2]
    type = dirichlet
    boundary = n2
    field = temperature
    value = 1
    function = hot
  []
  [t3]
    type = dirichlet
    boundary = n3
    field = temperature
    value = 1
    function = cold
  []
  [t4]
    type = dirichlet
    boundary = n4
    field = temperature
    value = 1
    function = warm
  []
[]
[Executioner]
  type = transient
  end_time = 0.1
  initial_time_step = 0.1
  minimum_time_step = 0.1
  maximum_time_step = 0.1
  growth_factor = 1
  cutback_factor = 0.5
  maximum_cutbacks = 0
  load_ramp_time = 0
[]
[Solver]
  linear_solver = direct
  direct_factorization = mumps
  absolute_tolerance = 1e-10
  relative_tolerance = 1e-10
  maximum_iterations = 40
[]
[Outputs]
  console = true
  csv = transient_b121_cax8rt_finite_probe_summary.csv
  exodus = transient_b121_cax8rt_finite_probe_results.e
[]
