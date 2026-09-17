[Case]
  version = 3
  problem = transient
  geometry = generalized_plane_strain
[]
[Mesh]
  type = exodus
  file = rectangle.e
[]
[TimeFunctions]
  [temperature]
    type = piecewise_linear
    times = 0 1
    values = 300 400
  []
  [ramp]
    type = piecewise_linear
    times = 0 1
    values = 0 1
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
        reference_temperature = 300
      []
    []
  []
[]
[Regions]
  [solid]
    element = cpeg8t
    block = solid
    material = solid
    strain = finite
    initial_temperature = 300
    volumetric_heat_source = 0
  []
[]
[GeneralizedPlaneStrain]
  [solid]
    blocks = solid
    initial_thickness = 0.1
    u3 = 0.01
    u3_function = ramp
    rotation_x = 0.2
    rotation_x_function = ramp
    rotation_y = -0.3
    rotation_y_function = ramp
  []
[]
[BoundaryConditions]
  [fix_x]
    type = dirichlet
    boundary = field
    field = displacement_x
    value = 0
  []
  [fix_y]
    type = dirichlet
    boundary = field
    field = displacement_y
    value = 0
  []
  [temperature]
    type = dirichlet
    boundary = corners
    field = temperature
    value = 1
    function = temperature
  []
[]
[Executioner]
  type = transient
  end_time = 1
  initial_time_step = 1
  minimum_time_step = 1
  maximum_time_step = 1
  growth_factor = 1
  cutback_factor = 0.5
  maximum_cutbacks = 0
  load_ramp_time = 0
[]
[Solver]
  absolute_tolerance = 1e-10
  relative_tolerance = 1e-11
  step_tolerance = 1e-14
  maximum_iterations = 20
  linear_solver = direct
  direct_factorization = mumps
[]
[Outputs]
  console = true
  exodus = finite_prescribed_results.e
  history = finite_prescribed_history.csv
  checkpoint = finite_prescribed.chk
[]
