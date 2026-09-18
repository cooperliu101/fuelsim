[Case]
version = 3
problem = transient
geometry = generalized_plane_strain
[]
[Mesh]
type = exodus
file = contact_law.e
[]
[TimeFunctions]
  [ramp]
  type = piecewise_linear
  times = 0 1
  values = 0 1
  []
  [press]
  type = piecewise_linear
  times = 0 0.25 0.5 0.75 1
  values = 0 -0.0003 -0.0005 0 -0.0007
  []
  [hot]
  type = piecewise_linear
  times = 0 1
  values = 300 400
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
    [plasticity]
    function = linear_isotropic_hardening
    yield_stress = 500
    hardening_modulus = 10000
    []
    [creep]
    function = norton
    coefficient = 0.0001
    reference_stress = 1000
    stress_exponent = 1
    []
  []
[]
[Regions]
  [lower]
  element = cpeg8t
  block = lower
  material = solid
  strain = small
  initial_temperature = 300
  volumetric_heat_source = 0
  []
  [upper]
  element = cpeg8t
  block = upper
  material = solid
  strain = small
  initial_temperature = 300
  volumetric_heat_source = 0
  []
[]
[GeneralizedPlaneStrain]
  [lower]
  blocks = lower
  initial_thickness = 0.1
  u3 = 0.003
  u3_function = ramp
  rotation_x = 0
  rotation_y = 0
  []
  [upper]
  blocks = upper
  initial_thickness = 0.1
  u3 = 0.003
  u3_function = ramp
  rotation_x = 0
  rotation_y = 0
  []
[]
[Contact]
  [interface]
  primary = primary
  secondary = secondary
    [thermal]
    law = affine
    conductance = 1000
    []
    [mechanical]
    formulation = penalty
    penalty = 1e9
    discretization = surface_to_surface
    sliding = finite
    []
  []
[]
[BoundaryConditions]
  [field_x]
  type = dirichlet
  boundary = field
  field = displacement_x
  value = 0
  []
  [bottom_y]
  type = dirichlet
  boundary = bottom
  field = displacement_y
  value = 0
  []
  [top_y]
  type = dirichlet
  boundary = top
  field = displacement_y
  value = 1
  function = press
  []
  [hot]
  type = dirichlet
  boundary = top
  field = temperature
  value = 1
  function = hot
  []
  [cold]
  type = dirichlet
  boundary = bottom
  field = temperature
  value = 300
  []
[]
[Executioner]
type = transient
end_time = 1
initial_time_step = 0.0625
minimum_time_step = 0.0625
maximum_time_step = 0.0625
growth_factor = 1
cutback_factor = 0.5
maximum_cutbacks = 0
load_ramp_time = 0
include_thermal_time_term = true
restart = coupled_contact_split.chk
[]
[Solver]
absolute_tolerance = 1e-11
relative_tolerance = 1e-12
step_tolerance = 1e-14
maximum_iterations = 40
linear_solver = direct
direct_factorization = mumps
[]
[Outputs]
console = true
exodus = coupled_contact_restart_results.e
history = coupled_contact_restart_history.csv
checkpoint = coupled_contact_restart.chk
checkpoint_interval = 1
[]
