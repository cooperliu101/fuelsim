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
  [motion]
  type = piecewise_linear
  times = 0 0.5 1
  values = 0 0.0003 0
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
  u3 = 0
  rotation_x = 0
  rotation_y = 0
  []
  [upper]
  blocks = upper
  initial_thickness = 0.1
  u3 = 0
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
    clearance_derivative = -1e6
    temperature_derivative = 2
    reference_temperature = 300
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
  [lower_y]
  type = dirichlet
  boundary = lower
  field = displacement_y
  value = 0
  []
  [upper_y]
  type = dirichlet
  boundary = upper
  field = displacement_y
  value = 1
  function = motion
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
initial_time_step = 0.125
minimum_time_step = 0.125
maximum_time_step = 0.125
growth_factor = 1
cutback_factor = 0.5
maximum_cutbacks = 0
load_ramp_time = 0
include_thermal_time_term = false
[]
[Solver]
absolute_tolerance = 1e-9
relative_tolerance = 1e-11
step_tolerance = 1e-14
maximum_iterations = 40
linear_solver = direct
direct_factorization = mumps
[]
[Outputs]
console = true
exodus = thermal_clearance_results.e
history = thermal_clearance_history.csv
checkpoint = thermal_clearance.chk
[]
