[Case]
version = 3
problem = steady
geometry = generalized_plane_strain
[]
[Mesh]
type = exodus
file = contact_refine_4.e
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
    [mechanical]
    formulation = penalty
    penalty = 1e9
    discretization = surface_to_surface
    sliding = finite
    []
  []
[]
[BoundaryConditions]
  [lower_x]
  type = dirichlet
  boundary = lower
  field = displacement_x
  value = 0
  []
  [lower_y]
  type = dirichlet
  boundary = lower
  field = displacement_y
  value = 0
  []
  [top_x]
  type = dirichlet
  boundary = top
  field = displacement_x
  value = 0
  []
  [top_y]
  type = dirichlet
  boundary = top
  field = displacement_y
  value = -0.0005
  []
  [temperature]
  type = dirichlet
  boundary = field
  field = temperature
  value = 300
  []
[]
[Executioner]
type = steady
load_steps = 4
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
exodus = contact_refine_4_results.e
[]
