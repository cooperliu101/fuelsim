[Case]
 version = 3
 problem = steady
 geometry = cartesian_3d
[]
[SectionModes]
 count = 12
 [ends]
  count = 288
  lower_interface = lower_fine
  upper_interface = upper_fine
 []
 [transition]
  count = 192
  lower_interface = lower_interface
  upper_interface = upper_interface
 []
[]
[Mesh]
 type = exodus
 file = plate_local.e
[]
[Materials]
 [al]
  [thermal]
   function = constant_thermophysical
   conductivity = 1
   density = 2700
   specific_heat = 900
  []
  [elasticity]
   function = constant_isotropic
   young_modulus = 7e10
   poisson_ratio = 0.3
  []
 []
[]
[Regions]
 [plate]
  element = c3d20t
  block = plate
  material = al
  strain = small
  initial_temperature = 300
  volumetric_heat_source = 0
 []
[]
[BoundaryConditions]
 [root_x]
  type = dirichlet
  boundary = root
  field = displacement_x
  value = 0
 []
 [root_y]
  type = dirichlet
  boundary = root
  field = displacement_y
  value = 0
 []
 [root_z]
  type = dirichlet
  boundary = root
  field = displacement_z
  value = 0
 []
 [temperature]
  type = dirichlet
  boundary = all_nodes
  field = temperature
  value = 300
 []
 [couple_plus]
  type = traction
  boundary = end_plus
  field = displacement_z
  value = 1e4
  configuration = reference
 []
 [couple_minus]
  type = traction
  boundary = end_minus
  field = displacement_z
  value = -1e4
  configuration = reference
 []
[]
[Executioner]
 type = steady
 load_steps = 1
[]
[Solver]
 linear_solver = direct
 preconditioner = lu
 direct_factorization = mumps
 absolute_tolerance = 1e-8
 relative_tolerance = 1e-10
 maximum_iterations = 10
[]
[Outputs]
 console = true
 csv = bending_local_12_summary.csv
 history = bending_local_12_history.csv
[]
