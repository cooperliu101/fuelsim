[Case]
 version = 3
 problem = steady
 geometry = cartesian_3d
[]
[SectionModes]
 count = 12
 width_lines = width_0 width_1 width_2 width_3 width_4
[]
[Mesh]
 type = exodus
 file = plate_width.e
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
 [support_a_x]
  type = dirichlet
  boundary = support_a
  field = displacement_x
  value = 0
 []
 [support_a_y]
  type = dirichlet
  boundary = support_a
  field = displacement_y
  value = 0
 []
 [support_a_z]
  type = dirichlet
  boundary = support_a
  field = displacement_z
  value = 0
 []
 [support_b_y]
  type = dirichlet
  boundary = support_b
  field = displacement_y
  value = 0
 []
 [support_b_z]
  type = dirichlet
  boundary = support_b
  field = displacement_z
  value = 0
 []
 [support_c_z]
  type = dirichlet
  boundary = support_c
  field = displacement_z
  value = 0
 []
 [temperature]
  type = dirichlet
  boundary = all_nodes
  field = temperature
  value = 300
 []
 [pressure]
  type = pressure
  boundary = broad
  value = 10
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
 csv = transverse_width4_64_summary.csv
 history = transverse_width4_64_history.csv
[]
