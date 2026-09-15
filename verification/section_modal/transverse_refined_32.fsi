[Case]
 version = 3
 problem = steady
 geometry = cartesian_3d
[]
[SectionModes]
 count = 32
[]
[Mesh]
 type = exodus
 file = plate_refined.e
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
 [transverse]
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
 csv = transverse_refined_32_summary.csv
 history = transverse_refined_32_history.csv
[]
