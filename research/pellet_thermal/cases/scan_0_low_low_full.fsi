[Case]
  version = 3
  physics = thermal
  problem = steady
  geometry = cartesian_3d
[]
[Mesh]
  type = exodus
  file = concentric_gas.e
[]
[Materials]
  [fuel]
    [thermal]
      function = constant_thermophysical
      conductivity = 3
      density = 10000
      specific_heat = 300
    []
  []
  [clad]
    [thermal]
      function = constant_thermophysical
      conductivity = 16
      density = 6500
      specific_heat = 300
    []
  []
[]
[Regions]
  [pellet]
    block = pellet
    element = dc3d8
    material = fuel
    initial_temperature = 600
    volumetric_heat_source = 50000000
    pellet_response = full
  []
  [cladding]
    block = cladding
    element = dc3d8
    material = clad
    initial_temperature = 600
    volumetric_heat_source = 0
  []
[]
[Contact]
  [sector_0]
    primary = clad_inner
    secondary = pellet_sector_0
    [thermal]
      law = gas_gap
      gap_conductivity = 0.25
      minimum_gap = 0.000001
      discretization = surface_to_surface
    []
  []
  [sector_1]
    primary = clad_inner
    secondary = pellet_sector_1
    [thermal]
      law = gas_gap
      gap_conductivity = 0.25
      minimum_gap = 0.000001
      discretization = surface_to_surface
    []
  []
  [sector_2]
    primary = clad_inner
    secondary = pellet_sector_2
    [thermal]
      law = gas_gap
      gap_conductivity = 0.25
      minimum_gap = 0.000001
      discretization = surface_to_surface
    []
  []
  [sector_3]
    primary = clad_inner
    secondary = pellet_sector_3
    [thermal]
      law = gas_gap
      gap_conductivity = 0.25
      minimum_gap = 0.000001
      discretization = surface_to_surface
    []
  []
  [sector_4]
    primary = clad_inner
    secondary = pellet_sector_4
    [thermal]
      law = gas_gap
      gap_conductivity = 0.25
      minimum_gap = 0.000001
      discretization = surface_to_surface
    []
  []
  [sector_5]
    primary = clad_inner
    secondary = pellet_sector_5
    [thermal]
      law = gas_gap
      gap_conductivity = 0.25
      minimum_gap = 0.000001
      discretization = surface_to_surface
    []
  []
  [sector_6]
    primary = clad_inner
    secondary = pellet_sector_6
    [thermal]
      law = gas_gap
      gap_conductivity = 0.25
      minimum_gap = 0.000001
      discretization = surface_to_surface
    []
  []
  [sector_7]
    primary = clad_inner
    secondary = pellet_sector_7
    [thermal]
      law = gas_gap
      gap_conductivity = 0.25
      minimum_gap = 0.000001
      discretization = surface_to_surface
    []
  []
[]
[BoundaryConditions]
  [coolant]
    type = convection
    boundary = clad_outer
    heat_transfer_coefficient = 10000
    ambient_temperature = 550
  []
[]
[Executioner]
  type = steady
  load_steps = 1
[]
[Solver]
  linear_solver = direct
  direct_factorization = mumps
  absolute_tolerance = 1e-10
  relative_tolerance = 1e-10
  step_tolerance = 1e-12
  maximum_iterations = 30
[]
[Outputs]
  console = true
  exodus = scan_0_low_low_full_results.e
  csv = scan_0_low_low_full_history.csv
[]
