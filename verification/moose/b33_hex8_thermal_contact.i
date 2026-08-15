[Mesh]
  file = b33_hex8_contact_mesh.e
[]

[Variables]
  [T]
  []
[]

[ICs]
  [primary_temperature]
    type = ConstantIC
    variable = T
    block = primary
    value = 300
  []
  [secondary_temperature]
    type = ConstantIC
    variable = T
    block = secondary
    value = 400
  []
[]

[Kernels]
  [primary_heat]
    type = ADMatDiffusion
    variable = T
    block = primary
    diffusivity = thermal_conductivity
  []
  [secondary_heat]
    type = ADMatDiffusion
    variable = T
    block = secondary
    diffusivity = thermal_conductivity
  []
[]

[ThermalContact]
  [interface]
    type = GapHeatTransfer
    variable = T
    primary = primary_right
    secondary = secondary_left
    gap_conductivity = 0.2
    quadrature = true
    min_gap = 1e-5
    min_gap_order = 0
    max_gap = 1e6
    emissivity_primary = 0
    emissivity_secondary = 0
  []
[]

[BCs]
  [primary_temperature]
    type = DirichletBC
    variable = T
    boundary = primary_left
    value = 300
  []
  [secondary_temperature]
    type = DirichletBC
    variable = T
    boundary = secondary_right
    value = 400
  []
[]

[Materials]
  [primary_conductivity]
    type = ADGenericConstantMaterial
    prop_names = thermal_conductivity
    prop_values = 10
    block = primary
  []
  [secondary_conductivity]
    type = ADGenericConstantMaterial
    prop_names = thermal_conductivity
    prop_values = 10
    block = secondary
  []
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  nl_abs_tol = 1e-10
  nl_rel_tol = 1e-12
  nl_max_its = 20
  automatic_scaling = true
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[VectorPostprocessors]
  [all_nodes]
    type = NodalValueSampler
    block = 'primary secondary'
    variable = T
    sort_by = id
    use_displaced_mesh = false
  []
[]

[Outputs]
  file_base = b33_hex8_thermal_contact
  csv = true
[]
