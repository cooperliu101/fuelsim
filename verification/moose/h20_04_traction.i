[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 1
  ny = 1
  nz = 1
  elem_type = HEX20
  xmax = 1
  ymax = 1
  zmax = 1
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Variables]
  [T]
    initial_condition = 300
  []
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [solid]
        block = 0
        strain = SMALL
        add_variables = true
        temperature = T
        use_automatic_differentiation = false
      []
    []
  []
[]

[Kernels]
  [heat_conduction]
    type = HeatConduction
    variable = T
    block = 0
    thermal_conductivity = thermal_conductivity
  []
[]

[BCs]
  [temperature]
    type = DirichletBC
    variable = T
    boundary = left
    value = 300
  []
  [fix_x]
    type = ADDirichletBC
    variable = disp_x
    boundary = left
    value = 0
  []
  [fix_y]
    type = ADDirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
  []
  [fix_z]
    type = ADDirichletBC
    variable = disp_z
    boundary = back
    value = 0
  []
  [traction]
    type = ADFunctionNeumannBC
    variable = disp_x
    boundary = right
    function = 1e6
    use_displaced_mesh = false
  []
[]

[Materials]
  [thermal]
    type = GenericConstantMaterial
    prop_names = 'thermal_conductivity density specific_heat'
    prop_values = '10 6000 1000'
    block = 0
  []
  [elasticity]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 1e9
    poissons_ratio = 0.25
    block = 0
  []
  [stress]
    type = ComputeLinearElasticStress
    block = 0
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
  [Quadrature]
    order = THIRD
  []
[]

[VectorPostprocessors]
  [temperature_nodes]
    type = NodalValueSampler
    variable = T
    sort_by = id
    use_displaced_mesh = false
  []
  [displacement_nodes]
    type = NodalValueSampler
    variable = 'disp_x disp_y disp_z'
    sort_by = id
    use_displaced_mesh = false
  []
[]

[Outputs]
  file_base = h20_04_traction
  csv = true
  exodus = true
[]
