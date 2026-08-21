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
    initial_condition = 400
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
        eigenstrain_names = thermal_strain
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
  [temperature_left]
    type = DirichletBC
    variable = T
    boundary = left
    value = 400
  []
  [convection_right]
    type = ADConvectiveHeatFluxBC
    variable = T
    boundary = right
    T_infinity = 300
    heat_transfer_coefficient = 100
    use_displaced_mesh = false
  []
  [fix_x]
    type = ADDirichletBC
    variable = disp_x
    boundary = 'left right bottom top back front'
    value = 0
  []
  [fix_y]
    type = ADDirichletBC
    variable = disp_y
    boundary = 'left right bottom top back front'
    value = 0
  []
  [fix_z]
    type = ADDirichletBC
    variable = disp_z
    boundary = 'left right bottom top back front'
    value = 0
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
  [thermal_expansion]
    type = ComputeThermalExpansionEigenstrain
    temperature = T
    stress_free_temperature = 300
    thermal_expansion_coeff = 1e-5
    eigenstrain_name = thermal_strain
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
  file_base = h20_02_convection
  csv = true
  exodus = true
[]
