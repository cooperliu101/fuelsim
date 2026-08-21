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
  [disp_x]
    order = SECOND
  []
  [disp_y]
    order = SECOND
  []
  [disp_z]
    order = SECOND
  []
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
  [temperature]
    type = DirichletBC
    variable = T
    boundary = 'left right bottom top back front'
    value = 400
  []
  [fix_x]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0
  []
  [fix_y]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
  []
  [fix_z]
    type = DirichletBC
    variable = disp_z
    boundary = back
    value = 0
  []
[]

[Materials]
  [conductivity]
    type = GenericConstantMaterial
    prop_names = thermal_conductivity
    prop_values = 10
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
  file_base = b6_hex20_u2_t1
  csv = true
  exodus = true
[]
