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
        eigenstrain_names = thermal_strain
        use_automatic_differentiation = false
      []
    []
  []
[]

[Kernels]
  [conduction]
    type = HeatConduction
    variable = T
    block = 0
    thermal_conductivity = thermal_conductivity
  []
  [capacity]
    type = HeatConductionTimeDerivative
    variable = T
    block = 0
  []
  [source]
    type = BodyForce
    variable = T
    block = 0
    value = 6e6
  []
[]

[BCs]
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
  type = Transient
  scheme = implicit-euler
  solve_type = NEWTON
  dt = 0.1
  end_time = 1
  nl_abs_tol = 1e-10
  nl_rel_tol = 1e-12
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
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
  file_base = h20_06_transient
  csv = true
  exodus = true
[]
