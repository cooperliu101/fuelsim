[Mesh]
  [mesh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 40
    ny = 10
    xmin = 0.0
    xmax = 0.00412
    ymin = 0.0
    ymax = 0.010
    boundary_name_prefix = fuel
  []
  coord_type = RZ
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Variables]
  [T]
    initial_condition = 600.0
  []
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [fuel]
        strain = SMALL
        add_variables = true
        temperature = T
        eigenstrain_names = thermal_strain
      []
    []
  []
[]

[Kernels]
  [heat_conduction]
    type = HeatConduction
    variable = T
    use_displaced_mesh = false
  []
  [heat_source]
    type = BodyForce
    variable = T
    value = 2.0e8
    use_displaced_mesh = false
  []
[]

[BCs]
  [temperature_outer]
    type = DirichletBC
    variable = T
    boundary = fuel_right
    value = 600.0
  []
  [radial_axis]
    type = DirichletBC
    variable = disp_x
    boundary = fuel_left
    value = 0.0
  []
  [axial_bottom]
    type = DirichletBC
    variable = disp_y
    boundary = fuel_bottom
    value = 0.0
  []
[]

[Materials]
  [conductivity]
    type = ParsedMaterial
    property_name = thermal_conductivity
    expression = '3824.0 / T + 0.61'
    coupled_variables = T
  []
  [elasticity]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 2.0e11
    poissons_ratio = 0.316
  []
  [thermal_expansion]
    type = ComputeThermalExpansionEigenstrain
    temperature = T
    stress_free_temperature = 600.0
    thermal_expansion_coeff = 10.0e-6
    eigenstrain_name = thermal_strain
  []
  [stress]
    type = ComputeLinearElasticStress
  []
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  nl_abs_tol = 1.0e-10
  nl_rel_tol = 1.0e-10
  nl_max_its = 40
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Postprocessors]
  [temperature_center]
    type = PointValue
    point = '0 0.005 0'
    variable = T
  []
  [radial_displacement_outer_mid]
    type = PointValue
    point = '0.00412 0.005 0'
    variable = disp_x
  []
  [radial_displacement_outer_top]
    type = PointValue
    point = '0.00412 0.010 0'
    variable = disp_x
  []
  [axial_displacement_axis_top]
    type = PointValue
    point = '0 0.010 0'
    variable = disp_y
  []
[]

[Outputs]
  csv = true
[]
