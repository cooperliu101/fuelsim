[Mesh]
  file = b35_hex8_shared_meat_clad_mesh.e
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
      [meat]
        block = meat
        strain = SMALL
        add_variables = true
        temperature = T
        eigenstrain_names = meat_thermal_strain
        use_automatic_differentiation = true
      []
      [clad]
        block = clad
        strain = SMALL
        add_variables = true
        temperature = T
        eigenstrain_names = clad_thermal_strain
        use_automatic_differentiation = true
      []
    []
  []
[]

[Kernels]
  [meat_heat]
    type = ADMatDiffusion
    variable = T
    block = meat
    diffusivity = thermal_conductivity
  []
  [clad_heat]
    type = ADMatDiffusion
    variable = T
    block = clad
    diffusivity = thermal_conductivity
  []
[]

[BCs]
  [temperature_left]
    type = DirichletBC
    variable = T
    boundary = plate_left
    value = 300
  []
  [temperature_right]
    type = DirichletBC
    variable = T
    boundary = plate_right
    value = 600
  []
  [fix_x]
    type = ADDirichletBC
    variable = disp_x
    boundary = plate_left
    value = 0
  []
  [fix_y]
    type = ADDirichletBC
    variable = disp_y
    boundary = plate_left
    value = 0
  []
  [fix_z]
    type = ADDirichletBC
    variable = disp_z
    boundary = plate_left
    value = 0
  []
[]

[Materials]
  [meat_conductivity]
    type = ADGenericConstantMaterial
    prop_names = thermal_conductivity
    prop_values = 10
    block = meat
  []
  [clad_conductivity]
    type = ADGenericConstantMaterial
    prop_names = thermal_conductivity
    prop_values = 20
    block = clad
  []
  [meat_elasticity]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 1e9
    poissons_ratio = 0.25
    block = meat
  []
  [clad_elasticity]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 5e9
    poissons_ratio = 0.3
    block = clad
  []
  [meat_thermal_expansion]
    type = ADComputeThermalExpansionEigenstrain
    temperature = T
    stress_free_temperature = 300
    thermal_expansion_coeff = 1e-5
    eigenstrain_name = meat_thermal_strain
    block = meat
  []
  [clad_thermal_expansion]
    type = ADComputeThermalExpansionEigenstrain
    temperature = T
    stress_free_temperature = 300
    thermal_expansion_coeff = 5e-6
    eigenstrain_name = clad_thermal_strain
    block = clad
  []
  [meat_stress]
    type = ADComputeLinearElasticStress
    block = meat
  []
  [clad_stress]
    type = ADComputeLinearElasticStress
    block = clad
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
    variable = 'T disp_x disp_y disp_z'
    sort_by = id
    use_displaced_mesh = false
  []
[]

[Outputs]
  file_base = b35_hex8_shared_meat_clad
  csv = true
[]
