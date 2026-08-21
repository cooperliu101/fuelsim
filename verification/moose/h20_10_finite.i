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
  [T] initial_condition = 300 []
[]
[Physics/SolidMechanics/QuasiStatic]
  [solid]
    strain = FINITE
    incremental = true
    add_variables = true
    temperature = T
    eigenstrain_names = thermal_strain
    use_automatic_differentiation = false
  []
[]
[Kernels]
  [heat] type = MatDiffusion variable = T diffusivity = thermal_conductivity []
[]
[BCs]
  [temperature] type = DirichletBC variable = T boundary = left value = 300 []
  [fix_x] type = DirichletBC variable = disp_x boundary = left value = 0 []
  [fix_y] type = DirichletBC variable = disp_y boundary = bottom value = 0 []
  [fix_z] type = DirichletBC variable = disp_z boundary = back value = 0 []
  [traction] type = FunctionNeumannBC variable = disp_x boundary = right function = 1e6 use_displaced_mesh = true []
[]
[Materials]
  [thermal] type = GenericConstantMaterial prop_names = 'thermal_conductivity density specific_heat' prop_values = '10 6000 1000' []
  [elasticity] type = ComputeIsotropicElasticityTensor youngs_modulus = 1e9 poissons_ratio = 0.25 []
  [stress] type = ComputeFiniteStrainElasticStress []
  [thermal_expansion]
    type = ComputeThermalExpansionEigenstrain
    temperature = T
    stress_free_temperature = 300
    thermal_expansion_coeff = 1e-5
    eigenstrain_name = thermal_strain
  []
[]
[Executioner]
  type = Steady
  solve_type = NEWTON
  nl_abs_tol = 1e-10
  nl_rel_tol = 1e-12
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]
[VectorPostprocessors]
  [temperature_nodes] type = NodalValueSampler variable = T sort_by = id use_displaced_mesh = false []
  [displacement_nodes] type = NodalValueSampler variable = 'disp_x disp_y disp_z' sort_by = id use_displaced_mesh = false []
[]
[Outputs]
  file_base = h20_10_finite
  csv = true
  exodus = true
[]
