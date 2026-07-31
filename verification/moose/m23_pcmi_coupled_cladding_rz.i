# M2.3 PCMI verification: elastic fuel thermal expansion closes the gap and
# loads a cladding with coupled Norton creep and J2 plasticity.

[Mesh]
  [fuel_mesh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 6
    ny = 4
    xmin = 0.0
    xmax = 0.00412
    ymin = 0.0
    ymax = 0.010
    boundary_name_prefix = fuel
  []
  [clad_mesh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 2
    ny = 4
    xmin = 0.004121
    xmax = 0.004692
    ymin = 0.0
    ymax = 0.010020
    boundary_name_prefix = clad
    boundary_id_offset = 10
  []
  [clad_id]
    type = RenameBlockGenerator
    input = clad_mesh
    old_block = '0'
    new_block = '1'
  []
  [combined]
    type = MeshCollectionGenerator
    inputs = 'fuel_mesh clad_id'
  []
  [rename_block]
    type = RenameBlockGenerator
    input = combined
    old_block = '0 1'
    new_block = 'fuel clad'
  []
  coord_type = RZ
  patch_update_strategy = iteration
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Variables]
  [T]
    initial_condition = 600.0
  []
[]

[Physics/SolidMechanics/QuasiStatic]
  [fuel]
    block = fuel
    strain = SMALL
    incremental = true
    add_variables = true
    temperature = T
    eigenstrain_names = 'fuel_thermal_strain'
    use_automatic_differentiation = true
  []
  [clad]
    block = clad
    strain = SMALL
    incremental = true
    add_variables = true
    temperature = T
    eigenstrain_names = 'clad_thermal_strain'
    use_automatic_differentiation = true
    generate_output = 'vonmises_stress'
  []
[]

[Functions]
  [heat_ramp]
    type = ParsedFunction
    expression = 't'
  []
[]

[Kernels]
  [heat_time]
    type = HeatConductionTimeDerivative
    variable = T
    block = 'fuel clad'
    use_displaced_mesh = false
  []
  [heat_conduction]
    type = HeatConduction
    variable = T
    block = 'fuel clad'
    use_displaced_mesh = false
  []
  [heat_source_fuel]
    type = BodyForce
    variable = T
    block = fuel
    value = 1.0e7
    function = heat_ramp
    use_displaced_mesh = false
  []
[]

[ThermalContact]
  [gap]
    type = GapHeatTransfer
    variable = T
    primary = clad_left
    secondary = fuel_right
    gap_conductivity = 0.4
    quadrature = true
    gap_geometry_type = CYLINDER
    min_gap = 1.0e-6
    min_gap_order = 0
    max_gap = 1.0e6
    emissivity_primary = 0
    emissivity_secondary = 0
  []
[]

[Contact]
  [mechanical]
    primary = clad_left
    secondary = fuel_right
    formulation = penalty
    model = frictionless
    penalty = 1.0e14
    normalize_penalty = true
  []
[]

[BCs]
  [u_center]
    type = DirichletBC
    variable = disp_x
    boundary = fuel_left
    value = 0.0
  []
  [u_bottom_f]
    type = DirichletBC
    variable = disp_y
    boundary = fuel_bottom
    value = 0.0
  []
  [u_bottom_c]
    type = DirichletBC
    variable = disp_y
    boundary = clad_bottom
    value = 0.0
  []
  [T_outer]
    type = DirichletBC
    variable = T
    boundary = clad_right
    value = 600.0
  []
[]

[Materials]
  [fuel_tk]
    type = ParsedMaterial
    block = fuel
    property_name = thermal_conductivity
    expression = '3824.0 / T + 0.61'
    coupled_variables = T
  []
  [fuel_thermal_mass]
    type = GenericConstantMaterial
    block = fuel
    prop_names = 'density specific_heat'
    prop_values = '10970.0 300.0'
  []
  [fuel_elasticity]
    type = ADComputeIsotropicElasticityTensor
    block = fuel
    youngs_modulus = 2.0e11
    poissons_ratio = 0.316
  []
  [fuel_expansion]
    type = ADComputeThermalExpansionEigenstrain
    block = fuel
    temperature = T
    stress_free_temperature = 600.0
    thermal_expansion_coeff = 10.0e-6
    eigenstrain_name = fuel_thermal_strain
  []
  [fuel_stress]
    type = ADComputeMultipleInelasticStress
    block = fuel
    inelastic_models = ''
    perform_finite_strain_rotations = false
  []

  [clad_thermal]
    type = GenericConstantMaterial
    block = clad
    prop_names = 'thermal_conductivity density specific_heat'
    prop_values = '16.0 6500.0 330.0'
  []
  [clad_elasticity]
    type = ADComputeIsotropicElasticityTensor
    block = clad
    youngs_modulus = 75.0e9
    poissons_ratio = 0.3
  []
  [clad_expansion]
    type = ADComputeThermalExpansionEigenstrain
    block = clad
    temperature = T
    stress_free_temperature = 600.0
    thermal_expansion_coeff = 0.0
    eigenstrain_name = clad_thermal_strain
  []
  [clad_stress]
    type = ADComputeMultipleInelasticStress
    block = clad
    inelastic_models = 'clad_creep clad_plasticity'
    perform_finite_strain_rotations = false
    max_iterations = 100
    relative_tolerance = 1.0e-12
    absolute_tolerance = 1.0e-6
  []
  [clad_creep]
    type = ADPowerLawCreepStressUpdate
    block = clad
    coefficient = 8.0e-26
    n_exponent = 3.0
    m_exponent = 0.0
    activation_energy = 0.0
  []
  [clad_plasticity]
    type = ADIsotropicPlasticityStressUpdate
    block = clad
    yield_stress = 5.0e6
    hardening_constant = 2.0e9
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  line_search = basic
  dt = 1.0
  end_time = 20.0
  nl_max_its = 80
  nl_abs_tol = 1.0e-8
  nl_rel_tol = 1.0e-8
  abort_on_solve_fail = true
  automatic_scaling = true
  compute_scaling_once = false
  off_diagonals_in_auto_scaling = true
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Postprocessors]
  [T_center]
    type = PointValue
    point = '0 0.005 0'
    variable = T
    use_displaced_mesh = false
  []
  [T_fuel_surface]
    type = PointValue
    point = '0.00412 0.005 0'
    variable = T
    use_displaced_mesh = false
  []
  [T_clad_inner]
    type = PointValue
    point = '0.004121 0.00501 0'
    variable = T
    use_displaced_mesh = false
  []
  [u_fuel_outer]
    type = PointValue
    point = '0.00412 0.005 0'
    variable = disp_x
    use_displaced_mesh = false
  []
  [u_clad_inner]
    type = PointValue
    point = '0.004121 0.00501 0'
    variable = disp_x
    use_displaced_mesh = false
  []
  [u_clad_outer]
    type = PointValue
    point = '0.004692 0.00501 0'
    variable = disp_x
    use_displaced_mesh = false
  []
  [u_z_fuel_top]
    type = PointValue
    point = '0 0.010 0'
    variable = disp_y
    use_displaced_mesh = false
  []
  [average_effective_plastic]
    type = ADElementAverageMaterialProperty
    mat_prop = effective_plastic_strain
    block = clad
  []
  [average_effective_creep]
    type = ADElementAverageMaterialProperty
    mat_prop = effective_creep_strain
    block = clad
  []
  [average_vonmises_stress]
    type = ElementAverageValue
    variable = vonmises_stress
    block = clad
  []
[]

[VectorPostprocessors]
  [fuel_surface]
    type = NodalValueSampler
    boundary = fuel_right
    variable = 'T disp_x disp_y contact_pressure nodal_area penetration'
    sort_by = y
    use_displaced_mesh = false
  []
  [clad_inner]
    type = NodalValueSampler
    boundary = clad_left
    variable = 'T disp_x disp_y'
    sort_by = y
    use_displaced_mesh = false
  []
[]

[Outputs]
  csv = true
[]
