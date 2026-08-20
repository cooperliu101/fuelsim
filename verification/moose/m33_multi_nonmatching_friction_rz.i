# M3.3 nonmatching frictional multi-contact verification: the three independent
# blocks use different axial subdivisions on the two radial contact interfaces.
# The same generated Exodus mesh is read by fuelsim for the comparison.

[Mesh]
  [pellet_mesh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 1
    ny = 2
    xmin = 0
    xmax = 0.001
    ymin = 0
    ymax = 0.001
    boundary_name_prefix = pellet
  []
  [inner_clad_mesh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 1
    ny = 3
    xmin = 0.001
    xmax = 0.00101
    ymin = 0
    ymax = 0.001
    boundary_name_prefix = inner_clad
    boundary_id_offset = 20
  []
  [outer_clad_mesh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 1
    ny = 4
    xmin = 0.00101
    xmax = 0.00102
    ymin = 0
    ymax = 0.001
    boundary_name_prefix = outer_clad
    boundary_id_offset = 40
  []
  [pellet_block]
    type = RenameBlockGenerator
    input = pellet_mesh
    old_block = 0
    new_block = 1
  []
  [inner_clad_block]
    type = RenameBlockGenerator
    input = inner_clad_mesh
    old_block = 0
    new_block = 2
  []
  [outer_clad_block]
    type = RenameBlockGenerator
    input = outer_clad_mesh
    old_block = 0
    new_block = 3
  []
  [combined]
    type = MeshCollectionGenerator
    inputs = 'pellet_block inner_clad_block outer_clad_block'
  []
  [blocks]
    type = RenameBlockGenerator
    input = combined
    old_block = '1 2 3'
    new_block = 'pellet inner_clad outer_clad'
  []
  coord_type = RZ
  patch_update_strategy = iteration
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[Variables]
  [T]
    initial_condition = 300
  []
[]

[AuxVariables]
  [pellet_tangential_force_y]
    family = LAGRANGE
    order = FIRST
  []
  [inner_clad_tangential_force_y]
    family = LAGRANGE
    order = FIRST
  []
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [pellet]
        block = pellet
        strain = SMALL
        incremental = true
        add_variables = true
        temperature = T
        eigenstrain_names = pellet_thermal_strain
      []
      [inner_clad]
        block = inner_clad
        strain = SMALL
        incremental = true
        add_variables = true
        temperature = T
        eigenstrain_names = inner_clad_thermal_strain
      []
      [outer_clad]
        block = outer_clad
        strain = SMALL
        incremental = true
        add_variables = true
        temperature = T
        eigenstrain_names = outer_clad_thermal_strain
      []
    []
  []
[]

[Kernels]
  [heat]
    type = HeatConduction
    variable = T
    block = 'pellet inner_clad outer_clad'
    use_displaced_mesh = false
  []
[]

[Contact]
  [pellet_to_inner_clad]
    primary = inner_clad_left
    secondary = pellet_right
    formulation = penalty
    model = coulomb
    friction_coefficient = 0.001
    penalty = 1e14
    normalize_penalty = true
  []
  [inner_to_outer_clad]
    primary = outer_clad_left
    secondary = inner_clad_right
    formulation = penalty
    model = coulomb
    friction_coefficient = 0.001
    penalty = 1e14
    normalize_penalty = true
  []
[]

[BCs]
  [pellet_axis]
    type = DirichletBC
    variable = disp_x
    boundary = pellet_left
    value = 0
  []
  [inner_clad_inner]
    type = DirichletBC
    variable = disp_x
    boundary = inner_clad_left
    value = 0
  []
  [outer_clad_inner]
    type = DirichletBC
    variable = disp_x
    boundary = outer_clad_left
    value = 0
  []
  [pellet_closure]
    type = DirichletBC
    variable = disp_x
    boundary = pellet_right
    value = 1.1e-5
  []
  [inner_clad_closure]
    type = DirichletBC
    variable = disp_x
    boundary = inner_clad_right
    value = 1.1e-5
  []
  [pellet_bottom]
    type = DirichletBC
    variable = disp_y
    boundary = pellet_bottom
    value = 0
  []
  [inner_clad_bottom]
    type = DirichletBC
    variable = disp_y
    boundary = inner_clad_bottom
    value = 0
  []
  [outer_clad_bottom]
    type = DirichletBC
    variable = disp_y
    boundary = outer_clad_bottom
    value = 0
  []
  [pellet_top]
    type = DirichletBC
    variable = disp_y
    boundary = pellet_top
    value = 5e-6
  []
  [inner_clad_top]
    type = DirichletBC
    variable = disp_y
    boundary = inner_clad_top
    value = 1e-5
  []
  [outer_clad_top]
    type = DirichletBC
    variable = disp_y
    boundary = outer_clad_top
    value = 1.5e-5
  []
  [pellet_temperature]
    type = DirichletBC
    variable = T
    boundary = pellet_left
    value = 300
  []
  [inner_clad_temperature]
    type = DirichletBC
    variable = T
    boundary = inner_clad_left
    value = 300
  []
  [outer_clad_temperature]
    type = DirichletBC
    variable = T
    boundary = outer_clad_left
    value = 300
  []
[]

[Materials]
  [conductivity]
    type = GenericConstantMaterial
    block = 'pellet inner_clad outer_clad'
    prop_names = thermal_conductivity
    prop_values = 10
  []
  [elasticity]
    type = ComputeIsotropicElasticityTensor
    block = 'pellet inner_clad outer_clad'
    youngs_modulus = 2e11
    poissons_ratio = 0.3
  []
  [pellet_expansion]
    type = ComputeThermalExpansionEigenstrain
    block = pellet
    temperature = T
    stress_free_temperature = 300
    thermal_expansion_coeff = 0
    eigenstrain_name = pellet_thermal_strain
  []
  [inner_clad_expansion]
    type = ComputeThermalExpansionEigenstrain
    block = inner_clad
    temperature = T
    stress_free_temperature = 300
    thermal_expansion_coeff = 0
    eigenstrain_name = inner_clad_thermal_strain
  []
  [outer_clad_expansion]
    type = ComputeThermalExpansionEigenstrain
    block = outer_clad
    temperature = T
    stress_free_temperature = 300
    thermal_expansion_coeff = 0
    eigenstrain_name = outer_clad_thermal_strain
  []
  [pellet_stress]
    type = ComputeMultipleInelasticStress
    block = pellet
    inelastic_models = ''
  []
  [inner_clad_stress]
    type = ComputeMultipleInelasticStress
    block = inner_clad
    inelastic_models = ''
  []
  [outer_clad_stress]
    type = ComputeMultipleInelasticStress
    block = outer_clad
    inelastic_models = ''
  []
[]

[AuxKernels]
  [pellet_tangential_force_y]
    type = PenetrationAux
    variable = pellet_tangential_force_y
    boundary = pellet_right
    paired_boundary = inner_clad_left
    quantity = tangential_force_y
    use_displaced_mesh = true
  []
  [inner_clad_tangential_force_y]
    type = PenetrationAux
    variable = inner_clad_tangential_force_y
    boundary = inner_clad_right
    paired_boundary = outer_clad_left
    quantity = tangential_force_y
    use_displaced_mesh = true
  []
[]

[Executioner]
  type = Transient
  end_time = 1
  dt = 1
  solve_type = NEWTON
  line_search = basic
  nl_max_its = 80
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-8
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[VectorPostprocessors]
  [all_nodes]
    type = NodalValueSampler
    variable = 'T disp_x disp_y'
    sort_by = id
    use_displaced_mesh = false
  []
  [pellet_contact]
    type = NodalValueSampler
    boundary = pellet_right
    variable = 'T disp_x disp_y contact_pressure nodal_area penetration pellet_tangential_force_y'
    sort_by = y
    use_displaced_mesh = false
  []
  [inner_clad_contact]
    type = NodalValueSampler
    boundary = inner_clad_right
    variable = 'T disp_x disp_y contact_pressure nodal_area penetration inner_clad_tangential_force_y'
    sort_by = y
    use_displaced_mesh = false
  []
[]

[Outputs]
  csv = true
[]
