[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 1
  ny = 1
  xmax = 0.001
  ymax = 0.001
  coord_type = RZ
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
[]

[AuxVariables]
  [effective_creep_strain]
    order = CONSTANT
    family = MONOMIAL
  []
  [effective_plastic_strain]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Physics/SolidMechanics/QuasiStatic]
  [all]
    strain = SMALL
    incremental = true
    add_variables = true
    use_automatic_differentiation = true
    generate_output = 'stress_xx stress_yy stress_zz creep_strain_xx creep_strain_yy creep_strain_zz plastic_strain_xx plastic_strain_yy plastic_strain_zz'
  []
[]

[Functions]
  [axial_pull]
    type = ParsedFunction
    expression = '0.001*0.002*t'
  []
[]

[AuxKernels]
  [effective_creep_strain]
    type = ADMaterialRealAux
    variable = effective_creep_strain
    property = effective_creep_strain
  []
  [effective_plastic_strain]
    type = ADMaterialRealAux
    variable = effective_plastic_strain
    property = effective_plastic_strain
  []
[]

[BCs]
  [axis]
    type = ADDirichletBC
    variable = disp_x
    boundary = left
    value = 0
  []
  [bottom]
    type = ADDirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
  []
  [top_pull]
    type = ADFunctionDirichletBC
    variable = disp_y
    boundary = top
    function = axial_pull
  []
[]

[Materials]
  [elasticity]
    type = ADComputeIsotropicElasticityTensor
    youngs_modulus = 2e11
    poissons_ratio = 0.3
  []
  [stress]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'creep plasticity'
    perform_finite_strain_rotations = false
    max_iterations = 100
    relative_tolerance = 1e-12
    absolute_tolerance = 1e-6
  []
  [creep]
    type = ADPowerLawCreepStressUpdate
    coefficient = 1e-28
    n_exponent = 3
    m_exponent = 0
    activation_energy = 0
  []
  [plasticity]
    type = ADIsotropicPlasticityStressUpdate
    yield_stress = 2e8
    hardening_constant = 2e9
  []
[]

[Postprocessors]
  [axial_stress]
    type = ElementAverageValue
    variable = stress_yy
  []
  [radial_stress]
    type = ElementAverageValue
    variable = stress_xx
  []
  [hoop_stress]
    type = ElementAverageValue
    variable = stress_zz
  []
  [effective_creep]
    type = ElementAverageValue
    variable = effective_creep_strain
  []
  [effective_plastic]
    type = ElementAverageValue
    variable = effective_plastic_strain
  []
  [axial_creep]
    type = ElementAverageValue
    variable = creep_strain_yy
  []
  [radial_creep]
    type = ElementAverageValue
    variable = creep_strain_xx
  []
  [hoop_creep]
    type = ElementAverageValue
    variable = creep_strain_zz
  []
  [axial_plastic]
    type = ElementAverageValue
    variable = plastic_strain_yy
  []
  [radial_plastic]
    type = ElementAverageValue
    variable = plastic_strain_xx
  []
  [hoop_plastic]
    type = ElementAverageValue
    variable = plastic_strain_zz
  []
  [axial_displacement]
    type = SideAverageValue
    variable = disp_y
    boundary = top
    use_displaced_mesh = false
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
  dt = 0.1
  end_time = 1
  nl_abs_tol = 1e-10
  nl_rel_tol = 1e-10
  automatic_scaling = true
  compute_scaling_once = false
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
[]

[Outputs]
  csv = true
[]
