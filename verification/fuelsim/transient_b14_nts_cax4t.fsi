[Case]
  version = 3
  problem = transient
  geometry = axisymmetric_rz
[]
[Mesh]
  type = exodus
  file = ../meshes/b14_nts_cax4t.e
[]
[Materials]
  [solid]
    [thermal]
      function = constant_thermophysical
      conductivity = 1e-6
      density = 1
      specific_heat = 1
    []
    [elasticity]
      function = constant_isotropic
      young_modulus = 1e6
      poisson_ratio = 0.3
    []
  []
[]
[Regions]
  [fuel]
    block = fuel
    element = cax4t
    material = solid
    strain = finite
    initial_temperature = 600
    volumetric_heat_source = 0
  []
  [cladding]
    block = clad
    element = cax4t
    material = solid
    strain = finite
    initial_temperature = 600
    volumetric_heat_source = 0
  []
[]
[Contact]
  [fuel_cladding]
    primary = clad_left
    secondary = fuel_right
    [thermal]
      law = affine
      conductance = 1000
      discretization = node_to_surface
    []
    [mechanical]
      formulation = penalty
      penalty = 1e12
      mu = 0.2
      slip_tolerance = 0.005
      discretization = node_to_surface
      sliding = finite
    []
  []
[]
[TimeFunctions]
  [sliding]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0 1e-5 1e-5 -1e-5
  []
  [radial_1]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0 0 0 0
  []
  [temperature_1]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 600 600 600 600
  []
  [radial_2]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.000167 0.000167 0.000167 0.000167
  []
  [temperature_2]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 600.6866666666666 600.6866666666666 600.6866666666666 600.6866666666666
  []
  [radial_3]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.000167 0.000167 0.000195125 0.000195125
  []
  [temperature_3]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 625.6866666666666 625.6866666666666 625.6866666666666 625.6866666666666
  []
  [radial_4]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0 0 0.000028125000000000003 0.000028125000000000003
  []
  [temperature_4]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 625 625 625 625
  []
  [radial_5]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.000334 0.000334 0.000334 0.000334
  []
  [temperature_5]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 601.3733333333333 601.3733333333333 601.3733333333333 601.3733333333333
  []
  [radial_6]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.000334 0.000334 0.000362125 0.000362125
  []
  [temperature_6]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 626.3733333333333 626.3733333333333 626.3733333333333 626.3733333333333
  []
  [radial_7]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.000501 0.000501 0.000501 0.000501
  []
  [temperature_7]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 602.06 602.06 602.06 602.06
  []
  [radial_8]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.000501 0.000501 0.000529125 0.000529125
  []
  [temperature_8]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 627.06 627.06 627.06 627.06
  []
  [radial_9]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.000668 0.000668 0.000668 0.000668
  []
  [temperature_9]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 602.7466666666667 602.7466666666667 602.7466666666667 602.7466666666667
  []
  [radial_10]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.000668 0.000668 0.0006961249999999999 0.0006961249999999999
  []
  [temperature_10]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 627.7466666666667 627.7466666666667 627.7466666666667 627.7466666666667
  []
  [radial_11]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.000835 0.000835 0.000835 0.000835
  []
  [temperature_11]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 603.4333333333333 603.4333333333333 603.4333333333333 603.4333333333333
  []
  [radial_12]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.000835 0.000835 0.000863125 0.000863125
  []
  [temperature_12]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 628.4333333333333 628.4333333333333 628.4333333333333 628.4333333333333
  []
  [radial_13]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001002 0.001002 0.001002 0.001002
  []
  [temperature_13]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 604.12 604.12 604.12 604.12
  []
  [radial_14]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001002 0.001002 0.0010301250000000002 0.0010301250000000002
  []
  [temperature_14]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 629.12 629.12 629.12 629.12
  []
  [radial_15]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.000167 0.000167 0.0002295 0.0002295
  []
  [temperature_15]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 650.6866666666666 650.6866666666666 650.6866666666666 650.6866666666666
  []
  [radial_16]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0 0 0.0000625 0.0000625
  []
  [temperature_16]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 650 650 650 650
  []
  [radial_17]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.000334 0.000334 0.0003965 0.0003965
  []
  [temperature_17]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 651.3733333333333 651.3733333333333 651.3733333333333 651.3733333333333
  []
  [radial_18]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.000501 0.000501 0.0005635000000000001 0.0005635000000000001
  []
  [temperature_18]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 652.06 652.06 652.06 652.06
  []
  [radial_19]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.000668 0.000668 0.0007305 0.0007305
  []
  [temperature_19]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 652.7466666666667 652.7466666666667 652.7466666666667 652.7466666666667
  []
  [radial_20]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.000835 0.000835 0.0008975000000000001 0.0008975000000000001
  []
  [temperature_20]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 653.4333333333333 653.4333333333333 653.4333333333333 653.4333333333333
  []
  [radial_21]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001002 0.001002 0.0010645 0.0010645
  []
  [temperature_21]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 654.12 654.12 654.12 654.12
  []
  [radial_22]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.000167 0.000167 0.00027012499999999997 0.00027012499999999997
  []
  [temperature_22]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 675.6866666666666 675.6866666666666 675.6866666666666 675.6866666666666
  []
  [radial_23]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0 0 0.000103125 0.000103125
  []
  [temperature_23]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 675 675 675 675
  []
  [radial_24]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.000334 0.000334 0.00043712499999999996 0.00043712499999999996
  []
  [temperature_24]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 676.3733333333333 676.3733333333333 676.3733333333333 676.3733333333333
  []
  [radial_25]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.000501 0.000501 0.0006041250000000001 0.0006041250000000001
  []
  [temperature_25]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 677.06 677.06 677.06 677.06
  []
  [radial_26]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.000668 0.000668 0.000771125 0.000771125
  []
  [temperature_26]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 677.7466666666667 677.7466666666667 677.7466666666667 677.7466666666667
  []
  [radial_27]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.000835 0.000835 0.000938125 0.000938125
  []
  [temperature_27]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 678.4333333333333 678.4333333333333 678.4333333333333 678.4333333333333
  []
  [radial_28]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001002 0.001002 0.001105125 0.001105125
  []
  [temperature_28]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 679.12 679.12 679.12 679.12
  []
  [radial_29]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.000167 0.000167 0.000317 0.000317
  []
  [temperature_29]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 700.6866666666666 700.6866666666666 700.6866666666666 700.6866666666666
  []
  [radial_30]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0 0 0.00015000000000000001 0.00015000000000000001
  []
  [temperature_30]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 700 700 700 700
  []
  [radial_31]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.000334 0.000334 0.000484 0.000484
  []
  [temperature_31]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 701.3733333333333 701.3733333333333 701.3733333333333 701.3733333333333
  []
  [radial_32]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.000501 0.000501 0.0006510000000000001 0.0006510000000000001
  []
  [temperature_32]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 702.06 702.06 702.06 702.06
  []
  [radial_33]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.000668 0.000668 0.000818 0.000818
  []
  [temperature_33]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 702.7466666666667 702.7466666666667 702.7466666666667 702.7466666666667
  []
  [radial_34]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.000835 0.000835 0.000985 0.000985
  []
  [temperature_34]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 703.4333333333333 703.4333333333333 703.4333333333333 703.4333333333333
  []
  [radial_35]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001002 0.001002 0.001152 0.001152
  []
  [temperature_35]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 704.12 704.12 704.12 704.12
  []
  [radial_36]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001 0.001 0.001 0.001
  []
  [temperature_36]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 499.5879 499.5879 499.5879 499.5879
  []
  [radial_37]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001 0.001 0.001 0.001
  []
  [temperature_37]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 499.55935 499.55935 499.55935 499.55935
  []
  [radial_38]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001 0.001 0.001022048008 0.001022048008
  []
  [temperature_38]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 505.57135 505.57135 505.57135 505.57135
  []
  [radial_39]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001 0.001 0.001022048008 0.001022048008
  []
  [temperature_39]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 505.5999 505.5999 505.5999 505.5999
  []
  [radial_40]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001 0.001 0.001 0.001
  []
  [temperature_40]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 499.5308 499.5308 499.5308 499.5308
  []
  [radial_41]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001 0.001 0.001022048008 0.001022048008
  []
  [temperature_41]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 505.5428 505.5428 505.5428 505.5428
  []
  [radial_42]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001 0.001 0.001048112032 0.001048112032
  []
  [temperature_42]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 511.58335 511.58335 511.58335 511.58335
  []
  [radial_43]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001 0.001 0.001048112032 0.001048112032
  []
  [temperature_43]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 511.6119 511.6119 511.6119 511.6119
  []
  [radial_44]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001 0.001 0.001048112032 0.001048112032
  []
  [temperature_44]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 511.5548 511.5548 511.5548 511.5548
  []
  [radial_45]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001 0.001 0.001078192072 0.001078192072
  []
  [temperature_45]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 517.5953499999999 517.5953499999999 517.5953499999999 517.5953499999999
  []
  [radial_46]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001 0.001 0.001078192072 0.001078192072
  []
  [temperature_46]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 517.6238999999999 517.6238999999999 517.6238999999999 517.6238999999999
  []
  [radial_47]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001 0.001 0.001078192072 0.001078192072
  []
  [temperature_47]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 517.5668 517.5668 517.5668 517.5668
  []
  [radial_48]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001 0.001 0.0011122881279999998 0.0011122881279999998
  []
  [temperature_48]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 523.60735 523.60735 523.60735 523.60735
  []
  [radial_49]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001 0.001 0.0011122881279999998 0.0011122881279999998
  []
  [temperature_49]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 523.6359 523.6359 523.6359 523.6359
  []
  [radial_50]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001 0.001 0.0011122881279999998 0.0011122881279999998
  []
  [temperature_50]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 523.5788 523.5788 523.5788 523.5788
  []
  [radial_51]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001 0.001 0.0011504002 0.0011504002
  []
  [temperature_51]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 529.6193499999999 529.6193499999999 529.6193499999999 529.6193499999999
  []
  [radial_52]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001 0.001 0.0011504002 0.0011504002
  []
  [temperature_52]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 529.6478999999999 529.6478999999999 529.6478999999999 529.6478999999999
  []
  [radial_53]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001 0.001 0.0011504002 0.0011504002
  []
  [temperature_53]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 600 529.5908 529.5908 529.5908 529.5908
  []
[]
[BoundaryConditions]
  [radial_1]
    type = dirichlet
    boundary = n_1
    field = radial_displacement
    value = 1
    function = radial_1
  []
  [axial_1]
    type = dirichlet
    boundary = n_1
    field = axial_displacement
    value = 1
    function = sliding
  []
  [temperature_1]
    type = dirichlet
    boundary = n_1
    field = temperature
    value = 1
    function = temperature_1
  []
  [radial_2]
    type = dirichlet
    boundary = n_2
    field = radial_displacement
    value = 1
    function = radial_2
  []
  [axial_2]
    type = dirichlet
    boundary = n_2
    field = axial_displacement
    value = 1
    function = sliding
  []
  [temperature_2]
    type = dirichlet
    boundary = n_2
    field = temperature
    value = 1
    function = temperature_2
  []
  [radial_3]
    type = dirichlet
    boundary = n_3
    field = radial_displacement
    value = 1
    function = radial_3
  []
  [axial_3]
    type = dirichlet
    boundary = n_3
    field = axial_displacement
    value = 1
    function = sliding
  []
  [temperature_3]
    type = dirichlet
    boundary = n_3
    field = temperature
    value = 1
    function = temperature_3
  []
  [radial_4]
    type = dirichlet
    boundary = n_4
    field = radial_displacement
    value = 1
    function = radial_4
  []
  [axial_4]
    type = dirichlet
    boundary = n_4
    field = axial_displacement
    value = 1
    function = sliding
  []
  [temperature_4]
    type = dirichlet
    boundary = n_4
    field = temperature
    value = 1
    function = temperature_4
  []
  [radial_5]
    type = dirichlet
    boundary = n_5
    field = radial_displacement
    value = 1
    function = radial_5
  []
  [axial_5]
    type = dirichlet
    boundary = n_5
    field = axial_displacement
    value = 1
    function = sliding
  []
  [temperature_5]
    type = dirichlet
    boundary = n_5
    field = temperature
    value = 1
    function = temperature_5
  []
  [radial_6]
    type = dirichlet
    boundary = n_6
    field = radial_displacement
    value = 1
    function = radial_6
  []
  [axial_6]
    type = dirichlet
    boundary = n_6
    field = axial_displacement
    value = 1
    function = sliding
  []
  [temperature_6]
    type = dirichlet
    boundary = n_6
    field = temperature
    value = 1
    function = temperature_6
  []
  [radial_7]
    type = dirichlet
    boundary = n_7
    field = radial_displacement
    value = 1
    function = radial_7
  []
  [axial_7]
    type = dirichlet
    boundary = n_7
    field = axial_displacement
    value = 1
    function = sliding
  []
  [temperature_7]
    type = dirichlet
    boundary = n_7
    field = temperature
    value = 1
    function = temperature_7
  []
  [radial_8]
    type = dirichlet
    boundary = n_8
    field = radial_displacement
    value = 1
    function = radial_8
  []
  [axial_8]
    type = dirichlet
    boundary = n_8
    field = axial_displacement
    value = 1
    function = sliding
  []
  [temperature_8]
    type = dirichlet
    boundary = n_8
    field = temperature
    value = 1
    function = temperature_8
  []
  [radial_9]
    type = dirichlet
    boundary = n_9
    field = radial_displacement
    value = 1
    function = radial_9
  []
  [axial_9]
    type = dirichlet
    boundary = n_9
    field = axial_displacement
    value = 1
    function = sliding
  []
  [temperature_9]
    type = dirichlet
    boundary = n_9
    field = temperature
    value = 1
    function = temperature_9
  []
  [radial_10]
    type = dirichlet
    boundary = n_10
    field = radial_displacement
    value = 1
    function = radial_10
  []
  [axial_10]
    type = dirichlet
    boundary = n_10
    field = axial_displacement
    value = 1
    function = sliding
  []
  [temperature_10]
    type = dirichlet
    boundary = n_10
    field = temperature
    value = 1
    function = temperature_10
  []
  [radial_11]
    type = dirichlet
    boundary = n_11
    field = radial_displacement
    value = 1
    function = radial_11
  []
  [axial_11]
    type = dirichlet
    boundary = n_11
    field = axial_displacement
    value = 1
    function = sliding
  []
  [temperature_11]
    type = dirichlet
    boundary = n_11
    field = temperature
    value = 1
    function = temperature_11
  []
  [radial_12]
    type = dirichlet
    boundary = n_12
    field = radial_displacement
    value = 1
    function = radial_12
  []
  [axial_12]
    type = dirichlet
    boundary = n_12
    field = axial_displacement
    value = 1
    function = sliding
  []
  [temperature_12]
    type = dirichlet
    boundary = n_12
    field = temperature
    value = 1
    function = temperature_12
  []
  [radial_13]
    type = dirichlet
    boundary = n_13
    field = radial_displacement
    value = 1
    function = radial_13
  []
  [axial_13]
    type = dirichlet
    boundary = n_13
    field = axial_displacement
    value = 1
    function = sliding
  []
  [temperature_13]
    type = dirichlet
    boundary = n_13
    field = temperature
    value = 1
    function = temperature_13
  []
  [radial_14]
    type = dirichlet
    boundary = n_14
    field = radial_displacement
    value = 1
    function = radial_14
  []
  [axial_14]
    type = dirichlet
    boundary = n_14
    field = axial_displacement
    value = 1
    function = sliding
  []
  [temperature_14]
    type = dirichlet
    boundary = n_14
    field = temperature
    value = 1
    function = temperature_14
  []
  [radial_15]
    type = dirichlet
    boundary = n_15
    field = radial_displacement
    value = 1
    function = radial_15
  []
  [axial_15]
    type = dirichlet
    boundary = n_15
    field = axial_displacement
    value = 1
    function = sliding
  []
  [temperature_15]
    type = dirichlet
    boundary = n_15
    field = temperature
    value = 1
    function = temperature_15
  []
  [radial_16]
    type = dirichlet
    boundary = n_16
    field = radial_displacement
    value = 1
    function = radial_16
  []
  [axial_16]
    type = dirichlet
    boundary = n_16
    field = axial_displacement
    value = 1
    function = sliding
  []
  [temperature_16]
    type = dirichlet
    boundary = n_16
    field = temperature
    value = 1
    function = temperature_16
  []
  [radial_17]
    type = dirichlet
    boundary = n_17
    field = radial_displacement
    value = 1
    function = radial_17
  []
  [axial_17]
    type = dirichlet
    boundary = n_17
    field = axial_displacement
    value = 1
    function = sliding
  []
  [temperature_17]
    type = dirichlet
    boundary = n_17
    field = temperature
    value = 1
    function = temperature_17
  []
  [radial_18]
    type = dirichlet
    boundary = n_18
    field = radial_displacement
    value = 1
    function = radial_18
  []
  [axial_18]
    type = dirichlet
    boundary = n_18
    field = axial_displacement
    value = 1
    function = sliding
  []
  [temperature_18]
    type = dirichlet
    boundary = n_18
    field = temperature
    value = 1
    function = temperature_18
  []
  [radial_19]
    type = dirichlet
    boundary = n_19
    field = radial_displacement
    value = 1
    function = radial_19
  []
  [axial_19]
    type = dirichlet
    boundary = n_19
    field = axial_displacement
    value = 1
    function = sliding
  []
  [temperature_19]
    type = dirichlet
    boundary = n_19
    field = temperature
    value = 1
    function = temperature_19
  []
  [radial_20]
    type = dirichlet
    boundary = n_20
    field = radial_displacement
    value = 1
    function = radial_20
  []
  [axial_20]
    type = dirichlet
    boundary = n_20
    field = axial_displacement
    value = 1
    function = sliding
  []
  [temperature_20]
    type = dirichlet
    boundary = n_20
    field = temperature
    value = 1
    function = temperature_20
  []
  [radial_21]
    type = dirichlet
    boundary = n_21
    field = radial_displacement
    value = 1
    function = radial_21
  []
  [axial_21]
    type = dirichlet
    boundary = n_21
    field = axial_displacement
    value = 1
    function = sliding
  []
  [temperature_21]
    type = dirichlet
    boundary = n_21
    field = temperature
    value = 1
    function = temperature_21
  []
  [radial_22]
    type = dirichlet
    boundary = n_22
    field = radial_displacement
    value = 1
    function = radial_22
  []
  [axial_22]
    type = dirichlet
    boundary = n_22
    field = axial_displacement
    value = 1
    function = sliding
  []
  [temperature_22]
    type = dirichlet
    boundary = n_22
    field = temperature
    value = 1
    function = temperature_22
  []
  [radial_23]
    type = dirichlet
    boundary = n_23
    field = radial_displacement
    value = 1
    function = radial_23
  []
  [axial_23]
    type = dirichlet
    boundary = n_23
    field = axial_displacement
    value = 1
    function = sliding
  []
  [temperature_23]
    type = dirichlet
    boundary = n_23
    field = temperature
    value = 1
    function = temperature_23
  []
  [radial_24]
    type = dirichlet
    boundary = n_24
    field = radial_displacement
    value = 1
    function = radial_24
  []
  [axial_24]
    type = dirichlet
    boundary = n_24
    field = axial_displacement
    value = 1
    function = sliding
  []
  [temperature_24]
    type = dirichlet
    boundary = n_24
    field = temperature
    value = 1
    function = temperature_24
  []
  [radial_25]
    type = dirichlet
    boundary = n_25
    field = radial_displacement
    value = 1
    function = radial_25
  []
  [axial_25]
    type = dirichlet
    boundary = n_25
    field = axial_displacement
    value = 1
    function = sliding
  []
  [temperature_25]
    type = dirichlet
    boundary = n_25
    field = temperature
    value = 1
    function = temperature_25
  []
  [radial_26]
    type = dirichlet
    boundary = n_26
    field = radial_displacement
    value = 1
    function = radial_26
  []
  [axial_26]
    type = dirichlet
    boundary = n_26
    field = axial_displacement
    value = 1
    function = sliding
  []
  [temperature_26]
    type = dirichlet
    boundary = n_26
    field = temperature
    value = 1
    function = temperature_26
  []
  [radial_27]
    type = dirichlet
    boundary = n_27
    field = radial_displacement
    value = 1
    function = radial_27
  []
  [axial_27]
    type = dirichlet
    boundary = n_27
    field = axial_displacement
    value = 1
    function = sliding
  []
  [temperature_27]
    type = dirichlet
    boundary = n_27
    field = temperature
    value = 1
    function = temperature_27
  []
  [radial_28]
    type = dirichlet
    boundary = n_28
    field = radial_displacement
    value = 1
    function = radial_28
  []
  [axial_28]
    type = dirichlet
    boundary = n_28
    field = axial_displacement
    value = 1
    function = sliding
  []
  [temperature_28]
    type = dirichlet
    boundary = n_28
    field = temperature
    value = 1
    function = temperature_28
  []
  [radial_29]
    type = dirichlet
    boundary = n_29
    field = radial_displacement
    value = 1
    function = radial_29
  []
  [axial_29]
    type = dirichlet
    boundary = n_29
    field = axial_displacement
    value = 1
    function = sliding
  []
  [temperature_29]
    type = dirichlet
    boundary = n_29
    field = temperature
    value = 1
    function = temperature_29
  []
  [radial_30]
    type = dirichlet
    boundary = n_30
    field = radial_displacement
    value = 1
    function = radial_30
  []
  [axial_30]
    type = dirichlet
    boundary = n_30
    field = axial_displacement
    value = 1
    function = sliding
  []
  [temperature_30]
    type = dirichlet
    boundary = n_30
    field = temperature
    value = 1
    function = temperature_30
  []
  [radial_31]
    type = dirichlet
    boundary = n_31
    field = radial_displacement
    value = 1
    function = radial_31
  []
  [axial_31]
    type = dirichlet
    boundary = n_31
    field = axial_displacement
    value = 1
    function = sliding
  []
  [temperature_31]
    type = dirichlet
    boundary = n_31
    field = temperature
    value = 1
    function = temperature_31
  []
  [radial_32]
    type = dirichlet
    boundary = n_32
    field = radial_displacement
    value = 1
    function = radial_32
  []
  [axial_32]
    type = dirichlet
    boundary = n_32
    field = axial_displacement
    value = 1
    function = sliding
  []
  [temperature_32]
    type = dirichlet
    boundary = n_32
    field = temperature
    value = 1
    function = temperature_32
  []
  [radial_33]
    type = dirichlet
    boundary = n_33
    field = radial_displacement
    value = 1
    function = radial_33
  []
  [axial_33]
    type = dirichlet
    boundary = n_33
    field = axial_displacement
    value = 1
    function = sliding
  []
  [temperature_33]
    type = dirichlet
    boundary = n_33
    field = temperature
    value = 1
    function = temperature_33
  []
  [radial_34]
    type = dirichlet
    boundary = n_34
    field = radial_displacement
    value = 1
    function = radial_34
  []
  [axial_34]
    type = dirichlet
    boundary = n_34
    field = axial_displacement
    value = 1
    function = sliding
  []
  [temperature_34]
    type = dirichlet
    boundary = n_34
    field = temperature
    value = 1
    function = temperature_34
  []
  [radial_35]
    type = dirichlet
    boundary = n_35
    field = radial_displacement
    value = 1
    function = radial_35
  []
  [axial_35]
    type = dirichlet
    boundary = n_35
    field = axial_displacement
    value = 1
    function = sliding
  []
  [temperature_35]
    type = dirichlet
    boundary = n_35
    field = temperature
    value = 1
    function = temperature_35
  []
  [radial_36]
    type = dirichlet
    boundary = n_36
    field = radial_displacement
    value = 1
    function = radial_36
  []
  [axial_36]
    type = dirichlet
    boundary = n_36
    field = axial_displacement
    value = 0
  []
  [temperature_36]
    type = dirichlet
    boundary = n_36
    field = temperature
    value = 1
    function = temperature_36
  []
  [radial_37]
    type = dirichlet
    boundary = n_37
    field = radial_displacement
    value = 1
    function = radial_37
  []
  [axial_37]
    type = dirichlet
    boundary = n_37
    field = axial_displacement
    value = 0
  []
  [temperature_37]
    type = dirichlet
    boundary = n_37
    field = temperature
    value = 1
    function = temperature_37
  []
  [radial_38]
    type = dirichlet
    boundary = n_38
    field = radial_displacement
    value = 1
    function = radial_38
  []
  [axial_38]
    type = dirichlet
    boundary = n_38
    field = axial_displacement
    value = 0
  []
  [temperature_38]
    type = dirichlet
    boundary = n_38
    field = temperature
    value = 1
    function = temperature_38
  []
  [radial_39]
    type = dirichlet
    boundary = n_39
    field = radial_displacement
    value = 1
    function = radial_39
  []
  [axial_39]
    type = dirichlet
    boundary = n_39
    field = axial_displacement
    value = 0
  []
  [temperature_39]
    type = dirichlet
    boundary = n_39
    field = temperature
    value = 1
    function = temperature_39
  []
  [radial_40]
    type = dirichlet
    boundary = n_40
    field = radial_displacement
    value = 1
    function = radial_40
  []
  [axial_40]
    type = dirichlet
    boundary = n_40
    field = axial_displacement
    value = 0
  []
  [temperature_40]
    type = dirichlet
    boundary = n_40
    field = temperature
    value = 1
    function = temperature_40
  []
  [radial_41]
    type = dirichlet
    boundary = n_41
    field = radial_displacement
    value = 1
    function = radial_41
  []
  [axial_41]
    type = dirichlet
    boundary = n_41
    field = axial_displacement
    value = 0
  []
  [temperature_41]
    type = dirichlet
    boundary = n_41
    field = temperature
    value = 1
    function = temperature_41
  []
  [radial_42]
    type = dirichlet
    boundary = n_42
    field = radial_displacement
    value = 1
    function = radial_42
  []
  [axial_42]
    type = dirichlet
    boundary = n_42
    field = axial_displacement
    value = 0
  []
  [temperature_42]
    type = dirichlet
    boundary = n_42
    field = temperature
    value = 1
    function = temperature_42
  []
  [radial_43]
    type = dirichlet
    boundary = n_43
    field = radial_displacement
    value = 1
    function = radial_43
  []
  [axial_43]
    type = dirichlet
    boundary = n_43
    field = axial_displacement
    value = 0
  []
  [temperature_43]
    type = dirichlet
    boundary = n_43
    field = temperature
    value = 1
    function = temperature_43
  []
  [radial_44]
    type = dirichlet
    boundary = n_44
    field = radial_displacement
    value = 1
    function = radial_44
  []
  [axial_44]
    type = dirichlet
    boundary = n_44
    field = axial_displacement
    value = 0
  []
  [temperature_44]
    type = dirichlet
    boundary = n_44
    field = temperature
    value = 1
    function = temperature_44
  []
  [radial_45]
    type = dirichlet
    boundary = n_45
    field = radial_displacement
    value = 1
    function = radial_45
  []
  [axial_45]
    type = dirichlet
    boundary = n_45
    field = axial_displacement
    value = 0
  []
  [temperature_45]
    type = dirichlet
    boundary = n_45
    field = temperature
    value = 1
    function = temperature_45
  []
  [radial_46]
    type = dirichlet
    boundary = n_46
    field = radial_displacement
    value = 1
    function = radial_46
  []
  [axial_46]
    type = dirichlet
    boundary = n_46
    field = axial_displacement
    value = 0
  []
  [temperature_46]
    type = dirichlet
    boundary = n_46
    field = temperature
    value = 1
    function = temperature_46
  []
  [radial_47]
    type = dirichlet
    boundary = n_47
    field = radial_displacement
    value = 1
    function = radial_47
  []
  [axial_47]
    type = dirichlet
    boundary = n_47
    field = axial_displacement
    value = 0
  []
  [temperature_47]
    type = dirichlet
    boundary = n_47
    field = temperature
    value = 1
    function = temperature_47
  []
  [radial_48]
    type = dirichlet
    boundary = n_48
    field = radial_displacement
    value = 1
    function = radial_48
  []
  [axial_48]
    type = dirichlet
    boundary = n_48
    field = axial_displacement
    value = 0
  []
  [temperature_48]
    type = dirichlet
    boundary = n_48
    field = temperature
    value = 1
    function = temperature_48
  []
  [radial_49]
    type = dirichlet
    boundary = n_49
    field = radial_displacement
    value = 1
    function = radial_49
  []
  [axial_49]
    type = dirichlet
    boundary = n_49
    field = axial_displacement
    value = 0
  []
  [temperature_49]
    type = dirichlet
    boundary = n_49
    field = temperature
    value = 1
    function = temperature_49
  []
  [radial_50]
    type = dirichlet
    boundary = n_50
    field = radial_displacement
    value = 1
    function = radial_50
  []
  [axial_50]
    type = dirichlet
    boundary = n_50
    field = axial_displacement
    value = 0
  []
  [temperature_50]
    type = dirichlet
    boundary = n_50
    field = temperature
    value = 1
    function = temperature_50
  []
  [radial_51]
    type = dirichlet
    boundary = n_51
    field = radial_displacement
    value = 1
    function = radial_51
  []
  [axial_51]
    type = dirichlet
    boundary = n_51
    field = axial_displacement
    value = 0
  []
  [temperature_51]
    type = dirichlet
    boundary = n_51
    field = temperature
    value = 1
    function = temperature_51
  []
  [radial_52]
    type = dirichlet
    boundary = n_52
    field = radial_displacement
    value = 1
    function = radial_52
  []
  [axial_52]
    type = dirichlet
    boundary = n_52
    field = axial_displacement
    value = 0
  []
  [temperature_52]
    type = dirichlet
    boundary = n_52
    field = temperature
    value = 1
    function = temperature_52
  []
  [radial_53]
    type = dirichlet
    boundary = n_53
    field = radial_displacement
    value = 1
    function = radial_53
  []
  [axial_53]
    type = dirichlet
    boundary = n_53
    field = axial_displacement
    value = 0
  []
  [temperature_53]
    type = dirichlet
    boundary = n_53
    field = temperature
    value = 1
    function = temperature_53
  []
[]
[Executioner]
  type = transient
  end_time = 4
  initial_time_step = 1
  minimum_time_step = 1
  maximum_time_step = 1
  growth_factor = 1
  cutback_factor = 0.5
  maximum_cutbacks = 0
  load_ramp_time = 0
  include_thermal_time_term = false
[]
[Solver]
  linear_solver = direct
  direct_factorization = mumps
  line_search = backtracking
  maximum_iterations = 100
  absolute_tolerance = 1e-10
  relative_tolerance = 1e-11
  step_tolerance = 1e-16
[]
[Outputs]
  console = true
  csv = transient_b14_nts_cax4t_summary.csv
  exodus = transient_b14_nts_cax4t_results.e
[]
