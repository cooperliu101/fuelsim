[Case]
  version = 3
  problem = transient
  geometry = axisymmetric_rz
[]
[Mesh]
  type = exodus
  file = ../meshes/b14_nts_cax8t.e
[]
[Materials]
  [solid]
    [thermal]
      function = constant_thermophysical
      conductivity = 1e-12
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
    element = cax8t
    material = solid
    strain = finite
    initial_temperature = 600
    volumetric_heat_source = 0
  []
  [cladding]
    block = clad
    element = cax8t
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
  [radial_54]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.0000835 0.0000835 0.0000835 0.0000835
  []
  [radial_55]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.000167 0.000167 0.00018028125 0.00018028125
  []
  [radial_56]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.0000835 0.0000835 0.00011162499999999999 0.00011162499999999999
  []
  [radial_57]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0 0 0.000013281250000000001 0.000013281250000000001
  []
  [radial_58]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.0002505 0.0002505 0.0002505 0.0002505
  []
  [radial_59]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.000334 0.000334 0.00034728125 0.00034728125
  []
  [radial_60]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.0002505 0.0002505 0.000278625 0.000278625
  []
  [radial_61]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.00041749999999999996 0.00041749999999999996 0.00041749999999999996 0.00041749999999999996
  []
  [radial_62]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.000501 0.000501 0.0005142812500000001 0.0005142812500000001
  []
  [radial_63]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.00041749999999999996 0.00041749999999999996 0.00044562499999999995 0.00044562499999999995
  []
  [radial_64]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.0005845000000000001 0.0005845000000000001 0.0005845000000000001 0.0005845000000000001
  []
  [radial_65]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.000668 0.000668 0.00068128125 0.00068128125
  []
  [radial_66]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.0005845000000000001 0.0005845000000000001 0.000612625 0.000612625
  []
  [radial_67]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.0007515 0.0007515 0.0007515 0.0007515
  []
  [radial_68]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.000835 0.000835 0.0008482812500000001 0.0008482812500000001
  []
  [radial_69]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.0007515 0.0007515 0.0007796249999999999 0.0007796249999999999
  []
  [radial_70]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.0009185 0.0009185 0.0009185 0.0009185
  []
  [radial_71]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001002 0.001002 0.00101528125 0.00101528125
  []
  [radial_72]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.0009185 0.0009185 0.000946625 0.000946625
  []
  [radial_73]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.000167 0.000167 0.00021153125 0.00021153125
  []
  [radial_74]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.0000835 0.0000835 0.000146 0.000146
  []
  [radial_75]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0 0 0.000044531249999999995 0.000044531249999999995
  []
  [radial_76]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.000334 0.000334 0.00037853124999999996 0.00037853124999999996
  []
  [radial_77]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.0002505 0.0002505 0.000313 0.000313
  []
  [radial_78]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.000501 0.000501 0.00054553125 0.00054553125
  []
  [radial_79]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.00041749999999999996 0.00041749999999999996 0.00047999999999999996 0.00047999999999999996
  []
  [radial_80]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.000668 0.000668 0.00071253125 0.00071253125
  []
  [radial_81]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.0005845000000000001 0.0005845000000000001 0.0006470000000000001 0.0006470000000000001
  []
  [radial_82]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.000835 0.000835 0.00087953125 0.00087953125
  []
  [radial_83]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.0007515 0.0007515 0.000814 0.000814
  []
  [radial_84]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001002 0.001002 0.0010465312500000002 0.0010465312500000002
  []
  [radial_85]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.0009185 0.0009185 0.000981 0.000981
  []
  [radial_86]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.000167 0.000167 0.00024903125 0.00024903125
  []
  [radial_87]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.0000835 0.0000835 0.000186625 0.000186625
  []
  [radial_88]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0 0 0.00008203125 0.00008203125
  []
  [radial_89]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.000334 0.000334 0.00041603125 0.00041603125
  []
  [radial_90]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.0002505 0.0002505 0.000353625 0.000353625
  []
  [radial_91]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.000501 0.000501 0.00058303125 0.00058303125
  []
  [radial_92]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.00041749999999999996 0.00041749999999999996 0.000520625 0.000520625
  []
  [radial_93]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.000668 0.000668 0.00075003125 0.00075003125
  []
  [radial_94]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.0005845000000000001 0.0005845000000000001 0.0006876250000000001 0.0006876250000000001
  []
  [radial_95]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.000835 0.000835 0.0009170312500000001 0.0009170312500000001
  []
  [radial_96]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.0007515 0.0007515 0.000854625 0.000854625
  []
  [radial_97]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001002 0.001002 0.0010840312500000001 0.0010840312500000001
  []
  [radial_98]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.0009185 0.0009185 0.001021625 0.001021625
  []
  [radial_99]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.000167 0.000167 0.00029278125 0.00029278125
  []
  [radial_100]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.0000835 0.0000835 0.0002335 0.0002335
  []
  [radial_101]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0 0 0.00012578125000000002 0.00012578125000000002
  []
  [radial_102]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.000334 0.000334 0.00045978125 0.00045978125
  []
  [radial_103]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.0002505 0.0002505 0.00040050000000000003 0.00040050000000000003
  []
  [radial_104]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.000501 0.000501 0.0006267812500000001 0.0006267812500000001
  []
  [radial_105]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.00041749999999999996 0.00041749999999999996 0.0005675 0.0005675
  []
  [radial_106]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.000668 0.000668 0.00079378125 0.00079378125
  []
  [radial_107]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.0005845000000000001 0.0005845000000000001 0.0007345000000000001 0.0007345000000000001
  []
  [radial_108]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.000835 0.000835 0.00096078125 0.00096078125
  []
  [radial_109]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.0007515 0.0007515 0.0009015000000000001 0.0009015000000000001
  []
  [radial_110]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001002 0.001002 0.00112778125 0.00112778125
  []
  [radial_111]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.0009185 0.0009185 0.0010685 0.0010685
  []
  [radial_112]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001 0.001 0.001 0.001
  []
  [radial_113]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001 0.001 0.001010522002 0.001010522002
  []
  [radial_114]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001 0.001 0.001022048008 0.001022048008
  []
  [radial_115]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001 0.001 0.001010522002 0.001010522002
  []
  [radial_116]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001 0.001 0.001 0.001
  []
  [radial_117]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001 0.001 0.001010522002 0.001010522002
  []
  [radial_118]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001 0.001 0.001022048008 0.001022048008
  []
  [radial_119]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001 0.001 0.001034578018 0.001034578018
  []
  [radial_120]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001 0.001 0.001048112032 0.001048112032
  []
  [radial_121]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001 0.001 0.001034578018 0.001034578018
  []
  [radial_122]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001 0.001 0.001034578018 0.001034578018
  []
  [radial_123]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001 0.001 0.001048112032 0.001048112032
  []
  [radial_124]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001 0.001 0.00106265005 0.00106265005
  []
  [radial_125]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001 0.001 0.001078192072 0.001078192072
  []
  [radial_126]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001 0.001 0.00106265005 0.00106265005
  []
  [radial_127]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001 0.001 0.00106265005 0.00106265005
  []
  [radial_128]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001 0.001 0.001078192072 0.001078192072
  []
  [radial_129]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001 0.001 0.001094738098 0.001094738098
  []
  [radial_130]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001 0.001 0.0011122881279999998 0.0011122881279999998
  []
  [radial_131]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001 0.001 0.001094738098 0.001094738098
  []
  [radial_132]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001 0.001 0.001094738098 0.001094738098
  []
  [radial_133]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001 0.001 0.0011122881279999998 0.0011122881279999998
  []
  [radial_134]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001 0.001 0.001130842162 0.001130842162
  []
  [radial_135]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001 0.001 0.0011504002 0.0011504002
  []
  [radial_136]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001 0.001 0.001130842162 0.001130842162
  []
  [radial_137]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001 0.001 0.001130842162 0.001130842162
  []
  [radial_138]
    type = piecewise_linear
    times = 0 1 2 3 4
    values = 0 0.001 0.001 0.0011504002 0.0011504002
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
  [radial_54]
    type = dirichlet
    boundary = n_54
    field = radial_displacement
    value = 1
    function = radial_54
  []
  [axial_54]
    type = dirichlet
    boundary = n_54
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_55]
    type = dirichlet
    boundary = n_55
    field = radial_displacement
    value = 1
    function = radial_55
  []
  [axial_55]
    type = dirichlet
    boundary = n_55
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_56]
    type = dirichlet
    boundary = n_56
    field = radial_displacement
    value = 1
    function = radial_56
  []
  [axial_56]
    type = dirichlet
    boundary = n_56
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_57]
    type = dirichlet
    boundary = n_57
    field = radial_displacement
    value = 1
    function = radial_57
  []
  [axial_57]
    type = dirichlet
    boundary = n_57
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_58]
    type = dirichlet
    boundary = n_58
    field = radial_displacement
    value = 1
    function = radial_58
  []
  [axial_58]
    type = dirichlet
    boundary = n_58
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_59]
    type = dirichlet
    boundary = n_59
    field = radial_displacement
    value = 1
    function = radial_59
  []
  [axial_59]
    type = dirichlet
    boundary = n_59
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_60]
    type = dirichlet
    boundary = n_60
    field = radial_displacement
    value = 1
    function = radial_60
  []
  [axial_60]
    type = dirichlet
    boundary = n_60
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_61]
    type = dirichlet
    boundary = n_61
    field = radial_displacement
    value = 1
    function = radial_61
  []
  [axial_61]
    type = dirichlet
    boundary = n_61
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_62]
    type = dirichlet
    boundary = n_62
    field = radial_displacement
    value = 1
    function = radial_62
  []
  [axial_62]
    type = dirichlet
    boundary = n_62
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_63]
    type = dirichlet
    boundary = n_63
    field = radial_displacement
    value = 1
    function = radial_63
  []
  [axial_63]
    type = dirichlet
    boundary = n_63
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_64]
    type = dirichlet
    boundary = n_64
    field = radial_displacement
    value = 1
    function = radial_64
  []
  [axial_64]
    type = dirichlet
    boundary = n_64
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_65]
    type = dirichlet
    boundary = n_65
    field = radial_displacement
    value = 1
    function = radial_65
  []
  [axial_65]
    type = dirichlet
    boundary = n_65
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_66]
    type = dirichlet
    boundary = n_66
    field = radial_displacement
    value = 1
    function = radial_66
  []
  [axial_66]
    type = dirichlet
    boundary = n_66
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_67]
    type = dirichlet
    boundary = n_67
    field = radial_displacement
    value = 1
    function = radial_67
  []
  [axial_67]
    type = dirichlet
    boundary = n_67
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_68]
    type = dirichlet
    boundary = n_68
    field = radial_displacement
    value = 1
    function = radial_68
  []
  [axial_68]
    type = dirichlet
    boundary = n_68
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_69]
    type = dirichlet
    boundary = n_69
    field = radial_displacement
    value = 1
    function = radial_69
  []
  [axial_69]
    type = dirichlet
    boundary = n_69
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_70]
    type = dirichlet
    boundary = n_70
    field = radial_displacement
    value = 1
    function = radial_70
  []
  [axial_70]
    type = dirichlet
    boundary = n_70
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_71]
    type = dirichlet
    boundary = n_71
    field = radial_displacement
    value = 1
    function = radial_71
  []
  [axial_71]
    type = dirichlet
    boundary = n_71
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_72]
    type = dirichlet
    boundary = n_72
    field = radial_displacement
    value = 1
    function = radial_72
  []
  [axial_72]
    type = dirichlet
    boundary = n_72
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_73]
    type = dirichlet
    boundary = n_73
    field = radial_displacement
    value = 1
    function = radial_73
  []
  [axial_73]
    type = dirichlet
    boundary = n_73
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_74]
    type = dirichlet
    boundary = n_74
    field = radial_displacement
    value = 1
    function = radial_74
  []
  [axial_74]
    type = dirichlet
    boundary = n_74
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_75]
    type = dirichlet
    boundary = n_75
    field = radial_displacement
    value = 1
    function = radial_75
  []
  [axial_75]
    type = dirichlet
    boundary = n_75
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_76]
    type = dirichlet
    boundary = n_76
    field = radial_displacement
    value = 1
    function = radial_76
  []
  [axial_76]
    type = dirichlet
    boundary = n_76
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_77]
    type = dirichlet
    boundary = n_77
    field = radial_displacement
    value = 1
    function = radial_77
  []
  [axial_77]
    type = dirichlet
    boundary = n_77
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_78]
    type = dirichlet
    boundary = n_78
    field = radial_displacement
    value = 1
    function = radial_78
  []
  [axial_78]
    type = dirichlet
    boundary = n_78
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_79]
    type = dirichlet
    boundary = n_79
    field = radial_displacement
    value = 1
    function = radial_79
  []
  [axial_79]
    type = dirichlet
    boundary = n_79
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_80]
    type = dirichlet
    boundary = n_80
    field = radial_displacement
    value = 1
    function = radial_80
  []
  [axial_80]
    type = dirichlet
    boundary = n_80
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_81]
    type = dirichlet
    boundary = n_81
    field = radial_displacement
    value = 1
    function = radial_81
  []
  [axial_81]
    type = dirichlet
    boundary = n_81
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_82]
    type = dirichlet
    boundary = n_82
    field = radial_displacement
    value = 1
    function = radial_82
  []
  [axial_82]
    type = dirichlet
    boundary = n_82
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_83]
    type = dirichlet
    boundary = n_83
    field = radial_displacement
    value = 1
    function = radial_83
  []
  [axial_83]
    type = dirichlet
    boundary = n_83
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_84]
    type = dirichlet
    boundary = n_84
    field = radial_displacement
    value = 1
    function = radial_84
  []
  [axial_84]
    type = dirichlet
    boundary = n_84
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_85]
    type = dirichlet
    boundary = n_85
    field = radial_displacement
    value = 1
    function = radial_85
  []
  [axial_85]
    type = dirichlet
    boundary = n_85
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_86]
    type = dirichlet
    boundary = n_86
    field = radial_displacement
    value = 1
    function = radial_86
  []
  [axial_86]
    type = dirichlet
    boundary = n_86
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_87]
    type = dirichlet
    boundary = n_87
    field = radial_displacement
    value = 1
    function = radial_87
  []
  [axial_87]
    type = dirichlet
    boundary = n_87
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_88]
    type = dirichlet
    boundary = n_88
    field = radial_displacement
    value = 1
    function = radial_88
  []
  [axial_88]
    type = dirichlet
    boundary = n_88
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_89]
    type = dirichlet
    boundary = n_89
    field = radial_displacement
    value = 1
    function = radial_89
  []
  [axial_89]
    type = dirichlet
    boundary = n_89
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_90]
    type = dirichlet
    boundary = n_90
    field = radial_displacement
    value = 1
    function = radial_90
  []
  [axial_90]
    type = dirichlet
    boundary = n_90
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_91]
    type = dirichlet
    boundary = n_91
    field = radial_displacement
    value = 1
    function = radial_91
  []
  [axial_91]
    type = dirichlet
    boundary = n_91
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_92]
    type = dirichlet
    boundary = n_92
    field = radial_displacement
    value = 1
    function = radial_92
  []
  [axial_92]
    type = dirichlet
    boundary = n_92
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_93]
    type = dirichlet
    boundary = n_93
    field = radial_displacement
    value = 1
    function = radial_93
  []
  [axial_93]
    type = dirichlet
    boundary = n_93
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_94]
    type = dirichlet
    boundary = n_94
    field = radial_displacement
    value = 1
    function = radial_94
  []
  [axial_94]
    type = dirichlet
    boundary = n_94
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_95]
    type = dirichlet
    boundary = n_95
    field = radial_displacement
    value = 1
    function = radial_95
  []
  [axial_95]
    type = dirichlet
    boundary = n_95
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_96]
    type = dirichlet
    boundary = n_96
    field = radial_displacement
    value = 1
    function = radial_96
  []
  [axial_96]
    type = dirichlet
    boundary = n_96
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_97]
    type = dirichlet
    boundary = n_97
    field = radial_displacement
    value = 1
    function = radial_97
  []
  [axial_97]
    type = dirichlet
    boundary = n_97
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_98]
    type = dirichlet
    boundary = n_98
    field = radial_displacement
    value = 1
    function = radial_98
  []
  [axial_98]
    type = dirichlet
    boundary = n_98
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_99]
    type = dirichlet
    boundary = n_99
    field = radial_displacement
    value = 1
    function = radial_99
  []
  [axial_99]
    type = dirichlet
    boundary = n_99
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_100]
    type = dirichlet
    boundary = n_100
    field = radial_displacement
    value = 1
    function = radial_100
  []
  [axial_100]
    type = dirichlet
    boundary = n_100
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_101]
    type = dirichlet
    boundary = n_101
    field = radial_displacement
    value = 1
    function = radial_101
  []
  [axial_101]
    type = dirichlet
    boundary = n_101
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_102]
    type = dirichlet
    boundary = n_102
    field = radial_displacement
    value = 1
    function = radial_102
  []
  [axial_102]
    type = dirichlet
    boundary = n_102
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_103]
    type = dirichlet
    boundary = n_103
    field = radial_displacement
    value = 1
    function = radial_103
  []
  [axial_103]
    type = dirichlet
    boundary = n_103
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_104]
    type = dirichlet
    boundary = n_104
    field = radial_displacement
    value = 1
    function = radial_104
  []
  [axial_104]
    type = dirichlet
    boundary = n_104
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_105]
    type = dirichlet
    boundary = n_105
    field = radial_displacement
    value = 1
    function = radial_105
  []
  [axial_105]
    type = dirichlet
    boundary = n_105
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_106]
    type = dirichlet
    boundary = n_106
    field = radial_displacement
    value = 1
    function = radial_106
  []
  [axial_106]
    type = dirichlet
    boundary = n_106
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_107]
    type = dirichlet
    boundary = n_107
    field = radial_displacement
    value = 1
    function = radial_107
  []
  [axial_107]
    type = dirichlet
    boundary = n_107
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_108]
    type = dirichlet
    boundary = n_108
    field = radial_displacement
    value = 1
    function = radial_108
  []
  [axial_108]
    type = dirichlet
    boundary = n_108
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_109]
    type = dirichlet
    boundary = n_109
    field = radial_displacement
    value = 1
    function = radial_109
  []
  [axial_109]
    type = dirichlet
    boundary = n_109
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_110]
    type = dirichlet
    boundary = n_110
    field = radial_displacement
    value = 1
    function = radial_110
  []
  [axial_110]
    type = dirichlet
    boundary = n_110
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_111]
    type = dirichlet
    boundary = n_111
    field = radial_displacement
    value = 1
    function = radial_111
  []
  [axial_111]
    type = dirichlet
    boundary = n_111
    field = axial_displacement
    value = 1
    function = sliding
  []
  [radial_112]
    type = dirichlet
    boundary = n_112
    field = radial_displacement
    value = 1
    function = radial_112
  []
  [axial_112]
    type = dirichlet
    boundary = n_112
    field = axial_displacement
    value = 0
  []
  [radial_113]
    type = dirichlet
    boundary = n_113
    field = radial_displacement
    value = 1
    function = radial_113
  []
  [axial_113]
    type = dirichlet
    boundary = n_113
    field = axial_displacement
    value = 0
  []
  [radial_114]
    type = dirichlet
    boundary = n_114
    field = radial_displacement
    value = 1
    function = radial_114
  []
  [axial_114]
    type = dirichlet
    boundary = n_114
    field = axial_displacement
    value = 0
  []
  [radial_115]
    type = dirichlet
    boundary = n_115
    field = radial_displacement
    value = 1
    function = radial_115
  []
  [axial_115]
    type = dirichlet
    boundary = n_115
    field = axial_displacement
    value = 0
  []
  [radial_116]
    type = dirichlet
    boundary = n_116
    field = radial_displacement
    value = 1
    function = radial_116
  []
  [axial_116]
    type = dirichlet
    boundary = n_116
    field = axial_displacement
    value = 0
  []
  [radial_117]
    type = dirichlet
    boundary = n_117
    field = radial_displacement
    value = 1
    function = radial_117
  []
  [axial_117]
    type = dirichlet
    boundary = n_117
    field = axial_displacement
    value = 0
  []
  [radial_118]
    type = dirichlet
    boundary = n_118
    field = radial_displacement
    value = 1
    function = radial_118
  []
  [axial_118]
    type = dirichlet
    boundary = n_118
    field = axial_displacement
    value = 0
  []
  [radial_119]
    type = dirichlet
    boundary = n_119
    field = radial_displacement
    value = 1
    function = radial_119
  []
  [axial_119]
    type = dirichlet
    boundary = n_119
    field = axial_displacement
    value = 0
  []
  [radial_120]
    type = dirichlet
    boundary = n_120
    field = radial_displacement
    value = 1
    function = radial_120
  []
  [axial_120]
    type = dirichlet
    boundary = n_120
    field = axial_displacement
    value = 0
  []
  [radial_121]
    type = dirichlet
    boundary = n_121
    field = radial_displacement
    value = 1
    function = radial_121
  []
  [axial_121]
    type = dirichlet
    boundary = n_121
    field = axial_displacement
    value = 0
  []
  [radial_122]
    type = dirichlet
    boundary = n_122
    field = radial_displacement
    value = 1
    function = radial_122
  []
  [axial_122]
    type = dirichlet
    boundary = n_122
    field = axial_displacement
    value = 0
  []
  [radial_123]
    type = dirichlet
    boundary = n_123
    field = radial_displacement
    value = 1
    function = radial_123
  []
  [axial_123]
    type = dirichlet
    boundary = n_123
    field = axial_displacement
    value = 0
  []
  [radial_124]
    type = dirichlet
    boundary = n_124
    field = radial_displacement
    value = 1
    function = radial_124
  []
  [axial_124]
    type = dirichlet
    boundary = n_124
    field = axial_displacement
    value = 0
  []
  [radial_125]
    type = dirichlet
    boundary = n_125
    field = radial_displacement
    value = 1
    function = radial_125
  []
  [axial_125]
    type = dirichlet
    boundary = n_125
    field = axial_displacement
    value = 0
  []
  [radial_126]
    type = dirichlet
    boundary = n_126
    field = radial_displacement
    value = 1
    function = radial_126
  []
  [axial_126]
    type = dirichlet
    boundary = n_126
    field = axial_displacement
    value = 0
  []
  [radial_127]
    type = dirichlet
    boundary = n_127
    field = radial_displacement
    value = 1
    function = radial_127
  []
  [axial_127]
    type = dirichlet
    boundary = n_127
    field = axial_displacement
    value = 0
  []
  [radial_128]
    type = dirichlet
    boundary = n_128
    field = radial_displacement
    value = 1
    function = radial_128
  []
  [axial_128]
    type = dirichlet
    boundary = n_128
    field = axial_displacement
    value = 0
  []
  [radial_129]
    type = dirichlet
    boundary = n_129
    field = radial_displacement
    value = 1
    function = radial_129
  []
  [axial_129]
    type = dirichlet
    boundary = n_129
    field = axial_displacement
    value = 0
  []
  [radial_130]
    type = dirichlet
    boundary = n_130
    field = radial_displacement
    value = 1
    function = radial_130
  []
  [axial_130]
    type = dirichlet
    boundary = n_130
    field = axial_displacement
    value = 0
  []
  [radial_131]
    type = dirichlet
    boundary = n_131
    field = radial_displacement
    value = 1
    function = radial_131
  []
  [axial_131]
    type = dirichlet
    boundary = n_131
    field = axial_displacement
    value = 0
  []
  [radial_132]
    type = dirichlet
    boundary = n_132
    field = radial_displacement
    value = 1
    function = radial_132
  []
  [axial_132]
    type = dirichlet
    boundary = n_132
    field = axial_displacement
    value = 0
  []
  [radial_133]
    type = dirichlet
    boundary = n_133
    field = radial_displacement
    value = 1
    function = radial_133
  []
  [axial_133]
    type = dirichlet
    boundary = n_133
    field = axial_displacement
    value = 0
  []
  [radial_134]
    type = dirichlet
    boundary = n_134
    field = radial_displacement
    value = 1
    function = radial_134
  []
  [axial_134]
    type = dirichlet
    boundary = n_134
    field = axial_displacement
    value = 0
  []
  [radial_135]
    type = dirichlet
    boundary = n_135
    field = radial_displacement
    value = 1
    function = radial_135
  []
  [axial_135]
    type = dirichlet
    boundary = n_135
    field = axial_displacement
    value = 0
  []
  [radial_136]
    type = dirichlet
    boundary = n_136
    field = radial_displacement
    value = 1
    function = radial_136
  []
  [axial_136]
    type = dirichlet
    boundary = n_136
    field = axial_displacement
    value = 0
  []
  [radial_137]
    type = dirichlet
    boundary = n_137
    field = radial_displacement
    value = 1
    function = radial_137
  []
  [axial_137]
    type = dirichlet
    boundary = n_137
    field = axial_displacement
    value = 0
  []
  [radial_138]
    type = dirichlet
    boundary = n_138
    field = radial_displacement
    value = 1
    function = radial_138
  []
  [axial_138]
    type = dirichlet
    boundary = n_138
    field = axial_displacement
    value = 0
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
  csv = transient_b14_nts_cax8t_summary.csv
  exodus = transient_b14_nts_cax8t_results.e
[]
