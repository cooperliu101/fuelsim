!include m1_fuel_cladding_gap_rz.i

# Distort only interior nodes. All axis, outer, top, bottom, primary, and
# secondary boundaries retain the original cylindrical geometry. The resulting
# Quad4 blocks cannot be represented as a tensor product of unique r and z
# coordinates.
[Mesh]
  [distort_interior]
    type = ParsedNodeTransformGenerator
    input = rename_block
    constant_names = 'pi'
    constant_expressions = '3.141592653589793238462643383279502884'
    x_function = 'if(x < 0.004121, x + 1.5e-5*sin(pi*x/0.00412)*sin(2*pi*y/0.010), x + 1.0e-5*sin(pi*(x-0.004122)/0.00057)*sin(2*pi*y/0.010020))'
    y_function = 'if(x < 0.004121, y + 5.0e-5*sin(pi*x/0.00412)*sin(pi*y/0.010), y + 5.0e-5*sin(pi*(x-0.004122)/0.00057)*sin(pi*y/0.010020))'
  []
[]
