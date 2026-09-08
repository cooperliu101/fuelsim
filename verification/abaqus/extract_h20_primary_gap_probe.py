from __future__ import print_function
import sys
from odbAccess import openOdb
odb=openOdb(sys.argv[1],readOnly=True)
try:
 labels=sorted(node.label for node in odb.rootAssembly.instances.values()[0].nodeSets['SECONDARY_FACE'].nodes)
 with open(sys.argv[2],'w') as out:
  out.write('step,node,gap\n')
  for name in ['BASE','PRIMARY_PLUS','PRIMARY_MINUS']:
   frame=odb.steps[name].frames[-1]
   keys=[key for key in frame.fieldOutputs.keys() if key.strip().startswith('COPEN')]
   if len(keys)!=1:raise RuntimeError('Expected one contact gap field')
   values=dict((v.nodeLabel,v.dataDouble) for v in frame.fieldOutputs[keys[0]].values)
   for node in labels:out.write('%s,%d,%.17g\n'%(name,node,values[node]))
finally:odb.close()
