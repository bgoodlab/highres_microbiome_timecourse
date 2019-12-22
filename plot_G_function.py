import barcode_utils
import numpy
import pylab

Dbins, raw_G_function_map = barcode_utils.parse_raw_G_function_map()

G_function_map = barcode_utils.parse_G_function_map()

print raw_G_function_map.keys()
print G_function_map.keys()

Ds = numpy.logspace(1,4,30)

pylab.figure(figsize=(5,3))

for sample_name in G_function_map:

    Gs = [G_function_map[sample_name](D) for D in Ds]
    
    pylab.loglog(Ds,Gs,'-')
  
pylab.xlabel('D')
pylab.ylabel('G(D)')  
pylab.savefig('G_functions.pdf',bbox_inches='tight')