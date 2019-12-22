import numpy
import pylab
from scipy.special import erfc, log_ndtr
from math import log,exp,sqrt

print exp(log(2)+log_ndtr(-2))
print exp(log(2)+log_ndtr(-(0.6*sqrt(2*13))))


Dhs = numpy.arange(5,100)

df = 0.6

variances = 1.0/2/Dhs
sigmas = numpy.sqrt(variances)

zs = df/sigmas

exact_log_Ps = -1*((log(2)+log_ndtr(-zs)))/log(10)
approx_log_Ps = -1*(-0.5*numpy.square(zs)-0.5*numpy.log(2*3.1415/4*numpy.square(zs)))/log(10)

pylab.plot(Dhs,exact_log_Ps,'b-')
pylab.plot(Dhs,approx_log_Ps,'r:')
pylab.show()



