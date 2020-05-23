#!/usr/bin/env python3
import numpy
from matplotlib import pyplot
from ba_probability import read_probability_file, verify_derivatives
pyplot.rc('font', size=8)
pyplot.rc('mathtext', fontset='cm')


     
data = read_probability_file('p.txt')
Y = 0.1*numpy.linspace(data.xlo, data.xhi, data.M)
X = numpy.linspace(data.ylo, data.yhi, data.N)
pyplot.figure(figsize=(4,4))
ax = pyplot.subplot(111)
vv = numpy.linspace(-10, 1, 74)
pd = numpy.log10(data.p.T)
jmax, imax = numpy.unravel_index(pd.argmax(axis=None), pd.shape)
print(pd.max(), 'at', X[imax], Y[jmax])
ax.contour(X, Y, pd, vv, cmap='Purples', linewidths=0.4)
ax.set_xlim(0, 180.0)
#ax.set_ylim(0.22, 0.27)
pyplot.xticks([0, 30, 60, 90.0, 120.0, 150.0, 180.0], 
              ['0', r'$\frac{\pi}{6}$', r'$\frac{\pi}{3}$',
                r'$\frac{\pi}{2}$', r'$\frac{2\pi}{3}$', 
               r'$\frac{5\pi}{6}$', r'$\pi$'])

ax.set_xlabel(r'$\theta$ (rad)')
ax.set_ylabel(r'$l$ (nm)')
pyplot.tight_layout()
pyplot.savefig('p.png', dpi=300)
pyplot.close()

dx = (data.xhi - data.xlo) / (data.M - 1)
dy = (data.yhi - data.ylo) / (data.N - 1)
verify_derivatives(data.p, data.dpdl, data.dpdq, data.dpdql, dx, dy, 'p-derivatives.png')
