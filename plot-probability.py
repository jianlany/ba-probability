#!/usr/bin/env python3
import numpy
from matplotlib import pyplot

pyplot.rc('font', size=8)
pyplot.rc('mathtext', fontset='cm')

p = numpy.genfromtxt('p.txt')
X = numpy.linspace(70.0, 180.0, p.shape[1])
Y = numpy.linspace(0.20, 0.29, p.shape[0])

pyplot.figure(figsize=(3,3))

ax = pyplot.subplot(111)
ax.contourf(X, Y, p, 500, cmap='OrRd', vmin=1e-5, vmax=0.014)
ax.set_xlim(90.0, 180.0)
ax.set_ylim(0.22, 0.27)
pyplot.xticks([90.0, 120.0, 150.0, 180.0], 
              [r'$\frac{\pi}{2}$', r'$\frac{2\pi}{3}$', 
               r'$\frac{5\pi}{6}$', r'$\pi$'])

ax.set_xlabel(r'$\theta$ (rad)')
ax.set_ylabel(r'$l$ (nm)')
pyplot.tight_layout()
pyplot.savefig('p.png', dpi=300)
