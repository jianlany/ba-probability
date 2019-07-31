#!/usr/bin/env python3
import numpy
from matplotlib import pyplot
pyplot.rc('font', size=8)
pyplot.rc('mathtext', fontset='cm')


def read_probability_file(path):
    class Data:  pass
    data = Data()
    fid = open('p.txt', 'rb')
    fid.readline()   # tells BA file (skip).
    data.M, data.N = [int(s) for s in fid.readline().split()[-2:]]
    data.xlo, data.xhi = [float(s) for s in fid.readline().split()[-2:]]
    data.ylo, data.yhi = [float(s) for s in fid.readline().split()[-2:]]
    fid.readline()  # should be p 
    
    data.p = numpy.genfromtxt(fid, max_rows=data.N)
    fid.readline()  # should be dpdl
    data.dpdl = numpy.genfromtxt(fid, max_rows=data.N)
    fid.readline()  # should be dpdq
    data.dpdq = numpy.genfromtxt(fid, max_rows=data.N)
    fid.readline()  # should be dpdql
    data.dpdql = numpy.genfromtxt(fid, max_rows=data.N)
    return data
     
data = read_probability_file('p.txt')
Y = 0.1*numpy.linspace(data.xlo, data.xhi, data.M)
X = numpy.linspace(data.ylo, data.yhi, data.N)
pyplot.figure(figsize=(3,3))
ax = pyplot.subplot(111)
ax.contourf(X, Y, data.p.T, 500, cmap='OrRd', vmin=1e-5, vmax=0.1)
ax.set_xlim(90.0, 180.0)
ax.set_ylim(0.22, 0.27)
pyplot.xticks([90.0, 120.0, 150.0, 180.0], 
              [r'$\frac{\pi}{2}$', r'$\frac{2\pi}{3}$', 
               r'$\frac{5\pi}{6}$', r'$\pi$'])

ax.set_xlabel(r'$\theta$ (rad)')
ax.set_ylabel(r'$l$ (nm)')
pyplot.tight_layout()
pyplot.savefig('p.png', dpi=300)
pyplot.close()

# Test that dpdl is accurate.
dl = (data.xhi - data.xlo) / (data.M - 1)
dpdl_num = (data.p[:,2:] - data.p[:,:-2]) / (2.0*dl)
e1 = numpy.linalg.norm(dpdl_num - data.dpdl[:,1:-1]) 
e1 /= numpy.linalg.norm(data.dpdl[:,1:-1]) 

dl = (data.xhi - data.xlo) / (data.M - 1)
dpdl_num = (data.p[:,2:] - data.p[:,:-2]) / (2.0*dl)
e1 = numpy.linalg.norm(dpdl_num - data.dpdl[:,1:-1]) 
e1 /= numpy.linalg.norm(data.dpdl[:,1:-1]) 
print('dpdl error norm: {:.5f}'.format(e1))

dq = (data.yhi - data.ylo) / (data.N - 1)
dpdq_num = (data.p[2:,:] - data.p[:-2,:]) / (2.0*dq)
e2 = numpy.linalg.norm(dpdq_num - data.dpdq[1:-1,:]) 
e2 /= numpy.linalg.norm(data.dpdq[1:-1,:]) 
print('dpdq error norm: {:.5f}'.format(e2))


dpdql_num = (dpdl_num[2:,:] - dpdl_num[:-2,:]) / (2.0*dq)
e3 = numpy.linalg.norm(dpdql_num - data.dpdql[1:-1,1:-1]) 
e3 /= numpy.linalg.norm(data.dpdql[1:-1,1:-1]) 
print('dpdql error norm: {:.5f}'.format(e3))

'''
pyplot.figure(figsize=(4,6))
vv = numpy.linspace(-1, 1, 20)
ax = pyplot.subplot(321)
ax.contour(dpdl_num.T, vv, vmin=-1, vmax=1)
ax = pyplot.subplot(322)
ax.contour(data.dpdl[:,1:-1].T, vv, vmin=-1, vmax=1)

ax = pyplot.subplot(323)
ax.contour(dpdq_num.T, 0.01*vv, vmin=-0.01, vmax=0.01)
ax = pyplot.subplot(324)
ax.contour(data.dpdq[1:-1,:].T, 0.01*vv, vmin=-0.01, vmax=0.01)

ax = pyplot.subplot(323)
ax.contour(dpdq_num.T, 0.01*vv, vmin=-0.01, vmax=0.01)
ax = pyplot.subplot(324)
ax.contour(data.dpdq[1:-1,:].T, 0.01*vv, vmin=-0.01, vmax=0.01)

ax = pyplot.subplot(325)
ax.contour(dpdql_num.T, 0.05*vv, vmin=-0.05, vmax=0.05)
ax = pyplot.subplot(326)
ax.contour(data.dpdql[1:-1,1:-1:].T, 0.05*vv, vmin=-0.05, vmax=0.05)

pyplot.tight_layout()
pyplot.savefig('derivatives.png', dpi=300)
'''
