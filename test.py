#!/usr/bin/env python3
'''Test code of ba-probability'''
import unittest
import sys
import numpy
from numpy.linalg import norm
from subprocess import check_call
from scipy.integrate import simps

rad2deg = 180.0/numpy.pi
deg2rad = 1.0/rad2deg



class Test_ba(unittest.TestCase):

    ba_values = numpy.array([[2.5053, 2.47447, 163.684],
                             [2.4744, 2.47886, 129.934],
                             [2.4788, 2.69089, 134.065]])

    cmd = '{} ba-test-sample.txt --output-file badf-test.txt'\
          ' --num_theta 500 --num_length 300 --d_range 1.5 3.5'\
          ' --theta_range 0.0 180.0 --d_width_factor 3.0 --q_width_factor 5.0'

    cmd_s = '{} ba-test-sample.txt --output-file badf-test-entropy.txt'\
            ' --num_theta 500 --num_length 300 --d_range 1.5 3.5'\
            ' --theta_range 0.0 180.0 --d_width_factor 3.0 --q_width_factor 5.0 --entropy'
    ba_exe = './ba-probability'

    lmin = 1.5 # in Anstrom
    lmax = 3.5
    qmin = 0.0 # in degree
    qmax = 180.0
    nl   = 300
    nq   = 500
    dl   = (lmax-lmin)/(nl-1)
    dq   = (qmax-qmin)/(nq-1)
    wl   = 3.0 * dl
    wq   = 5.0 * dq

    @classmethod
    def setUpClass(cls):
        print('Setting up tests...')
        numpy.savetxt('ba-test-sample.txt', cls.ba_values)
        try:
            check_call(cls.cmd.format(cls.ba_exe).split())
            check_call(cls.cmd_s.format(cls.ba_exe).split())
        except:
            print("################")
            print("Call ba-probability failed, note that must run test with GPU")
            print("################")
            sys.exit(1)
        cls.badf = cls.read_badf_file('badf-test.txt')
        cls.badf_s = cls.read_badf_file('badf-test-entropy.txt')

        l = numpy.concatenate((cls.ba_values[:,0], cls.ba_values[:,1]))
        q = numpy.concatenate((cls.ba_values[:,2], cls.ba_values[:,2]))
        l_grid, q_grid = numpy.meshgrid(numpy.linspace(cls.lmin,cls.lmax,cls.nl),
                                        numpy.linspace(cls.qmin,cls.qmax,cls.nq))
        cls.p     = numpy.zeros((cls.nq, cls.nl))
        cls.dpdl  = numpy.zeros((cls.nq, cls.nl))
        cls.dpdq  = numpy.zeros((cls.nq, cls.nl))
        cls.dpdlq = numpy.zeros((cls.nq, cls.nl))
        for l_mean, q_mean in zip(l, q):
            cls.p     +=     cls.get_p(l_grid, q_grid, cls.wl, cls.wq, l_mean, q_mean)
            cls.dpdl  +=  cls.get_dpdl(l_grid, q_grid, cls.wl, cls.wq, l_mean, q_mean)
            cls.dpdq  +=  cls.get_dpdq(l_grid, q_grid, cls.wl, cls.wq, l_mean, q_mean)
            cls.dpdlq += cls.get_dpdlq(l_grid, q_grid, cls.wl, cls.wq, l_mean, q_mean)
        cls.p    /= len(l)
        cls.dpdl /= len(l)
        cls.dpdq /= len(l)
        cls.dpdlq/= len(l)

        s_conv, ds_conv = cls.sin_convolution(cls)
        s_conv  = numpy.tile(s_conv, (cls.nl, 1)).T
        ds_conv  = numpy.tile(ds_conv, (cls.nl, 1)).T
        assert(s_conv.shape == ds_conv.shape == cls.p.shape)
        # Calculate p with entropy correction
        cls.p_s     = cls.p    /s_conv
        cls.dpdl_s  = cls.dpdl /s_conv
        cls.dpdq_s  = cls.dpdq /s_conv - cls.p   *ds_conv/s_conv**2
        cls.dpdlq_s = cls.dpdlq/s_conv - cls.dpdl*ds_conv/s_conv**2
        print("Setup done.")


    def setUp(self):
        pass


    def test_p(self):
        '''test the l2 norm error of p is smaller than 1e-5'''
        self.assertTrue(norm(self.badf.p-self.p)/norm(self.p) < 1e-5)


    def test_dpdl(self):
        '''test the l2 norm error of dpdl is smaller than 1e-5'''
        self.assertTrue(norm(self.badf.dpdl-self.dpdl)/norm(self.dpdl) < 1e-5)


    def test_dpdq(self):
        '''test the l2 norm error of dpdq is smaller than 1e-5'''
        self.assertTrue(norm(self.badf.dpdq-self.dpdq)/norm(self.dpdq) < 1e-5)


    def test_dpdlq(self):
        '''test the l2 norm error of dpdlq is smaller than 1e-5'''
        self.assertTrue(norm(self.badf.dpdlq-self.dpdlq)/norm(self.dpdlq) < 1e-5)


    def test_p_s(self):
        '''test the l2 norm error of p with entropy factor is smaller than 1e-5'''
        self.assertTrue(norm(self.badf_s.p-self.p_s)/norm(self.p_s) < 1e-5)

    def test_dpdl_s(self):
        '''test the l2 norm error of dpdl with entropy factor is smaller than 1e-5'''
        self.assertTrue(norm(self.badf_s.dpdl-self.dpdl_s)/norm(self.dpdl_s) < 1e-5)


    def test_dpdq_s(self):
        '''test the l2 norm error of dpdq with entropy factor is smaller than 1e-5'''
        self.assertTrue(norm(self.badf_s.dpdq-self.dpdq_s)/norm(self.dpdq_s) < 1e-5)


    def test_dpdlq_s(self):
        '''test the l2 norm error of dpdlq with entropy factor is smaller than 1e-5'''
        self.assertTrue(norm(self.badf_s.dpdlq-self.dpdlq_s)/norm(self.dpdlq_s) < 1e-5)


    def sin_convolution(self):
        ''' Calculate the convolution of sin(q) with gausian kernel over 0 to pi'''
        q = numpy.linspace(self.qmin, self.qmax, self.nq)
        q_rad = numpy.deg2rad(q)
        q_dense_rad = numpy.deg2rad(numpy.linspace(self.qmin, self.qmax, self.nq*10))
        normalizer = 1.0/numpy.sqrt(2*numpy.pi*(numpy.deg2rad(self.wq))**2)
        kernel = lambda q: numpy.exp(-0.5*(q**2/(numpy.deg2rad(self.wq))**2))*normalizer
        s    = numpy.zeros(len(q))
        dsdq = numpy.zeros(len(q))
        # calculate the convolution
        for i, qq in enumerate(q_rad):
            dsdq[i] = simps(numpy.cos(q_dense_rad) * kernel(qq - q_dense_rad), x = q_dense_rad)
            s[i]    = simps(numpy.sin(q_dense_rad) * kernel(qq - q_dense_rad), x = q_dense_rad)
        # convert unit from 1/rad to 1/deg
        dsdq = dsdq/rad2deg
        return s, dsdq


    def read_badf_file(path, var = False):
        ''' Reads the probability file generated by ba-probability. '''
        class Data:  pass
        data = Data()
        fid = open(path, 'rb')
        fid.readline()   # tells BA file (skip).
        data.nl, data.nq = [int(s) for s in fid.readline().split()[-2:]]
        data.llo, data.lhi = [float(s) for s in fid.readline().split()[-2:]]
        data.qlo, data.qhi = [float(s) for s in fid.readline().split()[-2:]]
        data.p = numpy.genfromtxt(fid, max_rows=data.nq)
        data.dpdl = numpy.genfromtxt(fid, max_rows=data.nq)
        data.dpdq = numpy.genfromtxt(fid, max_rows=data.nq)
        data.dpdlq = numpy.genfromtxt(fid, max_rows=data.nq)
        if var:
            data.var = numpy.genfromtxt(fid, max_rows = data.nq)
        fid.close()
        if len(data.p) != data.nq or len(data.p[0]) != data.nl or \
           len(data.dpdl) != data.nq or len(data.dpdl[0]) != data.nl or \
           len(data.dpdq) != data.nq or len(data.dpdq[0]) != data.nl or \
           len(data.dpdlq) != data.nq or len(data.dpdlq[0]) != data.nl:
            print("Dimension of the data read does not agree with the header in badf file!\n{}" \
                  .format(path))
            sys.exit(1)
        if var:
            if len(data.var) != data.nq or len(data.var[0]) != data.nl:
                print("Dimension of the data read does not agree with the header in badf file!\n{}" \
                      .format(path))
        return data


    def get_p(l, q, wl, wq, lmean, qmean):
        '''return the pdf of a 2d gaussian distribution'''
        normalizer = 1.0/(2*numpy.pi*wl*wq)
        return numpy.exp(-0.5*((l-lmean)**2/wl**2 +
                               (q-qmean)**2/wq**2))*normalizer


    def get_dpdl(l, q, wl, wq, lmean, qmean):
        '''return the dpdl of a 2d gaussian distribution function'''
        normalizer = 1.0/(2*numpy.pi*wl*wq)
        deri = -(l - lmean)/wl**2
        return deri*numpy.exp(-0.5*((l-lmean)**2/wl**2 +
                                    (q-qmean)**2/wq**2))*normalizer


    def get_dpdq(l, q, wl, wq, lmean, qmean):
        '''return the dpdq of a 2d gaussian distribution function'''
        normalizer = 1.0/(2*numpy.pi*wl*wq)
        deri = -(q - qmean)/wq**2
        return deri*numpy.exp(-0.5*((l-lmean)**2/wl**2 +
                                    (q-qmean)**2/wq**2))*normalizer


    def get_dpdlq(l, q, wl, wq, lmean, qmean):
        '''return the dpdlq of a 2d gaussian distribution function'''
        normalizer = 1.0/(2*numpy.pi*wl*wq)
        deri = (l - lmean)*(q - qmean)/(wq**2*wl**2)
        return deri*numpy.exp(-0.5*((l-lmean)**2/wl**2 +
                                    (q-qmean)**2/wq**2))*normalizer


if __name__ == "__main__":
    unittest.main()
