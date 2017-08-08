# distutils: language = c++
# distutils: sources = csource/src/synchrotron_impulsive_continuous.cxx

from libcpp.vector cimport vector
cimport numpy as np

cdef extern from "incl/synchrotron_impulsive_continuous.hh" namespace "emission":
    cdef cppclass SIC_PowerLaw:
        #SIC_PowerLaw() except +
        SIC_PowerLaw(double) except +
        void set_power_law(double, double, double, double)
        void initialize_electrons()
        void reset()
        vector[double] photons(double *, int, int, double)
        vector[double] electrons()
        vector[double] gamma_grid()


cdef class pySIC_Powerlaw:
    cdef SIC_PowerLaw* c_sic_pl # hold a c++ ref
    def __cinit__(self, double max_gamma):

        self.c_sic_pl = new SIC_PowerLaw(max_gamma)
        self.c_sic_pl.initialize_electrons()
    def __dealloc__(self):

        del self.c_sic_pl

        
        
    def set_power_law(self, ne, gamma_min, gamma_max, p):

        self.c_sic_pl.set_power_law(ne,gamma_min,gamma_max,p)

    def photons(self, np.ndarray[double, ndim=1, mode="c"] ene not None, double B, int steps):

        return self.c_sic_pl.photons(<double*> np.PyArray_DATA(ene),ene.shape[0], steps, B)

    @property
    def electrons(self):
        
        return self.c_sic_pl.electrons()
    @property
    def gamma(self):

        return self.c_sic_pl.gamma_grid()
    
    def reset(self):

        self.c_sic_pl.reset()
