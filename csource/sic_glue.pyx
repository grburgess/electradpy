# distutils: language = c++
# distutils: sources = src/synchrotron_impulsive_continuous.cxx

cdef extern from "incl/synchrotron_impulsive_continuous.hh" namespace "emission":
    cdef cppclass SIC_PowerLaw:
        SIC_PowerLaw() except +
        SIC_PowerLaw(double) except +
        void set_power_law(double, double, double, double)
