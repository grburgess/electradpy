#ifndef RADIATOR_H
#define RADIATOR_H

#include <math.h>
#include <stdlib.h>
#include <vector>
#include <array>
#include <algorithm>

//#include "gsl/gsl_integration.h"
#include "gsl/gsl_sf.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_sf_result.h"
#include <iostream>


namespace emission {

  class Radiator {
  public:
    Radiator() {std::cout << "Empty constructor" << std::endl;};
    Radiator(double maxg);
    ~Radiator();
    
    std::vector<double> photons(double * energy, int n_photon_energies_now, int steps, double this_B);
    std::vector<double> electrons();
    std::vector<double> gamma_grid();
    void reset();
    void initialize_electrons();
    virtual double source_function(){return 0.;};

  protected:

    // number of electron grid points
    static const int n_grid_points = 300;
    
    int n_photon_energies;
    
    double DT,cool,B;
    
    // electron grids
    std::array <double, n_grid_points> gamma, gamma2, fgamma, fgammatp1, G, coolg, coolg2;
    std::array <double, n_grid_points-1> delta_grid;
    
   
  private:
    
    // the total emission
    std::vector<double> *radiation;
    //std::array <double, n_photons_max> radiation;

    
    double **synchrotron_kernel;
    double MAX_GAMMA;
    static constexpr double Bcritical = 4.14E13; // Gauss
    bool emitted, kernel_exists;
 
    
    // intialize the grid and the distribution function
    
    virtual void emission_internal_computations(){};
    void compute_synchrotron_kernel(double * energy);
    void synchrotron_spectrum(double * energy);

    // inherited evolution solver
    virtual void chang_cooper(){};

    
    gsl_sf_result result;
    
    
  
    //double cool, step;
  
  };

}


#endif
