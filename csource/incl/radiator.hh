#include <math.h>
#include <stdlib.h>
#include <vector>
#include <array>
#include <algorithm>

//#include "gsl/gsl_integration.h"
#include "gsl/gsl_sf.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_sf_result.h"

#include <boost/range/combine.hpp>


namespace emission {

  class Radiator {
  public:
    Radiator();
    Radiator(double maxg);
    ~Radiator();
    
    std::vector<double> emission(double * energy, int n_photon_energies, int steps);
    std::vector<double> electrons();
    void reset();
    virtual double source_function(){return 0.;};

  protected:

    // number of electron grid points
    static const int n_grid_points = 300;
    int n_photon_energies;
    
    double DT,cool,B;
    
    // electron grids
    std::array <double, n_grid_points> gamma, gamma2, fgamma, fgammatp1, G, coolg, coolg2;
    

   
  private:
    
    // the total emission
    std::vector<double> radiation;

    
    double **synchrotron_kernel;
    double MAX_GAMMA;
    static constexpr double Bcritical = 4.14E13; // Gauss
    bool emitted, kernel_exists;
 
    
    // intialize the grid and the distribution function
    void initialize_electrons();
    virtual void emission_internal_computations(){};
    void compute_synchrotron_kernel(double * energy);
    void synchrotron_spectrum(double * energy);

    // inherited evolution solver
    virtual void chang_cooper(){};

    
    gsl_sf_result result;
    
    
  
    //double cool, step;
  
  };

}
