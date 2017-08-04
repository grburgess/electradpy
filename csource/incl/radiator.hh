#include <math.h>
#include <stdlib.h>
#include <vector>


namespace emission {

  class Radiator {
  public:
    Radiator();
    Radiator(int n_grid_points);
    ~Radiator();
    std::vector<double> emission(double * energy);
    double source_function();
    
  private:
    int n_grid_points;
    void initialize_electrons();
    double * chang_cooper();


    double *gamma = new double[n_grid_points];
    double *gamma2= new double[n_grid_points];
    double *fgamma = new double[n_grid_points];
    double *fgammatp1 = new double[n_grid_points];
    double *G = new double [n_grid_points+1];

    double **synchrotron_kernel;
  
    //double cool, step;
  


}
