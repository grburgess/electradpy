#include "radiator.hh"

namespace emission
{
  
  class Synchrotron_Impulsive_Continuous : public virtual Radiator
  {
  public:

    virtual double source_function(double index){return index;};
    
  protected:

    
    std::array <double, n_grid_points> source_function_eval, V2, V3, sf_v2_ratio, V2_3_ratio;
    

    // precomputed factors for all CC rounds
    double factor1, DT; 
    
  private:
    
    void chang_cooper();
    virtual void emission_internal_computations(){};

    
  };




  class SIC_PowerLaw : public virtual Synchrotron_Impulsive_Continuous
  {

  public:

    void set_power_law(double ne, double gamma_min, double gamma_max, double p){this_ne=ne; this_gamma_min=gamma_min; this_gamma_max=gamma_max; this_p=p;};
    
    double source_function(double x, double ne ,double gamma_min,double gamma_max, double p);
    
  private:

    void emission_internal_computations();
    
    double this_ne, this_gamma_min, this_gamma_max, this_p;

  };

}

    
