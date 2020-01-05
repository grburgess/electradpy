#ifndef SYNCHROTRON_STOCHASTIC_H
#define SYNCHROTRON_STOCHASTIC_H

#include "radiator.hh"

namespace emission
{
  
  class Synchrotron_Stochastic : public Radiator
  {
  public:

    Synchrotron_Stochastic(double maxg);


    virtual double source_function(double index){return index;};
    
  protected:

    
    std::array <double, n_grid_points> source_function_eval, V2, V3, sf_v2_ratio, V2_3_ratio;
    

    // precomputed factors for all CC rounds
    double factor1, DT; 
    
  private:
    
    void chang_cooper();
    virtual void emission_internal_computations(){};

    
  };

 #endif
