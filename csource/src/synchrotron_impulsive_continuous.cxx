#include "synchrotron_impulsive_continuous.hh"
#include <iostream>

namespace emission {

  //  Synchrotron_Impulsive_Continuous::Synchrotron_Impulsive_Continuous() {}
  
  Synchrotron_Impulsive_Continuous::Synchrotron_Impulsive_Continuous(double maxg):Radiator(maxg){}

  void Synchrotron_Impulsive_Continuous::chang_cooper()
  {
    double delta_gamma;
    
    
    int j;
    
    
    // We will walk backwards so we must set the end points
    
    fgammatp1[n_grid_points-1]=fgamma[n_grid_points-1]*factor1;
    
    for(j = n_grid_points-2; j>=1; j--)
      {
	// we have precomputed the contant factors
	fgammatp1[j] = fgamma[j]/V2[j] + sf_v2_ratio[j] + fgammatp1[j+1] * V2_3_ratio[j];




	//	std::cout<<fgammatp1[j] <<"\n";
	
      }
    
    // Set the end point 
    delta_gamma = G[1]-G[0];
    fgammatp1[0]=fgamma[0]+(DT* coolg2[1]* fgammatp1[1])/ delta_gamma;
    
    
    // Set fgamma to the soultion of the energy
    // equation from this run.
  
    std::copy(std::begin(fgammatp1),std::end(fgammatp1),std::begin(fgamma) );
    // for(j=0;j<n_grid_points;j++)
    //   {
    //     fgamma[j]=fgammatp1[j];
    //   }
 
  }

  //////////////////////////////////////////////////////////////////////////////////


  /////// SIC POWER LAW

  //  SIC_PowerLaw::SIC_PowerLaw() {}
  
  SIC_PowerLaw::SIC_PowerLaw(double maxg) : Synchrotron_Impulsive_Continuous(maxg) {}

  

  void SIC_PowerLaw::emission_internal_computations()
  {

    double delta_gamma;
    int j;
    
    // compute the source function

    for(std::size_t i = 0; i != gamma.size(); ++i)
      {
	source_function_eval[i] =source_function(gamma[i], this_ne , this_gamma_min, this_gamma_max, this_p);
	
      }

    DT = cool/this_gamma_max;

    
    // compute terms for cc
    
    
    // factor 1
    delta_gamma = G[n_grid_points-1]-G[n_grid_points-2];
    
    factor1 = 1./(1. + (DT * cool * gamma2[n_grid_points-1] )/ delta_gamma);
  
    for(j = n_grid_points-2; j>=1; j--)
      {
      
	delta_gamma = 0.5*(G[j]-G[j-1]); //Half steps are at j+.5 and j-.5
      
	// gdotp=coolg2[j+1]; // Forward  step cooling 
	// gdotm=coolg2[j];     // Backward step cooling

      
	V3[j] = (DT * cool*gamma2[j+1] )/ delta_gamma;     // Tridiagonal coeff.
	V2[j] = 1. + (DT*cool*gamma2[j])/ delta_gamma;// Tridiagonal coeff.

	// the ratio
	V2_3_ratio[j] = V3[j]/V2[j];

	sf_v2_ratio[j] = source_function_eval[j]*DT/V2[j];

      
	// Solve for forward electron distribution
    
      	  



      }
  
  
  

  }


  double SIC_PowerLaw::source_function(double x, double ne ,double gamma_min,double gamma_max, double p)
  {
    /*

      the source spectrum of the electrons. it is a power law
    
      arguments:

      - x: the energy grid
      - ne: the number of electrons
      - gamma_min: minimum electron injection energy
      - gamma_max: maximum electron energy

    */

  
    if((x<gamma_min) || (x>gamma_max)){

      return 0.;

    }

    else
      {
      
	return ne * pow(x,-p) * (1-p)* 1./ (pow(gamma_max,1-p)-pow(gamma_min,1-p));
      }
  
  
  }


}
