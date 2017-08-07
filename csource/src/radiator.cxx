#include "radiator.hh"


namespace emission{ 
    
  Radiator::Radiator() {}
  
  Radiator::Radiator(double maxg){
    
    //    n_grid_points = ngp;
    MAX_GAMMA = maxg;

    emitted = false;
    kernel_exists = false;
    
  }

  Radiator::~Radiator() {

    int i;
    
    for (i=0;i<n_photon_energies;i++)
      {
	delete[] synchrotron_kernel[i];
      }
    
    delete[] synchrotron_kernel;
    
  }
  

  void Radiator::initialize_electrons() {
    
    // iterators
    

    double step, step_plus_one;
    
    // define the step size such that we have a grid slightly
    // aout gamma max to avoid boundary issues

  
  
    step = exp(1./n_grid_points*log(MAX_GAMMA));
    step_plus_one = step + 1.;
    
    // initalize the grid of electrons
    
    fgamma.fill(0.);
    
    for(std::size_t i = 0; i != gamma.size(); ++i)
      {
	
	gamma[i]=pow(step,i);
	gamma2[i]=gamma[i]*gamma[i];
	
	if(i<n_grid_points-1)
	  {
	    G[i]=0.5*gamma[i]*step_plus_one;
	  }
	else
	  {
	    G[n_grid_points-1]=0.5*gamma[i]*step_plus_one;
	  }
	

	
      }


    
    

    

    
    
  }
  

  std::vector<double> Radiator::emission(double * energy, int n_photon_energies, int steps)
  {
      // precompute the synchrotron kernel

    int i;

    // FILL THE RADIATION
    radiation.reserve(n_photon_energies);
    
    if (!kernel_exists)
      {
	compute_synchrotron_kernel(energy);
      }

    
  /*
    Define a factor that is dependent on the magnetic field B.
    The cooling time of an electron at energy gamma is 

    DT = 6 * pi * me*c / (sigma_T B^2 gamma)
    DT = (1.29234E-9 B^2 gamma)^-1 seconds 
    DT = (cool * gamma)^-1 seconds

    where B is in Gauss. 
   */



    emission_internal_computations();
  
    cool = 1.29234E-9*B*B;
    
  

    // precompute gamma*cool and gamm2*cool
    for(std::size_t j = 0; j != gamma.size(); ++j)
      {
	coolg[j] = gamma[j]*cool;
	coolg2[j] = gamma2[j]*cool;
	
      }
    
    
    
    
    // now set the starting grid point for the synchrotron integration
    
    // now cool the electrons
    
    for(i=0;i<=steps;i++)
      {
	
	// run C&C for this round
	chang_cooper();
	
	// compute the emitted spectrum for this round
	// and add it to the total radiation
	synchrotron_spectrum(energy);
	
	// now fgamma is  cooled in the round	
      }
    
    return radiation;
  }
  
  void Radiator::reset()
  {
    // should reset all arrays

  }

  
  void Radiator::synchrotron_spectrum(double * energy)
  {
    /*
      compute the synchrotron emission spectrum in phts/(s keV cm^2) for an array of electron energies
      
      arguments:
      
      - energies: array of energies to compute spectrum at
      - B: magnetic field strength
      - gamma: electron energy
      - gamma2: electron energy squared
      - fgamma: electron distribution
      - n_grid_points: the number of electron grid points
      - n_photon_energies: number of photon energies
      
    */  
    
    //std::vector<double> out_val(n_photon_energies);
    double h,s, y[n_grid_points];
    double sum;
    int i,j,k;
    
    // integrate (convolve) the synchrotron kernel
    // with the electron distribution via trap. rule
    
    
    for(i=0;i<n_photon_energies;i++)
      {
	for(j=1;j<n_grid_points;j++)
	  {
	    y[j] = synchrotron_kernel[i][j] * fgamma[j] * (gamma[j] - gamma[j-1]);
	  }
	
	sum=0.;
	
	for(k=1;k<n_grid_points;k++)
	  {
	    sum+=y[k];
	  }
	
	// add this to the total emission
	radiation[i] += sum/(2.* energy[i]);
	
      }

    
    
  }

  void Radiator::compute_synchrotron_kernel(double * energy)
  {
    
    double h,s, ec, y[n_grid_points];
    double syncArg, sum;
    int i,j,k, status;
    
    
    // initialize the out matrix pointers
    
    
    synchrotron_kernel = new double *[n_photon_energies];
    for (i=0;i<n_photon_energies;i++)
      {
	synchrotron_kernel[i] = new double[n_grid_points];
      }
    
    
    gsl_set_error_handler_off();
    
    ec = 1.5*B/Bcritical; //Put proper units in!
    
    // integrate (convolve) the synchrotron kernel
    // with the electron distribution via trap. rule
    
    for(i=0;i<n_photon_energies;i++)
      {
	for(j=0;j<n_grid_points;j++)
	  {
	    syncArg = energy[i]/(ec * gamma2[j] );
	    status = gsl_sf_synchrotron_1_e(syncArg, &result);
	    
	    if (status == GSL_ERANGE)
	      {
		synchrotron_kernel[i][j]=0.;
		
	      }
	    else
	      {
		synchrotron_kernel[i][j] = result.val;	
	      }
	    
	  }
	
      }
    
    kernel_exists = true;
    
  }
  
  
  
  
  
  
  // end namespace
}
