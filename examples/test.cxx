#include "synchrotron_impulsive_continuous.hh"
#include "radiator.hh"



int main()
{
  double a= 10000.;
  emission::SIC_PowerLaw *test = new emission::SIC_PowerLaw(a);
  //emission::Radiator *test = new emission::Radiator(a);
  
  
  test->initialize_electrons();
  test->set_power_law(10000.,1000.,a,2.2);

    

  

  
  return 0;
}
