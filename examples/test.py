#!/Users/jburgess/.environs/3ml/bin/python
# Scientific libraries
import numpy as np

# Graphic libraries


from electradpy import sic_glue
ene = np.logspace(1,4,20)




sic = sic_glue.pySIC_Powerlaw(1.E8)

sic.set_power_law(1.E5, 1.E5, 1.E8,2.2)

phts = sic.photons(ene,1.E2,5)


print phts



del sic
