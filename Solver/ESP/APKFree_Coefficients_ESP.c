#include "APKFree_Coefficients_ESP.h"

int Free_Coefficients_ESP(Base *base)
{	

	int i,j;
	int numZones, nZones;
	Zone *zone;

	nZones=(base->Number_of_zones);
	i=0;
	
	for(numZones=0;numZones<nZones;numZones++)
	{
		zone=base->CgnsZone+numZones;
		for(j=0;j<12;j++)
			free(zone->DeferredCorrCoeffPhi_nb[j]);
		
		for(j=0;j<18;j++)
		{
			free((zone->grid+i)->DiffusionCoeffPhi_nb[j]);
			free((zone->grid+i)->ConvDiffusionCoeffPhi_nb[j]);
			free((zone->grid+i)->DiffusionCoeff_Poten_nb[j]);
			free((zone->grid+i)->DiffusionCoeff_nb[j]);
		}
		for(j=0;j<6;j++)
		free((zone->grid+i)->Conv_Coeff_nb[j]);
		
		free((zone->grid+i)->DiffusionCoeff_P);
		free((zone->grid+i)->ConvDiffusionCoeffPhi_P);
		free((zone->grid+i)->DiffusionCoeffPhi_P);
		free((zone->grid+i)->DiffusionCoeff_Poten);
		free((zone->grid+i)->Conv_Coeff_P);
		free(zone->DeferredCorrCoeffPhi_P);
		for(j=0;j<6;j++)
			free((zone->grid+i)->Diffusion_Diag[j]);
		for(j=0;j<36;j++)
			free((zone->grid+i)->Diffusion_nb[j]);
		for(j=0;j<24;j++)
			free((zone->grid+i)->Diffusion_Nonortho[j]);
	}
	return (0);
}

