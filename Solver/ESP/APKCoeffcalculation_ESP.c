#include "APKCoeffcalculation_ESP.h"

int Coeffcalculation_ESP(Base *base)
{	

//	int *ptrSize,NI,NJ,NK,i,j;
// 	int numZones, nZones;
// 	Zone *zone;
// 
// 	nZones=(base->Number_of_zones);
// 	i=0;
// 	
// 	for(numZones=0;numZones<nZones;numZones++)
// 	{
// 		zone=base->CgnsZone+numZones;
// 		ptrSize = (zone->grid+0)->GSize;
// 		NI=*(ptrSize+0);
// 		NJ=*(ptrSize+1);
// 		NK=*(ptrSize+2);
// 	
// 		if(zone->Flowsol->RhoIon!=NULL)
// 		{
// 			zone->DeferredCorrCoeffRhoIon_P = (double *)malloc((NI-1)*(NJ-1)*(NK-1)*sizeof(double));
// 			for(j=0;j<12;j++)
// 			{
// 				zone->DeferredCorrCoeffRhoIon_nb[j] = (double *)malloc((NI-1)*(NJ-1)*(NK-1)*sizeof(double));
// 			}
// 			
// 			(zone->grid+i)->DiffusionCoeffRhoIon_P=(double *)malloc((NI-1)*(NJ-1)*(NK-1)*sizeof(double));
// 			for(j=0;j<18;j++)
// 			{
// 				(zone->grid+i)->DiffusionCoeffRhoIon_nb[j]=(double *)malloc((NI-1)*(NJ-1)*(NK-1)*sizeof(double));
// 			}
// 			(zone->grid+i)->ConvDiffusionCoeffRhoIon_P=(double *)malloc((NI-1)*(NJ-1)*(NK-1)*sizeof(double));
// 			for(j=0;j<18;j++)
// 			{
// 			(zone->grid+i)->ConvDiffusionCoeffRhoIon_nb[j]=(double *)malloc((NI-1)*(NJ-1)*(NK-1)*sizeof(double));
// 			}
// 		}
// 	}
	
	CDcoeffMemAlloc(base);
	Coeffcalculation_Potential(base);
	

// 	if(base->CgnsZone->Flowsol->RhoIon!=NULL)
// 	{
// 		CallCalcCoefficientsForVariables(base,0,9);
// 		///Remember that Char ion equation is not solvable for Transient diffusion coefficient (i.e., diffusion coefficient varying with time).
// 	}

	return (0);
}

