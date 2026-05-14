#include "APKDiffusionCoeffCalculation.h"



int DiffusionCoeffCalculation(Base *base)
{
	Zone *zone;
	int *ptrSize,i,j,NI,NJ,NK;
	int iZone,numZones;
	numZones=base->Number_of_zones;
//	int k,g,INDEX;
	
	if((base->MultiGridflag)==1)
	{
		for(iZone=0;iZone<numZones;iZone++)
		{
			zone=base->CgnsZone+iZone;
			for(i=0;i<=(zone->zonalCoarseGridIndex);i++)
			{
				ptrSize = (zone->grid+i)->GSize;
				NI=*(ptrSize+0);
				NJ=*(ptrSize+1);
				NK=*(ptrSize+2);
				(zone->grid+i)->DiffusionCoeff_P=(double *)malloc((NI-1)*(NJ-1)*(NK-1)*sizeof(double));
				for(j=0;j<18;j++)
				{
					(zone->grid+i)->DiffusionCoeff_nb[j]=(double *)malloc((NI-1)*(NJ-1)*(NK-1)*sizeof(double));
				}
			}
		}
		
		for(iZone=0;iZone<numZones;iZone++)
		{
			zone=base->CgnsZone+iZone;
			for(i=0;i<=(zone->zonalCoarseGridIndex);i++)
			{
				if(i == 0)
					CallCalcCoefficients(base,zone,i);
				else
					CoarseGridCoefficientsDiffusion(base,zone,i);
			}
		}
	}
        else if((base->MultiGridflag)==0)
	{
         
	
           	for(iZone=0;iZone<numZones;iZone++)
		{
			zone=base->CgnsZone+iZone;
			i = 0;
			ptrSize = (zone->grid+i)->GSize;
			NI=*(ptrSize+0);
			NJ=*(ptrSize+1);
			NK=*(ptrSize+2);
			(zone->grid+i)->DiffusionCoeff_P=(double *)malloc((NI-1)*(NJ-1)*(NK-1)*sizeof(double));
			for(j=0;j<18;j++)
			{
				(zone->grid+i)->DiffusionCoeff_nb[j]=(double *)malloc((NI-1)*(NJ-1)*(NK-1)*sizeof(double));
			}
		}
                 

		for(iZone=0;iZone<numZones;iZone++)
		{
                   	
			zone=base->CgnsZone+iZone;
			i = 0;
			CallCalcCoefficients_opti(base,zone,i);
		}
	}
		
	return 0;
}
