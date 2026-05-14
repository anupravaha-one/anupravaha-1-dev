#include"APKCoeffcalculation_Potential.h"

int Coeffcalculation_Potential(Base *base)
{
	int iZone,numZones,index,total_cells;
	Zone *zone;
	int *ptrSize,i,j,NI,NJ,NK;
	double *volumeDomain;
	double delt;
	
	DiffusionCoeffCalculation(base);
	
	numZones=base->Number_of_zones;
	delt = base->BaseParameters->TimeStep;
	for(iZone=0;iZone<numZones;iZone++)
	{
		zone=base->CgnsZone+iZone;
		i = 0;
		ptrSize = (zone->grid+i)->GSize;
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		total_cells = (NI-1) * (NJ-1) *(NK-1);
		volumeDomain = zone->ZGeom->Volume;
		
		(zone->grid+i)->DiffusionCoeff_Poten=(double *)malloc((NI-1)*(NJ-1)*(NK-1)*sizeof(double));
		zone->DeferredCorrCoeffDummy_P = (double *)malloc((NI-1)*(NJ-1)*(NK-1)*sizeof(double));
		if((zone->grid+i)->DiffusionCoeff_Poten==NULL)
		{
			printf("Memory allocation failure in APK SolverElectric.c\nExitign....\n");
			exit(1);
		}
		for(j=0;j<18;j++)
		{
			(zone->grid+i)->DiffusionCoeff_Poten_nb[j]=(double *)malloc((NI-1)*(NJ-1)*(NK-1)*sizeof(double));
			
			if((zone->grid+i)->DiffusionCoeff_Poten_nb[j]==NULL)
			{
				printf("Memory allocation failure in APK SolverElectric.c\nExitign....\n");
				exit(1);
			}
		}
		
		
		for (index=0;index<total_cells;index++)
		{
			*((zone->grid+i)->DiffusionCoeff_Poten+index)=*((zone->grid+i)->DiffusionCoeff_P+index)*delt*-1.0 + *(volumeDomain+index);
			*(zone->DeferredCorrCoeffDummy_P+index)=0.0;
			for(j=0;j<18;j++)
			{
				*((zone->grid+i)->DiffusionCoeff_Poten_nb[j]+index)=*((zone->grid+i)->DiffusionCoeff_nb[j]+index)*-1.0*delt;
				
			}
		}
	}
	return (0);	
}



