#include "APKCDcoeffMemAlloc.h"

/**
It allocates memory to convection-diffusion equation coefficients and deffered correction coefficients for all the
variables i.e., velocities, pressure, temperature, scalars.
*/

int CDcoeffMemAlloc(Base *base)
{

	int *ptrSize,nZones,iZone,NI,NJ,NK,i,j;
	Zone *zone;
	nZones=(base->Number_of_zones);

	i=0;
	for(iZone=0;iZone<nZones;iZone++)
	{
		zone=base->CgnsZone+iZone;
		ptrSize = (zone->grid+0)->GSize;
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		//zone->DeferredCorrCoeffPhi_P = (double *)malloc((NI-1)*(NJ-1)*(NK-1)*sizeof(double));
		zone->DeferredCorrCoeffPhi_P = (double *)calloc((NI-1)*(NJ-1)*(NK-1),sizeof(double));
		for(j=0;j<12;j++)
		{
			//zone->DeferredCorrCoeffPhi_nb[j] = (double *)malloc((NI-1)*(NJ-1)*(NK-1)*sizeof(double));
			zone->DeferredCorrCoeffPhi_nb[j] = (double *)calloc((NI-1)*(NJ-1)*(NK-1),sizeof(double));
		}
		
		///Diffusion coefficient memory allocation is also done here
		(zone->grid+i)->DiffusionCoeffPhi_P=(double *)calloc((NI-1)*(NJ-1)*(NK-1),sizeof(double));
		for(j=0;j<18;j++)
		{
			(zone->grid+i)->DiffusionCoeffPhi_nb[j]=(double *)calloc((NI-1)*(NJ-1)*(NK-1),sizeof(double));
		}
		//double *Diffusion_Daig[6],*Diffusion_nb[36],*Diffusion_Nonortho[24];
		///This Memory allocation is done to improve the calculations in case of cases with varying diffusion coefficient problems.

		for(j=0;j<6;j++)
		{
			(zone->grid+i)->Diffusion_Diag[j]=	(double *)calloc((NI-1)*(NJ-1)*(NK-1),sizeof(double));
			if((zone->grid+i)->Diffusion_Diag[j]==NULL)
			{
				printf("Memory can not be allocated in APKCDcoeffMemAlloc.c\nExiting.....\n");
				exit(76);
			}
		}
		for(j=0;j<36;j++)
		{
			(zone->grid+i)->Diffusion_nb[j]=	(double *)calloc((NI-1)*(NJ-1)*(NK-1),sizeof(double));
			if((zone->grid+i)->Diffusion_nb[j]==NULL)
			{
				printf("Memory can not be allocated in APKCDcoeffMemAlloc.c\nExiting.....\n");
				exit(76);
			}
		}
		for(j=0;j<24;j++)
		{
			(zone->grid+i)->Diffusion_Nonortho[j]=	(double *)calloc((NI-1)*(NJ-1)*(NK-1),sizeof(double));
			if((zone->grid+i)->Diffusion_Nonortho[j]==NULL)
			{
				printf("Memory can not be allocated in APKCDcoeffMemAlloc.c\nExiting.....\n");
				exit(76);
			}
		}

		
		///Convective coefficient memory allocation is also done here
		(zone->grid+i)->Conv_Coeff_P=(double *)calloc((NI-1)*(NJ-1)*(NK-1),sizeof(double));
		for(j=0;j<6;j++)
		{
			(zone->grid+i)->Conv_Coeff_nb[j]=(double *)calloc((NI-1)*(NJ-1)*(NK-1),sizeof(double));
		}
		
		
		
		if((zone->Flowsol->PotentialMHD!=NULL)|| (base->flagMPFV==1)||(base->flagDroplet==1)||(base->flagParticle==1))
			zone->DeferredCorrCoeffDummy_P = (double *)calloc((NI-1)*(NJ-1)*(NK-1),sizeof(double));
		
		
	}
	
	for(iZone=0;iZone<nZones;iZone++)
	{
		zone=base->CgnsZone+iZone;
		i = 0;
		ptrSize = (zone->grid+i)->GSize;
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		(zone->grid+i)->ConvDiffusionCoeffPhi_P=(double *)calloc((NI-1)*(NJ-1)*(NK-1),sizeof(double));
		for(j=0;j<18;j++)
		{
			(zone->grid+i)->ConvDiffusionCoeffPhi_nb[j]=(double *)calloc((NI-1)*(NJ-1)*(NK-1),sizeof(double));
		}
	}
	
	if((base->MultiGridflag)==1)
	{
		for(iZone=0;iZone<nZones;iZone++)
		{
			zone=base->CgnsZone+iZone;
			for(i=1;i<=(zone->zonalCoarseGridIndex);i++)
			{
				ptrSize = (zone->grid+i)->GSize;
				NI=*(ptrSize+0);
				NJ=*(ptrSize+1);
				NK=*(ptrSize+2);
				(zone->grid+i)->ConvDiffusionCoeffPhi_P=(double *)calloc((NI-1)*(NJ-1)*(NK-1),sizeof(double));
				for(j=0;j<18;j++)
				{
					(zone->grid+i)->ConvDiffusionCoeffPhi_nb[j]=(double *)calloc((NI-1)*(NJ-1)*(NK-1),sizeof(double));
				}
			}
		}
	}
	return 0;
}


