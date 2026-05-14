#include "APKCalcNewVolumeFluxDisperse.h"
/**
The function CalcNewVolumeFluxDisperse() calculates the volumeflux of Dispersed Phase for 3 dimensional geometry from the current computed values of Pressure and the old volumeflux. Similarly the function CalcNewVolumeFluxDisperse2d() is used for 2 dimensional geometries.
*/

int CalcNewVolumeFluxDisperse(Base *base)
{
	int numzones,cntZ,*ptrSize,NI,NJ,NK;
	int i,j,k,index[3],indexP,iface,flagrhoDisperse;
	Zone *zone;
	double pressureDiff[6],diffdia,*pressureDomain,*volumeDomain,*beta;
	double *slabPressure[6],*slabVolume[6];
	double delta_t,rhoDisperse,*diffusionCoeffdomain;
	
	delta_t = base->BaseParameters->TimeStep;
	numzones= base->Number_of_zones;
	
	for(cntZ=0;cntZ<numzones;cntZ++)
	{
		zone = (base->CgnsZone+cntZ);
		
		ptrSize = zone->Size;
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);

		beta = zone->ZGeom->Beta;
		pressureDomain=MakeDomainArray(zone,4);
		volumeDomain=MakeDomainArray(zone,5);

		MakeSlabPtrArray(zone,slabPressure,4);
		MakeSlabPtrArray(zone,slabVolume,5);
		diffusionCoeffdomain = MakeDiffusionCoeffDomain(zone,4);
		flagrhoDisperse = base->BaseParameters->flagRhoDisperse;
		
		indexP=0;
		for(k=0;k<NK-1;k++)  
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{		
			index[0]=k,index[1]=j,index[2]=i;
								
			rhoDisperse = *(zone->Flowsol->RhoDisperse+indexP*flagrhoDisperse);
					
			diffdia = Diffusion3dWithConstDiffCoeff_opti(index,pressureDomain,volumeDomain,ptrSize,slabPressure,slabVolume,diffusionCoeffdomain,beta,pressureDiff);
		
			for(iface=0;iface<6;iface++)
				*(zone->Flowsol->MassfluxDisperse+indexP*6+iface)=*(zone->Flowsol->MassfluxDisperse+indexP*6+iface) - (delta_t/rhoDisperse)*pressureDiff[iface];
					
			indexP++;
		}
		
	}
	return 0;
}

int CalcNewVolumeFluxDisperse2d(Base *base)
{
	int numzones,cntZ,*ptrSize,NI,NJ,NK;
	int i,j,k,index[3],indexP,iface,flagrhoDisperse;
	Zone *zone;
	double pressureDiff[6],diffdia,*pressureDomain,*volumeDomain,*beta;
	double *slabPressure[6],*slabVolume[6];
	double delta_t,rhoDisperse,*diffusionCoeffdomain;
	
	delta_t = base->BaseParameters->TimeStep;
	numzones= base->Number_of_zones;
	
	for(cntZ=0;cntZ<numzones;cntZ++)
	{
		zone = (base->CgnsZone+cntZ);
		
		ptrSize = zone->Size;
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);

		beta = zone->ZGeom->Beta;
		pressureDomain=MakeDomainArray(zone,4);
		volumeDomain=MakeDomainArray(zone,5);

		MakeSlabPtrArray(zone,slabPressure,4);
		MakeSlabPtrArray(zone,slabVolume,5);
		diffusionCoeffdomain = MakeDiffusionCoeffDomain(zone,4);
		flagrhoDisperse = base->BaseParameters->flagRhoDisperse;
		
		indexP=0;
		for(k=0;k<NK-1;k++)  
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{		
			index[0]=k,index[1]=j,index[2]=i;
								
			rhoDisperse = *(zone->Flowsol->RhoDisperse+indexP*flagrhoDisperse);
					
			diffdia = Diffusion3dWithConstDiffCoeff_opti(index,pressureDomain,volumeDomain,ptrSize,slabPressure,slabVolume,diffusionCoeffdomain,beta,pressureDiff);
		
			for(iface=0;iface<4;iface++)
			{
				*(zone->Flowsol->MassfluxDisperse+indexP*6+iface)=*(zone->Flowsol->MassfluxDisperse+indexP*6+iface) - (delta_t/rhoDisperse)*pressureDiff[iface];
			/*	if(i==(int)NI/2 && j==(int)NJ/2)
				{
					printf("\nFor%d,%d:\n",i,j);
					printf("total volume flux of dispersed phase for face = %d = %e\n\n",iface,*(zone->Flowsol->MassfluxDisperse+indexP*6+iface));
				}	*/
			}		
			indexP++;
		}
		
	}
	return 0;
}

