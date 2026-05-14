#include "APKCalcNewVolumeFlux.h"
/**
The function CalcNewVolumeFlux() calculates the volumeflux for 3 dimensional geometry
from the current computed values of Pressure and the old volumeflux.
Similarly the function CalcNewVolumeFlux2d() is used for 2 dimensional geometries.
See Equation(3.39) in Project Report No.1, Dec 2005.
*/


int CalcNewVolumeFlux(Base *base)
{
	int numzones,cntZ,*ptrSize,NI,NJ,NK;
	int i,j,k,index[3],indexP,iface,flagrho;
	Zone *zone;
	double pressureDiff[6],diffdia,*pressureDomain,*volumeDomain,*beta;
	double *slabPressure[6],*slabVolume[6];
	double delta_t,rho,*diffusionCoeffdomain;
	
	delta_t = base->BaseParameters->TimeStep;
	numzones= base->Number_of_zones;
	
	for(cntZ=0;cntZ<numzones;cntZ++)
	{
		zone = (base->CgnsZone+cntZ);
		if(zone->flagSolidFluid==0)
		{
			ptrSize = zone->Size;
			NI=*(ptrSize+0);
			NJ=*(ptrSize+1);
			NK=*(ptrSize+2);
			beta = zone->ZGeom->Beta;
			pressureDomain=MakeDomainArray(zone,4); //variableDomain stores variable pointer
			volumeDomain=MakeDomainArray(zone,5);

			MakeSlabPtrArray(zone,slabPressure,4);/*Called for making slabs of Pressure*/
			MakeSlabPtrArray(zone,slabVolume,5);/*Called for making slabs of Volumes*/
			diffusionCoeffdomain = MakeDiffusionCoeffDomain(zone,4);
			flagrho = zone->Flowsol->flagRho;
			
			if(base->flagMPFV==0)
			{	
				indexP=0;
				for(k=0;k<NK-1;k++)   /*Storing current solution in Previous solution*/
				for(j=0;j<NJ-1;j++)
				for(i=0;i<NI-1;i++)
				{
					
	
					index[0]=k,index[1]=j,index[2]=i;
					
					
					rho = *(zone->Flowsol->Rho+indexP*flagrho);
					
					diffdia = Diffusion3dWithConstDiffCoeff_opti(index,pressureDomain,volumeDomain,ptrSize,slabPressure,slabVolume,diffusionCoeffdomain,beta,pressureDiff);
				
					
					for(iface=0;iface<6;iface++)
						*(zone->Flowsol->Massflux+indexP*6+iface)=*(zone->Flowsol->Massflux+indexP*6+iface) - (delta_t/rho)*pressureDiff[iface];
					
					indexP++;
				}
			}
			else if(base->flagMPFV==1)
			{
				indexP=0;
				for(k=0;k<NK-1;k++)   /*Storing current solution in Previous solution*/
				for(j=0;j<NJ-1;j++)
				for(i=0;i<NI-1;i++)
				{
					
	
					index[0]=k,index[1]=j,index[2]=i;

					rho = 1.0;
					
					diffdia = Diffusion3dWithConstDiffCoeff_optiMPFV(index,zone,pressureDomain,volumeDomain,ptrSize,slabPressure,slabVolume,beta,pressureDiff);
					

					for(iface=0;iface<6;iface++)
						*(zone->Flowsol->Massflux+indexP*6+iface)=*(zone->Flowsol->Massflux+indexP*6+iface) - (delta_t/rho)*pressureDiff[iface];
					
					indexP++;
				}
			}
				
			
			
		}
	}
	return 0;
}

int CalcNewVolumeFlux2d(Base *base)
{
	int numzones,cntZ,*ptrSize,NI,NJ,NK;
	int i,j,k,index[3],indexP,iface,flagrho;
	Zone *zone;
	double pressureDiff[6],diffdia,*pressureDomain,*volumeDomain,*beta;
	double *slabPressure[6],*slabVolume[6],CoeffN[6];
	double delta_t,rho,*diffusionCoeffdomain;
	
	delta_t = base->BaseParameters->TimeStep;
	numzones= base->Number_of_zones;
	
	for(cntZ=0;cntZ<numzones;cntZ++)
	{
		zone = (base->CgnsZone+cntZ);
		if(zone->flagSolidFluid==0)
		{
			ptrSize = zone->Size;
			NI=*(ptrSize+0);
			NJ=*(ptrSize+1);
			NK=*(ptrSize+2);
			beta = zone->ZGeom->Beta;
			pressureDomain=MakeDomainArray(zone,4); //variableDomain stores variable pointer
			volumeDomain=MakeDomainArray(zone,5);

			MakeSlabPtrArray(zone,slabPressure,4);/*Called for making slabs of Pressure*/
			MakeSlabPtrArray(zone,slabVolume,5);/*Called for making slabs of Volumes*/
			diffusionCoeffdomain = MakeDiffusionCoeffDomain(zone,4);
			flagrho = zone->Flowsol->flagRho;
			
			if(base->flagMPFV==0)
			{	
				indexP=0;
				for(k=0;k<NK-1;k++)   /*Storing current solution in Previous solution*/
				for(j=0;j<NJ-1;j++)
				for(i=0;i<NI-1;i++)
				{
					
	
					index[0]=k,index[1]=j,index[2]=i;
					
					
					rho = *(zone->Flowsol->Rho+indexP*flagrho);
	
					diffdia = Diffusion2dWithConstDiffCoeff_opti(index,pressureDomain,volumeDomain,ptrSize,slabPressure,slabVolume,diffusionCoeffdomain,beta,pressureDiff,CoeffN);
					
					
					for(iface=0;iface<4;iface++)
					{
						*(zone->Flowsol->Massflux+indexP*6+iface)=*(zone->Flowsol->Massflux+indexP*6+iface) - (delta_t/rho)*pressureDiff[iface];
					//fprintf(value,"\nvalue[%d]=%e",iface,*(zone->Flowsol->Massflux+indexP*6+iface));
					}
					indexP++;
				}	
			}
			else if(base->flagMPFV==1)
			{	
				indexP=0;
				for(k=0;k<NK-1;k++)   /*Storing current solution in Previous solution*/
				for(j=0;j<NJ-1;j++)
				for(i=0;i<NI-1;i++)
				{
					
	
					index[0]=k,index[1]=j,index[2]=i;

					rho = 1.0;
	
					diffdia = Diffusion2dWithConstDiffCoeff_optiMPFV(index,zone,pressureDomain,volumeDomain,ptrSize,slabPressure,slabVolume,beta,pressureDiff,CoeffN);
				
					for(iface=0;iface<4;iface++)
					{
						*(zone->Flowsol->Massflux+indexP*6+iface)=*(zone->Flowsol->Massflux+indexP*6+iface) - (delta_t/rho)*pressureDiff[iface];
					//fprintf(value,"\nvalue[%d]=%e",iface,*(zone->Flowsol->Massflux+indexP*6+iface));
					}
					indexP++;
				}			
			}
		}
	}
	
	return 0;
}
