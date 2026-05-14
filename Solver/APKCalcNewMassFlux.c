#include "APKCalcNewMassFlux.h"
/**
The function CalcNewMassFlux() calculates the massflux for 3 dimensional geometry
from the current computed values of Pressure and the old Massflux.
See Equation(3.39) in Project Report No.1, Dec 2005.
*/
int CalcNewMassFlux(Base *base)
{
	int numzones,cntZ,*ptrSize,NI,NJ,NK;
	int i,j,k,index[3],indexP,iface;
	Zone *zone;
	double pressureDiff[6],diffdia,*pressureDomain,*volumeDomain,*beta;
	double *slabPressure[6],*slabVolume[6];
	double delta_t;
	
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
			
			for(k=0;k<NK-1;k++)   /*Storing current solution in Previous solution*/
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;

				index[0]=k,index[1]=j,index[2]=i;
				
				diffdia = Diffusion(index,pressureDomain,volumeDomain,ptrSize,slabPressure,slabVolume,beta,pressureDiff);
				for(iface=0;iface<6;iface++)
					*(zone->Flowsol->Massflux+indexP*6+iface)=*(zone->Flowsol->Massflux+indexP*6+iface) - (delta_t)*pressureDiff[iface];
			}
		}
	}
	return 0;
}

