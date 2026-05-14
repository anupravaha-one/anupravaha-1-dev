#include "APKCalcNewVolumefluxRestart.h"
/**
The function CalcNewVolumefluxRestart() recalculates the volumeflux from stored values of massvelocities
when the program is restarted. It uses the old time step "delta_t", which is stored in the file
"old_dt.dat" when the program is stopped. This function works for both 2D and 3D geometries.
*/
int CalcNewVolumefluxRestart(Base *base)
{
	int numzones,cntZ,*ptrSize,NI,NJ,NK;
	int i,j,k,index[3],indexP,iface,flagrho;
	Zone *zone;
	double pressureDiff[6],diffdia,*pressureDomain,*volumeDomain,*beta;
	double *slabPressure[6],*slabVolume[6];
	double coefn[6],delta_t=0.01,rho;
	FILE *mrt;
	numzones= base->Number_of_zones;
	mrt=fopen("Settings/old_dt.dat","r");
	fscanf(mrt,"%lf",&delta_t);
	fclose(mrt);
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
			flagrho = zone->Flowsol->flagRho;
			if(base->BaseParameters->flag2d==0)
			{
				for(k=0;k<NK-1;k++)   /*Storing current solution in Previous solution*/
				for(j=0;j<NJ-1;j++)
				for(i=0;i<NI-1;i++)
				{
					indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;

					index[0]=k,index[1]=j,index[2]=i;

					rho = *(zone->PrevFlowsol->Rho+flagrho*indexP);
					diffdia = Diffusion(index,pressureDomain,volumeDomain,ptrSize,slabPressure,slabVolume,beta,pressureDiff);
					for(iface=0;iface<6;iface++)
						*(zone->Flowsol->Massflux+indexP*6+iface)=*(zone->Flowsol->Massflux+indexP*6+iface) - (delta_t/rho)*pressureDiff[iface];
				}
			}
			if(base->BaseParameters->flag2d==1)
			{
				for(k=0;k<NK-1;k++)   /*Storing current solution in Previous solution*/
				for(j=0;j<NJ-1;j++)
				for(i=0;i<NI-1;i++)
				{
					indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;

					index[0]=k,index[1]=j,index[2]=i;

					rho = *(zone->PrevFlowsol->Rho+flagrho*indexP);
					diffdia = Diffusion2d(index,pressureDomain,volumeDomain,ptrSize,slabPressure,slabVolume,beta,pressureDiff,coefn);
					for(iface=0;iface<4;iface++)
					{
						*(zone->Flowsol->Massflux+indexP*6+iface)=*(zone->Flowsol->Massflux+indexP*6+iface)- (delta_t/rho)*pressureDiff[iface];
					}
				}
			}
		}
	}
	
	return 0;
}
