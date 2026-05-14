#include "APKPressureIterArtificial.h"

/**
This function solves for pressure in the current iteration after the mass velocities have been
converged. This function loops over all the zones and solves for pressure until the total residue
is less than the specified value. This function actually implements the artificial compressibility
method. This function inturn calls "UpdatePressureArtificial" or "UpdatePressureArtificial2D"
which actually does the G-S iterations.
*/

int PressureIterArtificial(Base *base)
{
	FILE *pressConvcri;
	Zone *zone;
	double residBase,residZone,*volumeDomain;
	int g,i,j,k,NI,NJ,NK,*ptrSize,indexP,iFace,flagrho;
	double *massFlux,sumMassFlux,betatau=1.0,rho,*prevpressureDomain;
	int totalBaseCells=0,iZone,numZones;
	numZones=base->Number_of_zones;
	totalBaseCells = base->TotalCellFluid;
	
	
	for(iZone=0;iZone<numZones;iZone++)
	{
		zone=base->CgnsZone+iZone;
		g = 0;
		ptrSize = (zone->grid+g)->GSize;
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		(zone->grid+g)->GRHS=(double *)malloc((NI-1)*(NJ-1)*(NK-1)*sizeof(double));
	}
	
	if(base->BaseParameters->flag2d == 1)
	{
		for(iZone=0;iZone<numZones;iZone++)
		{
			zone=base->CgnsZone+iZone;
			volumeDomain = zone->ZGeom->Volume;
			flagrho = zone->Flowsol->flagRho;
			prevpressureDomain = zone->PrevFlowsol->P;
			ptrSize = (zone->grid+0)->GSize;
			NI=*(ptrSize+0);NJ=*(ptrSize+1);NK=*(ptrSize+2);
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				rho = *(zone->Flowsol->Rho+indexP*flagrho);
				
				massFlux=zone->Flowsol->Massflux+indexP*6;
				sumMassFlux=0.0;
				for(iFace=0;iFace<4;iFace++)
				{
					sumMassFlux=sumMassFlux +(*(massFlux+iFace));
				}
				*((zone->grid)->GRHS+indexP) = (*(volumeDomain+indexP)/(rho*betatau))*(*(prevpressureDomain+indexP))-sumMassFlux;
			}
		}
		do
		{
			pressConvcri=fopen("Settings/pressureconvergCriteria.dat","r");
			fscanf(pressConvcri,"%lf",&(base->BaseParameters->PressLoopConv));
			fclose(pressConvcri);

			residBase=0.0;
			for(iZone=0;iZone<numZones;iZone++)
			{
				zone=base->CgnsZone+iZone;
				if(zone->flagSolidFluid==0)
				{
					ExtractVariables(4,base,zone);
					
					residZone=UpdatePressureArtificial2d(base,zone);
					
					residBase=residBase+residZone;
				}
			}
			residBase=sqrt(residBase/totalBaseCells);
			//printf("residual pressure=%12.8e\n",residBase);
		}
		while(residBase > base->BaseParameters->PressLoopConv);
	}
	else
	{
		for(iZone=0;iZone<numZones;iZone++)
		{
			zone=base->CgnsZone+iZone;
			volumeDomain = zone->ZGeom->Volume;
			flagrho = zone->Flowsol->flagRho;
			prevpressureDomain = zone->PrevFlowsol->P;
			ptrSize = (zone->grid+0)->GSize;
			NI=*(ptrSize+0);NJ=*(ptrSize+1);NK=*(ptrSize+2);
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
				rho = *(zone->Flowsol->Rho+indexP*flagrho);
				
				massFlux=zone->Flowsol->Massflux+indexP*6;
				sumMassFlux=0.0;
				for(iFace=0;iFace<6;iFace++)
				{
					sumMassFlux=sumMassFlux +(*(massFlux+iFace));
				}
				*((zone->grid)->GRHS+indexP) = (*(volumeDomain+indexP)/(rho*betatau))*(*(prevpressureDomain+indexP))-sumMassFlux;
			}
		}
		do
		{
			pressConvcri=fopen("Settings/pressureconvergCriteria.dat","r");
			fscanf(pressConvcri,"%lf",&(base->BaseParameters->PressLoopConv));
			fclose(pressConvcri);

			residBase=0.0;
			for(iZone=0;iZone<numZones;iZone++)
			{
				zone=base->CgnsZone+iZone;
				if(zone->flagSolidFluid==0)
				{
					ExtractVariables(4,base,zone);
					
					residZone=UpdatePressureArtificial(base,zone);
					
					residBase=residBase+residZone;
				}
			}
			residBase=sqrt(residBase/totalBaseCells);
			//printf("residual pressure=%12.8e\n",residBase);
		}
		while(residBase > base->BaseParameters->PressLoopConv);
	}
	
	RefreshSurfaceCells(4,base);
	
	return 0;
}
