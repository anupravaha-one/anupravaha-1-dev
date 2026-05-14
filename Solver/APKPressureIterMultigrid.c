#include "APKPressureIterMultigrid.h"

/**
This function solves for pressure in the current iteration after the mass velocities have been
converged. This function loops over all the zones and solves for pressure until the total residue
is less than the specified value.This function calls other functions for solving pressure
for different cases of density (constant or variable) and geometry (2D or 3D) in which the actual
G-S iterations are executed.
*/

int PressureIterMultigrid(Base *base)
{
//	printf("I am into pressureIterMultigrid\n");
	
	FILE *pressConvcri;
	Zone *zone;
	int pressureItercase;
//	int pressureIterCount=0;
	int *ptrSize,i,j,NI,NJ,NK;
	double residBase,residZone;
	int totalBaseCells=0,iZone,numZones;
	numZones=base->Number_of_zones;
	totalBaseCells = base->TotalCellFluid;
	
//	FILE *vn;
//	vn = fopen("iterations.dat","a");
	
//	int l,a,g,k,indexP,bc_no;
	int itr,gdi;
	
	if(base->BaseParameters->flag2d == 1 && base->flagVariableDensityProblem == 0)
		pressureItercase = 0; /*case for calling UpdatePressure2d*/
	
	if(base->BaseParameters->flag2d == 0 && base->flagVariableDensityProblem == 0)
		pressureItercase = 1; /*case for calling UpdatePressure, this is for normal 3d*/
	
	if(base->BaseParameters->flag2d == 1 && base->flagVariableDensityProblem == 1)
		pressureItercase = 2; /*case for calling UpdatePressure2dVarDensity*/
	
	if(base->BaseParameters->flag2d == 0 && base->flagVariableDensityProblem == 1)
		pressureItercase = 3; /*case for calling UpdatePressureVarDensity, this is for 3d variable Density case*/
	
	for(iZone=0;iZone<numZones;iZone++)
	{
		zone=base->CgnsZone+iZone;
		for(i=0;i<=(zone->zonalCoarseGridIndex);i++)
		{
			ptrSize = (zone->grid+i)->GSize;
			NI=*(ptrSize+0);
			NJ=*(ptrSize+1);
			NK=*(ptrSize+2);
			(zone->grid+i)->Gresidual=(double *)malloc((NI-1)*(NJ-1)*(NK-1)*sizeof(double));
			(zone->grid+i)->E=(double *)malloc((NI-1)*(NJ-1)*(NK-1)*sizeof(double));
			(zone->grid+i)->GRHS=(double *)malloc((NI-1)*(NJ-1)*(NK-1)*sizeof(double));
		}
	}
	for(iZone=0;iZone<numZones;iZone++)
	{
		zone=base->CgnsZone+iZone;
		for(i=1;i<=(zone->zonalCoarseGridIndex);i++)
		{
			ptrSize = (zone->grid+i)->GSize;
			NI=*(ptrSize+0);
			NJ=*(ptrSize+1);
			NK=*(ptrSize+2);
			(zone->grid+i)->V=(double *)malloc((NI-1)*(NJ-1)*(NK-1)*sizeof(double));
			for(j=0;j<6;j++)
			{
				if(j==0||j==1)
				(zone->grid+i)->SlabV[j]=(double *)malloc((NK-1)*(NJ-1)*sizeof(double));
				else if(j==2||j==3)
				(zone->grid+i)->SlabV[j]=(double *)malloc((NK-1)*(NI-1)*sizeof(double));
				else
				(zone->grid+i)->SlabV[j]=(double *)malloc((NJ-1)*(NI-1)*sizeof(double));
			}
		}
		for(i=0;i<=(zone->zonalCoarseGridIndex);i++)
		{
			ptrSize = (zone->grid+i)->GSize;
			NI=*(ptrSize+0);
			NJ=*(ptrSize+1);
			NK=*(ptrSize+2);
			(zone->grid+i)->DummyV=(double *)malloc((NI-1)*(NJ-1)*(NK-1)*sizeof(double));
			for(j=0;j<6;j++)
			{
				if(j==0||j==1)
				{
					(zone->grid+i)->DummySlabV[j]=(double *)malloc((NK-1)*(NJ-1)*sizeof(double));
				}
				else if(j==2||j==3)
				{
					(zone->grid+i)->DummySlabV[j]=(double *)malloc((NK-1)*(NI-1)*sizeof(double));
				}
				else
				{
					(zone->grid+i)->DummySlabV[j]=(double *)malloc((NJ-1)*(NI-1)*sizeof(double));
				}
			}
		}
	}
	
	do
	{
		pressConvcri=fopen("Settings/pressureconvergCriteria.dat","r");
		fscanf(pressConvcri,"%lf",&(base->BaseParameters->PressLoopConv));
		fclose(pressConvcri);
		
		residBase=0.0;
		gdi = 0;
		
		if(base->BaseParameters->flag2d == 1)
			Multigrid2d(base,gdi);
		else
			Multigrid(base,gdi);
		
//		pressureIterCount = pressureIterCount + 1;
			
		for(itr=0;itr<3;itr++)
		{
			residBase=0.0;
			for(iZone=0;iZone<numZones;iZone++)
			{
				zone=base->CgnsZone+iZone;
				if(zone->flagSolidFluid==0)
				{
					MultigridExtractVariables(4,base,zone,gdi);
					
					switch(pressureItercase)
					{
						case 0:
								residZone=UpdatePressure2d(base,zone);
						break;
						
						case 1:
								residZone=UpdatePressure(base,zone);
						break;
						
						case 2:
								residZone=UpdatePressure2dVarDensity(base,zone);
						break;
						
						case 3:
								residZone=UpdatePressureVarDensity(base,zone);
						break;
					}
					residBase=residBase+residZone;
				}
			}
			residBase=(residBase/numZones);
//			pressureIterCount++;
		}
//		printf("residual pressure=%12.8e\n",residBase);
	}
	while(residBase > base->BaseParameters->PressLoopConv);
	
	RefreshSurfaceCells(4,base);
	
//	fprintf(vn,"%d\n",pressureIterCount);
	
	for(iZone=0;iZone<numZones;iZone++)
	{
		zone=base->CgnsZone+iZone;	
		for(i=0;i<=(zone->zonalCoarseGridIndex);i++)
		{
			free((zone->grid+i)->DummyV);
			free((zone->grid+i)->Gresidual);
			free((zone->grid+i)->E);
			free((zone->grid+i)->GRHS);
			for(j=0;j<6;j++)
			{
				free((zone->grid+i)->DummySlabV[j]);
			}
		}
		for(i=1;i<=(zone->zonalCoarseGridIndex);i++)
		{
			free((zone->grid+i)->V);
			for(j=0;j<6;j++)
			{
				free((zone->grid+i)->SlabV[j]);
			}
		}
	}
//	fclose(vn);
	return 0;
}
