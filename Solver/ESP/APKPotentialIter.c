#include "APKPotentialIter.h"

int PotentialIter(Base *base, int variableTag)
{

	FILE *PotConv,*Tstep;
	Zone *zone;
	int PotItercase,stp;
	double residBase,residZone,PotConvCrit;
	int totalBaseCells=0,iZone,numZones;
	numZones=base->Number_of_zones;
	totalBaseCells = base->TotalCellFluid;
	
	if(base->BaseParameters->flag2d == 1 && base->flagVariableDensityProblem == 0)
		PotItercase = 0; /*case for calling UpdatePotential2d*/
	
	if(base->BaseParameters->flag2d == 0 && base->flagVariableDensityProblem == 0)
		PotItercase = 1; /*case for calling UpdatePotential, this is for normal 3d*/
	
	if(base->BaseParameters->flag2d == 1 && base->flagVariableDensityProblem == 1)
		PotItercase = 2; /*case for calling UpdatePotential2dVarDensity*/
	
	if(base->BaseParameters->flag2d == 0 && base->flagVariableDensityProblem == 1)
		PotItercase = 3; /*case for calling UpdatePotentialVarDensity, this is for 3d variable Density case*/
	
	//printf("PressureItercase = %d\n",PotItercase);
	
	do
	{
		residBase=0.0;
		PotConv=fopen("Settings/ESPInput/PotentialConvergence_TDMA.dat","r");
		fscanf(PotConv,"%lf",&PotConvCrit);
		fclose(PotConv);
		Tstep=fopen("Settings/ESPInput/EfieldTimeSteps.dat","r");
		fscanf(Tstep,"%d",&stp);
		fclose(Tstep);

		for(iZone=0;iZone<numZones;iZone++)
		{
			zone=base->CgnsZone+iZone;
			if(zone->flagSolidFluid==0)
			{
				if(zone->Number_of_BC_Patches>=3)
				{
					//RefreshSurfaceCells(4,base);
					ExtractVariables(8,base,zone);//4 - pressure
					
					switch(PotItercase)
					{
						case 0:
							residZone=UpdatePotential2d(base,zone);
							//residZone=UpdatePotentialArtificial2d(base,zone);
						break;
						/*
						case 1:
							residZone=UpdatePotential(base,zone);
						break;
						
						case 2:
							residZone=UpdatePotential2dVarDensity(base,zone);
						break;
						
						case 3:
							residZone=UpdatePotentialVarDensity(base,zone);
						break;
						*/
					}
					residBase=residBase+residZone;
				}
			}
		}
		
		for(iZone=0;iZone<numZones;iZone++)
		{
			zone=base->CgnsZone+iZone;
			if(zone->flagSolidFluid==0)
			{
				if(zone->Number_of_BC_Patches<3)
				{
					ExtractVariables(8,base,zone);
					
					switch(PotItercase)
					{
						case 0:
							residZone=UpdatePotential2d(base,zone);
							//residZone=UpdatePotentialArtificial2d(base,zone);
						break;
						/*
						case 1:
							residZone=UpdatePotential(base,zone);
						break;
						
						case 2:
							residZone=UpdatePotential2dVarDensity(base,zone);
						break;
						
						case 3:
							residZone=UpdatePotentialVarDensity(base,zone);
						break;
						*/
					}
					
					residBase=residBase+residZone;
				}
			}
		}
		
		residBase=(residBase/numZones);
		if((stp==1)||(stp==2))
			printf("\nresidual Potential=%12.8e",residBase);
		
	}
	while(residBase > PotConvCrit);
	RefreshSurfaceCells(8,base);
		
	return 0;
}
