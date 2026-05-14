#include "APKCalcVariablePotential.h"
/**
This function CalcVariablePotential() is similar to the CalcVariable() function,
and only solves for electric potential with a pseudo-transient term.

*/
int CalcVariablePotential(Base *base,int variableTag)
{
	int iteru,toggle=0,zoneCell;
	int totalcells,numBC;
	double totalresidual,iterresidual,ResidZone,ErrorZone = 1.0e-08,iter=0;
	int nZones;
	Zone *zone;
	int numZones,stp=0;
	
	nZones=(base->Number_of_zones);
	FILE *Tstep, *PotConv;
	iteru = 0;
	
	totalcells = base->TotalCellDomain;
	
#ifdef FILE_OPEN
{

	Tstep=fopen("Settings/ESPInput/EfieldTimeSteps.dat","r");
	fscanf(Tstep,"%d",&stp);
	fclose(Tstep);

	PotConv=fopen("Settings/ESPInput/PotentialConvergence_TDMA.dat","r");
	fscanf(PotConv,"%lf",&ErrorZone);
	fclose(PotConv);
}
#endif

	do/* Iteration loop starts here*/
	{		
		toggle=toggle+1;
		iteru = iteru+1;
		//totalcells = 0;
		totalresidual = 0.0;

		/* totalresidual is used to find RMS residual for entire geometry rather than adding RMS residual for individual zones*/
		for(numZones=0;numZones<nZones;numZones++)
		{
			zone = (base->CgnsZone+numZones);
	
			numBC = zone->Number_of_BC_Patches;
			//if(numBC>=3)
			{
				//do
				{
					ExtractVariables(variableTag,base,zone);//Misc  

					zoneCell=CountCells(zone);
					if(base->BaseParameters->flag2d==1) //2D
					iterresidual = UpdateVarCD2d(base,zone,variableTag);	
					else
					iterresidual = UpdateVarCD(base,zone,variableTag);	
					ResidZone = sqrt(iterresidual/zoneCell);
					iter++;
				}
				//while(ResidZone> ErrorZone);
				totalresidual = totalresidual+iterresidual;
			}
		}
		
// 		for(numZones=0;numZones<nZones;numZones++)
// 		{
// 			zone = (base->CgnsZone+numZones);
// 			numBC = zone->Number_of_BC_Patches;
// 				
// 			if(numBC<3)
// 			{
// 				do
// 				{	
// 					ExtractVariables(variableTag,base,zone);//Misc  
// 
// 					zoneCell=CountCells(zone);
// 		
// 					if(base->BaseParameters->flag2d==1) //2D
// 					iterresidual = UpdateVarCD2d(base,zone,variableTag);	
// 					else
// 					iterresidual = UpdateVarCD(base,zone,variableTag);	
// 
// 					ResidZone = sqrt(iterresidual/zoneCell);
// 				}
// 				while(ResidZone> ErrorZone);
// 				totalresidual = totalresidual+iterresidual;
// 			}
// 		}
		totalresidual = sqrt(totalresidual/totalcells);
		if((stp==1)||(stp==2))
			printf("\nPotential Total Residue = %e\n", totalresidual);
		
	}
	while(totalresidual> ErrorZone);
	
	RefreshSurfaceCells(variableTag,base);/*supply tag for Variable*/
	return 0;
}

