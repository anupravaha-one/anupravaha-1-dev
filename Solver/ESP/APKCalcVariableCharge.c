#include "APKCalcVariableCharge.h"

/**
This function CalcVariableCharge() is similar to the CalcVariable() function,
and only solves for Charge Density with a pseudo-transient term.
The variable in base data structure, base->ChargeConv is set to "1" initially
in SolverElectric.c which implies that charge is not converged or error between
calculated and Peek's electric field is morethan the given tolerance, so this
function in turn calls RefreshWireCharge() function for updating the value of
charge density on wire.

*/
int CalcVariableCharge(Base *base,int variableTag)
{
	int iter,iteru,toggle=0;
	int totalcells,numBC;
	double totalresidual,iterresidual,ErrorZone = 1.0e-10;
	int nZones;
	Zone *zone;
	int numZones,stp=0;

	nZones=(base->Number_of_zones);
	FILE *Tstep, *ChargeConv;
	iteru = 0;

#ifdef FILE_OPEN
{

	Tstep=fopen("Settings/ESPInput/EfieldTimeSteps.dat","r");
	fscanf(Tstep,"%d",&stp);
	fclose(Tstep);

	totalcells = base->TotalCellDomain;

	ChargeConv=fopen("Settings/ESPInput/ChargeConvergence_GS.dat","r");
	fscanf(ChargeConv,"%lf",&ErrorZone);
	fclose(ChargeConv);
}
#endif
	if(base->BaseParameters->ConvectionScheme==0)//CD
	{
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
					for(iter=0;iter<5;iter++)
					{
						ExtractVariables(variableTag,base,zone);

						if(base->BaseParameters->flag2d==1) //2D
						iterresidual = UpdateVarCD2d(base,zone,variableTag);	
						else
						iterresidual = UpdateVarCD(base,zone,variableTag);	
						//printf("Returned value =%e\n",iterresidual);
					}
					totalresidual = totalresidual+iterresidual;
				}
			}
			
// 			for(numZones=0;numZones<nZones;numZones++)
// 			{
// 				zone = (base->CgnsZone+numZones);
// 				
// 				if(zone->flagSolidFluid==0)
// 				{
// 					numBC = zone->Number_of_BC_Patches;
// 					if(numBC<3)
// 					{
// 						for(iter=0;iter<5;iter++)
// 						{
// 							ExtractVariables(variableTag,base,zone);
// 	
// 							if(base->BaseParameters->flag2d==1) //2D
// 							iterresidual = UpdateVarCD2d(base,zone,variableTag);
// 							else
// 							iterresidual = UpdateVarCD(base,zone,variableTag);
// 						}
// 						totalresidual = totalresidual+iterresidual;
// 					}
// 				}
// 			}
			totalresidual = sqrt(totalresidual/totalcells);
			if((stp==1)||(stp==2))
				printf("\nCharge Density Total Residue = %e\n", totalresidual);
				//exit(12);
		}
		while(fabs(totalresidual)> 1e-08);
	}
	else
	{
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
					for(iter=0;iter<5;iter++)
					{
						ExtractVariables(variableTag,base,zone);

						if(base->BaseParameters->flag2d==1) //2D
						iterresidual = UpdateVarQVK2d(base,zone,variableTag);	
						else
						iterresidual = UpdateVarQVK(base,zone,variableTag);	
					}
					totalresidual = totalresidual+iterresidual;
				}
			}
			
// 			for(numZones=0;numZones<nZones;numZones++)
// 			{
// 				zone = (base->CgnsZone+numZones);
// 				
// 				if(zone->flagSolidFluid==0)
// 				{
// 					numBC = zone->Number_of_BC_Patches;
// 					
// 					if(numBC<3)
// 					{
// 						for(iter=0;iter<5;iter++)
// 						{
// 							ExtractVariables(variableTag,base,zone);
// 	
// 							if(base->BaseParameters->flag2d==1) //2D
// 							iterresidual = UpdateVarQVK2d(base,zone,variableTag);
// 							else
// 							iterresidual = UpdateVarQVK(base,zone,variableTag);
// 						}
// 						totalresidual = totalresidual+iterresidual;
// 					}
// 				}
// 			}
			totalresidual = sqrt(totalresidual/totalcells);
			if((stp==1)||(stp==2))
				printf("\nCharge Density Total Residue = %e\n", totalresidual);
		}
		while(fabs(totalresidual)> 1e-08);
	}

	RefreshSurfaceCells(variableTag,base);/*supply tag for Variable*/

	if(base->ChargeConv==1) //Implies charge is not converged
		RefreshWireCharge(variableTag,base);
	
	//Now implement the wire charge boundary condition
	
	return 0;
}
