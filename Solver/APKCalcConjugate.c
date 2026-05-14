#include "APKCalcConjugate.h"
#include "APKsolverCoupled.h"
#include "APKCalcVariableConjugate.h"
#include "APKsolverSegregated.h"

/**
The function CalcVariableConjugate() is similar to CalcVariable(). The difference is that this function checks if the zone is a
solid zone or a fluid zone and if the zone is a solid it calls the function UpdatevarConjugate(). It calls The function
UpdateVar() if the zone is a fluid zone.
The actual functionalities of this function is similar to that of CalcVariable().
*/


int CalcConjugate(Base *base,int variableTag)
{
	Zone *zone;
	int numZones,totalcells,numBC,nZones,counter=0,iter=0;
	int Max_Potential_Flow_iter=100000, Max_Potential_Flow_TimeError=10000;
	double totalresidual,iterresidual,convcri,timeerror;

	nZones=(base->Number_of_zones);
	
	totalcells = base->TotalCellDomain;

	convcri = *(base->BaseParameters->convergenceCriterion+1);
	do{
		FillNewsolToPrevsol_T(base,variableTag);
		EvaluateEnergysource(base);
		iter=0;
		if(base->BaseParameters->flag2d==1) //2D
		{
			do
			{
				totalresidual = 0.0;
				
				// totalresidual is used to find RMS residual for entire geometry rather than adding RMS residual for individual zones
				for(numZones=0;numZones<nZones;numZones++)
				{
					zone = (base->CgnsZone+numZones);
					
					numBC = zone->Number_of_BC_Patches;

					ExtractVariables(variableTag,base,zone);

					if(!(variableTag==990 && zone->flagSolidFluid==1))
						iterresidual = UpdateVarConjugate2d(base,zone,variableTag);

					totalresidual = totalresidual+iterresidual;
				}
				totalresidual = sqrt(totalresidual/totalcells);
				totalresidual=totalresidual;
				iter++;
				if(iter%1000==0)
				printf("iter = %d \t Residue variableNew = %12.8e\n",iter,totalresidual);
				//exit(1);
				
				if(iter > Max_Potential_Flow_iter) totalresidual = 0.0;
			}
			while(totalresidual> convcri) ;
		}
		else if(base->BaseParameters->flag2d==0) //3D
		{
			do
			{
				totalresidual = 0.0;
				// totalresidual is used to find RMS residual for entire geometry rather than adding RMS residual for individual zones
				
				for(numZones=0;numZones<nZones;numZones++)
				{
					zone = (base->CgnsZone+numZones);
					numBC = zone->Number_of_BC_Patches;

					ExtractVariables(variableTag,base,zone);

					if(!(variableTag==990 && zone->flagSolidFluid==1))
						iterresidual = UpdateVarConjugate(base,zone,variableTag);

					totalresidual = totalresidual+iterresidual;
					
				}
				totalresidual = sqrt(totalresidual/totalcells);
				totalresidual=totalresidual;
				iter++;
				if(iter%1000==0)
					printf("iter = %d \t Residue variableNew = %12.8e\n",iter,totalresidual);
				if(iter > Max_Potential_Flow_iter) totalresidual = 0.0;
			}
			while(totalresidual> convcri) ;
		}
		
		TimeDepend(base->BaseParameters->Time,variableTag,base);
		RefreshSurfaceCells(variableTag,base);//supply tag for Variable
		timeerror=CalcVariableError(variableTag,base);
		
		if(counter%1000==0)
		printf("Counter = %d \t Time error = %e\n",counter,timeerror);
		if(counter > Max_Potential_Flow_TimeError) timeerror = 0.0;
		counter++;
	}while(timeerror>convcri);

	
	
	return 0;
}
