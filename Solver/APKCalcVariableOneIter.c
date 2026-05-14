#include "APKCalcVariableOneIter.h"
/**
The function CalcVariableOneIter() is exactly similar to the function CalcVariable() except for that
it does not pass the control to perform Gauss Seidel Iterations for more than once. That is it does not compute
the rms value of the residual.
The actual functionalities of this function is exactly similar to that of CalcVariable().
*/
double CalcVariableOneIter(Base *base,int variableTag)
{
	Zone *zone;
	int totalBaseCells,numZones,nZones,SolveSolidTag=1;
	double baseResidual,zoneResidual;
	
	nZones=(base->Number_of_zones);

	if(variableTag==1||variableTag==2||variableTag==3 || (base->flagTurb==1 && (variableTag==600||variableTag==601)))
	{
		SolveSolidTag = 0;
		totalBaseCells = base->TotalCellFluid;
	}
	else
	{
		SolveSolidTag = 1;
		totalBaseCells = base->TotalCellDomain;
	}
	
	baseResidual = 0.0;

	if(base->BaseParameters->ConvectionScheme==1)//QUICK
	{
		if(base->BaseParameters->flag2d==1) //2D
		{
			for(numZones=0;numZones<nZones;numZones++)
			{
				zone = (base->CgnsZone+numZones);

				if(zone->flagSolidFluid==0 || SolveSolidTag == 1)
				{
					ExtractVariables(variableTag,base,zone);

					zoneResidual = UpdateVarQVK2d(base,zone,variableTag);

					baseResidual = baseResidual+zoneResidual;
				}
			}

			baseResidual = sqrt(baseResidual/totalBaseCells);
		}
		else if(base->BaseParameters->flag2d==0) //3D
		{
			for(numZones=0;numZones<nZones;numZones++)
			{
				zone = (base->CgnsZone+numZones);

				if(zone->flagSolidFluid==0 || SolveSolidTag == 1)
				{
					ExtractVariables(variableTag,base,zone);

					zoneResidual = UpdateVarQVK(base,zone,variableTag);

					baseResidual = baseResidual+zoneResidual;
				}
			}

			baseResidual = sqrt(baseResidual/totalBaseCells);
		}
	}
	else if(base->BaseParameters->ConvectionScheme==0)//CD
	{
		if(base->BaseParameters->flag2d==1) //2D
		{
			for(numZones=0;numZones<nZones;numZones++)
			{
				zone = (base->CgnsZone+numZones);

				if(zone->flagSolidFluid==0 || SolveSolidTag == 1)
				{
					ExtractVariables(variableTag,base,zone);

					zoneResidual = UpdateVarCD2d(base,zone,variableTag);

					baseResidual = baseResidual+zoneResidual;
				}
			}

			baseResidual = sqrt(baseResidual/totalBaseCells);
		}
		else if(base->BaseParameters->flag2d==0) //3D
		{
			for(numZones=0;numZones<nZones;numZones++)
			{
				zone = (base->CgnsZone+numZones);

				if(zone->flagSolidFluid==0 || SolveSolidTag == 1)
				{
					ExtractVariables(variableTag,base,zone);
					
					zoneResidual = UpdateVarCD(base,zone,variableTag);
					
					baseResidual = baseResidual+zoneResidual;
				}
			}

			baseResidual = sqrt(baseResidual/totalBaseCells);
		}
	}
	//printf("base residual = %e\n",baseResidual);
	/* baseResidual is used to find RMS residual for entire geometry rather than adding RMS residual for individual zones*/
	RefreshSurfaceCells(variableTag,base);/*supply tag for Variable*/
	return (baseResidual);
}

double CalcVariableOneIter_porousMedia(Base *base,int variableTag)
{
	int totalBaseCells;
	double baseResidual,zoneResidual;
	int *ptrSize,nZones;
	Zone *zone;
	int numZones;
	double *tempVar,*prevTstepVar,*variableDomain;
	//printf("i am in CalcVariableOneIter\n");
	nZones=(base->Number_of_zones);

	if(variableTag==1||variableTag==2||variableTag==3)
		totalBaseCells = base->TotalCellFluid;
	else
		totalBaseCells = base->TotalCellDomain;
	
	baseResidual = 0.0;
	
	for(numZones=0;numZones<nZones;numZones++)
	{
		zone = (base->CgnsZone+numZones);
				
			ptrSize= zone->Size;
			variableDomain=MakeDomainArray(zone,variableTag); //variableDomain stores variable pointer
			prevTstepVar=MakeDomainArray(zone,variableTag*10);   //prevTstepVar stores prevFlowsol->variable 
			tempVar=(double*)malloc((*(ptrSize+3)) * (*(ptrSize+4)) * (*(ptrSize+5))*sizeof(double));
		
			ExtractVariables(variableTag,base,zone);

			fillnewtoprev(variableDomain,tempVar,ptrSize);/*tempVar stores variable of previous iteration (not time step)*/

			zoneResidual= UpdateScaler(base,zone,variableDomain,prevTstepVar,variableTag,tempVar);

			baseResidual = baseResidual+zoneResidual;

			free(tempVar);
		
	}
	baseResidual = sqrt(baseResidual/totalBaseCells);
	/* baseResidual is used to find RMS residual for entire geometry rather than adding RMS residual for individual zones*/
	RefreshSurfaceCells(variableTag,base);/*supply tag for Variable*/
	return (baseResidual);
}


