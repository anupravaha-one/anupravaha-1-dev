#include "APKCalcVariableOneIterConjugate.h"
/**
The function CalcVariableOneIterConjugate() is exactly similar to the function CalcVariableConjugate() except for that
it does not pass the control to perform Gauss Seidel Iterations for more than once. That is it does not compute
the rms value of the residual.
The actual functionalities of this function is exactly similar to that of CalcVariableConjugate().
*/
double CalcVariableOneIterConjugate(Base *base,int variableTag)
{
	Zone *zone;
	int totalBaseCells,nZones,numZones;
	double baseResidual,zoneResidual;
	
	nZones=(base->Number_of_zones);

	if(variableTag==1||variableTag==2||variableTag==3)
		totalBaseCells = base->TotalCellFluid;
	else
		totalBaseCells = base->TotalCellDomain;
	
	baseResidual = 0.0;

	/* baseResidual is used to find RMS residual for entire geometry rather than adding RMS residual for individual zones*/
	if(base->BaseParameters->ConvectionScheme==1)//QUICK
	{
		if(base->BaseParameters->flag2d==1) //2D
		{
			for(numZones=0;numZones<nZones;numZones++)
			{
				zone = (base->CgnsZone+numZones);

				ExtractVariables(variableTag,base,zone);
		
				if(zone->flagSolidFluid==0)
				{
					zoneResidual = UpdateVarQVK2d(base,zone,variableTag);
				}
				else
				{
					zoneResidual = UpdateVarConjugate2d(base,zone,variableTag);
				}

				baseResidual = baseResidual+zoneResidual;
			}
			baseResidual = sqrt(baseResidual/totalBaseCells);
		}
		else if(base->BaseParameters->flag2d==0) //3D
		{
			for(numZones=0;numZones<nZones;numZones++)
			{
				zone = (base->CgnsZone+numZones);

				ExtractVariables(variableTag,base,zone);
		
				if(zone->flagSolidFluid==0)
				{
					zoneResidual = UpdateVarQVK(base,zone,variableTag);
				}
				else
				{
					zoneResidual = UpdateVarConjugate(base,zone,variableTag);
				}

				baseResidual = baseResidual+zoneResidual;
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

				ExtractVariables(variableTag,base,zone);

				if(zone->flagSolidFluid==0)
				{
					zoneResidual = UpdateVarCD2d(base,zone,variableTag);
				}
				else
				{
					zoneResidual = UpdateVarConjugate2d(base,zone,variableTag);
				}

				baseResidual = baseResidual+zoneResidual;
			}
			baseResidual = sqrt(baseResidual/totalBaseCells);
		}
		else if(base->BaseParameters->flag2d==0) //3D
		{
			for(numZones=0;numZones<nZones;numZones++)
			{
				zone = (base->CgnsZone+numZones);

				ExtractVariables(variableTag,base,zone);
	
				if(zone->flagSolidFluid==0)
				{
					zoneResidual = UpdateVarCD(base,zone,variableTag);
				}
				else
				{
					zoneResidual = UpdateVarConjugate(base,zone,variableTag);
				}
				
				baseResidual = baseResidual+zoneResidual;
			}
			baseResidual = sqrt(baseResidual/totalBaseCells);
		}
	}

	return (baseResidual);
}

