#include "APKCalcVariableConjugate.h"
/**
The function CalcVariableConjugate() is similar to CalcVariable(). The difference is that this function checks if the zone is a
solid zone or a fluid zone and if the zone is a solid it calls the function UpdatevarConjugate(). It calls The function
UpdateVar() if the zone is a fluid zone.
The actual functionalities of this function is similar to that of CalcVariable().
*/
int CalcVariableConjugate(Base *base,int variableTag)
{
	Zone *zone;
	int numZones,totalcells,nZones,iteru=0;
	double totalresidual,iterresidual,convcri;
	nZones=(base->Number_of_zones);
	
	if(variableTag==1||variableTag==2||variableTag==3 ||base->flagMPFV==1 || (base->flagTurb==1 && (variableTag==600||variableTag==601)))
		totalcells = base->TotalCellFluid;
	else
		totalcells = base->TotalCellDomain;
	
	convcri =find_convcri(variableTag,base);

	if(base->BaseParameters->ConvectionScheme==1)//QUICK
	{
		if(base->BaseParameters->flag2d==1) //2D
		{
			do
			{
				iteru = iteru+1;
				totalresidual = 0.0;

				// totalresidual is used to find RMS residual for entire geometry rather than adding RMS residual for individual zones
				for(numZones=0;numZones<nZones;numZones++)
				{
					zone = (base->CgnsZone+numZones);
					
					ExtractVariables(variableTag,base,zone);

					if(zone->flagSolidFluid==0)
					{
						iterresidual = UpdateVarQVK2d(base,zone,variableTag);
					}
					else
					{
						iterresidual = UpdateVarConjugate2d(base,zone,variableTag);
					}
					
					totalresidual = totalresidual+iterresidual;
				}
				totalresidual = sqrt(totalresidual/totalcells);
				
				if(iteru%50 == 0)
					printf("GS IterCount=%d\tResidual of Variable[%d] = %e\n", iteru,variableTag,totalresidual);
			}
			while(totalresidual > convcri);
		}
		else if(base->BaseParameters->flag2d==0) //3D
		{
			do
			{
				iteru = iteru+1;
				totalresidual = 0.0;

				// totalresidual is used to find RMS residual for entire geometry rather than adding RMS residual for individual zones
				for(numZones=0;numZones<nZones;numZones++)
				{
					zone = (base->CgnsZone+numZones);
					
					ExtractVariables(variableTag,base,zone);
					
					if(zone->flagSolidFluid==0)
					{
						iterresidual = UpdateVarQVK(base,zone,variableTag);
					}
					else
					{
						iterresidual = UpdateVarConjugate(base,zone,variableTag);
					}
					
					totalresidual = totalresidual+iterresidual;
				}
				totalresidual = sqrt(totalresidual/totalcells);
				
				if(iteru%50 == 0)
					printf("GS IterCount=%d\tResidual of Variable[%d] = %e\n", iteru,variableTag,totalresidual);
			}
			while(totalresidual > convcri);
		}
	}
	else if(base->BaseParameters->ConvectionScheme==0)//CD
	{
		if(base->BaseParameters->flag2d==1) //2D
		{
			do
			{
				iteru = iteru+1;
				totalresidual = 0.0;

				//totalresidual is used to find RMS residual for entire geometry rather than adding RMS residual for individual zones
				for(numZones=0;numZones<nZones;numZones++)
				{
					zone = (base->CgnsZone+numZones);
					
					ExtractVariables(variableTag,base,zone);
					
					if(zone->flagSolidFluid==0)
					{
						iterresidual = UpdateVarCD2d(base,zone,variableTag);
					}
					else
					{
						iterresidual = UpdateVarConjugate2d(base,zone,variableTag);
					}
					
					//printf("%s Zone Residual of Variable[%d] = %e\n", zone->zonename,variableTag,iterresidual);
					totalresidual = totalresidual+iterresidual;
				}
				totalresidual = sqrt(totalresidual/totalcells);
				
				if(iteru%50 == 0)
					printf("GS IterCount=%d\tResidual of Variable[%d] = %e\n", iteru,variableTag,totalresidual);
			}
			while(totalresidual > convcri);
		}
		else if(base->BaseParameters->flag2d==0) //3D
		{
			do
			{
				iteru = iteru+1;
				totalresidual = 0.0;

				// totalresidual is used to find RMS residual for entire geometry rather than adding RMS residual for individual zones
				for(numZones=0;numZones<nZones;numZones++)
				{
					zone = (base->CgnsZone+numZones);
					
					ExtractVariables(variableTag,base,zone);
					
					if(zone->flagSolidFluid==0)
					{
						iterresidual = UpdateVarCD(base,zone,variableTag);
					}
					else
					{
						iterresidual = UpdateVarConjugate(base,zone,variableTag);
					}
					
					totalresidual = totalresidual+iterresidual;
				}
				totalresidual = sqrt(totalresidual/totalcells);
				
				if(iteru%50 == 0)
					printf("GS IterCount=%d\tResidual of Variable[%d] = %e\n", iteru,variableTag,totalresidual);
			}
			while(totalresidual > convcri);
		}
	}

	RefreshSurfaceCells(variableTag,base);
	
	return 0;
}
