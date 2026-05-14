
#include "APKCalcDiffusionCoefCommon.h"
/**
The function CalcDiffusionCoefCommon() takes the variable tag and returns its corresponding diffusion coefficient.
If the diffusion coefficient is a variable it returns the diffusion coefficient corresponding to that
cell index
*/

double CalcDiffusionCoefCommon(Base *base,Zone *zone,int variableTag,int index)
{
	double diffusionCoefficient=1.0;
	int flagviscoeff,flagdiffcoeff,flagscalcoeff,flagmagdiff,flageleccond,flagpermeability;//rajendra
	int flagdropdiffusivity; //mrunalini
	
	if(variableTag==4)//For Pressure Coeff
	{
		diffusionCoefficient = 1.0;
	}
	else if(variableTag==6)
	{
		if(base->flagTurb == 1 && zone->flagSolidFluid == 0)
		{
			if(base->flagTurbModel == 1 || base->flagTurbModel == 3)
			{
				diffusionCoefficient = *(zone->Flowsol->Tcoeff->TurbConductivity+index); 
			}
			else
			{
				diffusionCoefficient = *(zone->Flowsol->Tcoeff->EffTurbConductivity+index);
			}
		}
		else
		{
			flagdiffcoeff= zone->Flowsol->flagDiffcoeff;
			diffusionCoefficient = *(zone->Flowsol->Diffcoeff+(index*flagdiffcoeff));
		}
	}
	else if((variableTag==1)||(variableTag==2)||(variableTag==3)||(variableTag==100)||(variableTag==200)||(variableTag==300))
	{
		if(base->flagTurb == 1 && zone->flagSolidFluid == 0)
		{
			diffusionCoefficient = *(zone->Flowsol->Tcoeff->EffViscoeff+index);//Turbulent Viscosity should be called???
// 			diffusionCoefficient = *(zone->Flowsol->Tcoeff->TViscoeff+index);//Effective Viscosity should be called???
		}
		else
		{
			flagviscoeff=zone->Flowsol->flagViscoeff;
			diffusionCoefficient= *(zone->Flowsol->Viscoeff+(index*flagviscoeff));
		}
	}
	else if(variableTag==7)
	{
		flagdiffcoeff= zone->Flowsol->flagDiffcoeff;
		diffusionCoefficient = *(zone->Flowsol->Diffcoeff+(index*flagdiffcoeff));
	}
	else if((variableTag>= 600) && (variableTag<=609))
	{
		flagscalcoeff=zone->Flowsol->flagScalcoeff[variableTag-600];
		diffusionCoefficient=*(zone->Flowsol->Scalcoeff[variableTag-600]+(index*flagscalcoeff));
	}
	else if((variableTag==56)||(variableTag==57)||(variableTag==58)||(variableTag==560)||(variableTag==570)||(variableTag==580))
	{
	  flagpermeability = zone->Flowsol->flagpermeability;
          flageleccond= zone->Flowsol->flageleccond;
	
	  diffusionCoefficient= 1.0/(*(zone->Flowsol->Eleconductivity+(index*flageleccond)) * *(zone->Flowsol->Permeability+(index*flagpermeability)));//rajendra
	}
	else if((variableTag==8)||(variableTag==9))//For Electrical Parameters 
	{
		diffusionCoefficient = 1.0;
	}
	else if(variableTag==17)
	{
		flageleccond= zone->Flowsol->flageleccond;
		diffusionCoefficient = *(zone->Flowsol->Eleconductivity+(index*flageleccond));
	}
	else if(variableTag==59)//rajendra
		diffusionCoefficient = 1.0;
	else if(variableTag==68)//evaporated fuel mass-fraction
	{
		flagdropdiffusivity = zone->Flowsol->flagDropDiffcoeff;
		diffusionCoefficient = *(zone->Flowsol->DropDiffcoeff + (index*flagdropdiffusivity));
	}
	return(diffusionCoefficient);
}


