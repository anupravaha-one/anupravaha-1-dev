#include "APKCalcDiffusionCoef.h"

/**
The function CalcDiffusionCoef() takes the variable tag and returns its corresponding diffusion coefficient.
If the diffusion coefficient is a variable it returns the diffusion coefficient corresponding to that
cell index
*/
double CalcDiffusionCoef(Base *base,Zone *zone,int variableTag,int index)
{
	double diffusionCoefficient;
	int flagviscoeff,flagdiffcoeff,flagscalcoeff;
	
	if((variableTag==1)||(variableTag==2)||(variableTag==3)||(variableTag==100)||(variableTag==200)||(variableTag==300))
	{
		flagviscoeff=zone->Flowsol->flagViscoeff;
		diffusionCoefficient= *(zone->Flowsol->Viscoeff+(index*flagviscoeff));
	}
	else if(variableTag==6 || variableTag==7)
	{
		flagdiffcoeff= zone->Flowsol->flagDiffcoeff;
		diffusionCoefficient = *(zone->Flowsol->Diffcoeff+(index*flagdiffcoeff));
	}
	else if((variableTag>= 600) && (variableTag<=609))
	{
		flagscalcoeff=zone->Flowsol->flagScalcoeff[variableTag-600];
		diffusionCoefficient=*(zone->Flowsol->Scalcoeff[variableTag-600]+(index*flagscalcoeff));
	}
	return(diffusionCoefficient);
}


