#include "APKDiffusionCoeffCalculationForVariables.h"



int DiffusionCoeffCalculationForVariables(Base *base)
{
	int i,iScal;
	
	i = 0;
	if(base->flaginduced==1)//rajendra
	{
		CallCalcCoefficientsForVariables(base,i,56);
	}
	
	CallCalcCoefficientsForVariables(base,i,1);

	if(base->CgnsZone->Flowsol->T!=NULL	&&	base->flagTurbComb==0)
	{
	CallCalcCoefficientsForVariables(base,i,6);
	}

	if(base->CgnsZone->Flowsol->enthalpy!=NULL)
	{
	CallCalcCoefficientsForVariables(base,i,7);
	}

	if(base->CgnsZone->Flowsol->MixFrac!=NULL)
	{
		CallCalcCoefficientsForVariables(base,i,14);
	}
	

	for(iScal=0;iScal<10;iScal++)
	{
		if(base->CgnsZone->Flowsol->Scal[iScal]!=NULL)
		{
			CallCalcCoefficientsForVariables(base,i,600+iScal);
		}
	}

	return 0;
}
