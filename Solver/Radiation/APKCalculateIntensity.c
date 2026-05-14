#include "APKCalculateIntensity.h"

/** This function is used to call the two functions Solution2DConstRefIndex() and  Solution2DVariRefIndex() to solve RTE in constant refractive
index medium and variable refractive index medium respectively*/

int CalculateIntensity(Base *base, Zone *zone, int variableTag)
	{
	 double *variableDomain, *preVariableDomain, alpha;
	 
	 
	 variableDomain=MakeDomainArray(zone, variableTag);
	 preVariableDomain=MakeDomainArray(zone,variableTag*10);
	 
	 alpha=base->BaseParameters->RatioCenUpRadiation;
		 
	 if(base->BaseParameters->flag2d==1)
	 	{
		 if(zone->Flowsol->flagRefIndex==0)
		  Solution2DConstRefIndex(base, zone, variableDomain, preVariableDomain, alpha);
		 
		 else
		  {
		  //EvaluateVarCoefficients(base);
		 // Solution2DVariRefIndex(base, zone, variableDomain, preVariableDomain, alpha);
		  }
		}
	 if(base->BaseParameters->flag2d==0)
	 	{
		 if(zone->Flowsol->flagRefIndex==0)
		  Solution3DConstRefIndex(base, zone, variableDomain, preVariableDomain, alpha);
		 
		 else
		 {
		//  EvaluateVarCoefficients(base);
		 // Solution3DVariRefIndex(base, zone, variableDomain, prevVariableDomain, alpha);
		  }
		}
	 return 0;
	}
