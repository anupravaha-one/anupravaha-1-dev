#include"APKTURBfreeGamma.h"
/**
This function makes the allocated memory block to the function Gamma free. For details refer to the function CalcGamma().
*/

int freeGamma(Base *base)
{
	int nz,cz;
	Zone *zone;
	
	if(base->BaseParameters->Turb->numTurbScalars >= 2)
	{
		nz=(base->Number_of_zones);
		
		if(base->flagLogKEps == 0 || base->flagTurbModel == 2)
		{
			for(cz=0;cz<nz;cz++)
			{
				zone = base->CgnsZone+cz;
				if(zone->flagSolidFluid==0)
				{
					free(zone->Flowsol->Gamma);
	
				}
			}
		}
	}

	return(100);
}

