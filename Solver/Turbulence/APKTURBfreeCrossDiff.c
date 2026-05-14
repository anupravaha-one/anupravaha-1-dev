#include"APKTURBfreeCrossDiff.h"
/**
This function makes the allocated memory block to the variable 'CrossDiff' free. For details refer to the function CrossDiffusion().
*/

int freeCrossDiff(Base *base)
{
	int nz,cz;
	Zone *zone;

	nz=(base->Number_of_zones);

	for(cz=0;cz<nz;cz++)
	{
		zone = base->CgnsZone+cz;
		if(zone->flagSolidFluid==0)
		{
			free(zone->Flowsol->CrossDiff);
		}
	}

	return(100);
}

