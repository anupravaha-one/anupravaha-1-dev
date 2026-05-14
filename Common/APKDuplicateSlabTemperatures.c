#include "APKDuplicateSlabTemperatures.h"

int DuplicateSlabTemperatures(Zone *zone,double *slabTPtrArray[6],double *dummySlabTPtrArray[6])
{
	int dir;
	for(dir=0; dir<6; dir++)
		slabTPtrArray[dir] = dummySlabTPtrArray[dir];
	return 0;
}

