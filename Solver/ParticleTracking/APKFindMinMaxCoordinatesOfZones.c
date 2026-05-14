#include "APKFindMinMaxCoordinatesOfZones.h"
 
int FindMinMaxCoordinatesOfZones(Base *base)
{
	int nZones,izone,NI,NJ,NK,i,j,k,pointIndex;
	Zone *zone;
	double *coordX,*coordY,*coordZ,min[3],max[3];
	
	nZones=((base)->Number_of_zones);
	for(izone=0;izone<nZones;izone++)
	{
		zone = base->CgnsZone+izone;
		
		coordX = (zone->CgnsGrid)->X;
		coordY = (zone->CgnsGrid)->Y; 
		coordZ = (zone->CgnsGrid)->Z;
		
		max[0] = min[0] = *coordX;
		max[1] = min[1] = *coordY;
		max[2] = min[2] = *coordZ; 
		
		NI=*((zone)->Size+0);
		NJ=*((zone)->Size+1);
		NK=*((zone)->Size+2);
		
		for(k=0;k<NK;k++)
		for(j=0;j<NJ;j++)
		for(i=0;i<NI;i++)
		{
			pointIndex = NI*NJ*k+NI*j+i;
			if(*(coordX+pointIndex) < min[0])
				min[0] = *(coordX+pointIndex);
			else if(*(coordX+pointIndex) > max[0])
				max[0] = *(coordX+pointIndex);
			
			if(*(coordY+pointIndex) < min[1])
				min[1] = *(coordY+pointIndex);
			else if(*(coordY+pointIndex) > max[1])
				max[1] = *(coordY+pointIndex);
			
			if(*(coordZ+pointIndex) < min[2])
				min[2] = *(coordZ+pointIndex);
			else if(*(coordZ+pointIndex) > max[2])
				max[2] = *(coordZ+pointIndex);
		}
		for(i=0;i<3;i++)
		{
			(zone->CgnsGrid)->Min[i] = min[i];
			(zone->CgnsGrid)->Max[i] = max[i];
		}
		//printf("ZoneNumber=%d\tMinX = %12.8lf\tMinY = %12.8lf\tMinZ = %12.8lf\n",izone,min[0],min[1],min[2]);
		//printf("ZoneNumber=%d\tMaxX = %12.8lf\tMaxY = %12.8lf\tMaxZ = %12.8lf\n",izone,max[0],max[1],max[2]);	
	}
	return 0;
}


