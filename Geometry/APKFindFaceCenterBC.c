#include "APKFindFaceCenterBC.h"

/**
This function is used to find the coordinates of the center of a cell face. 

The inputs are the coordinates of the four vertices of the cell face. The face center
coordinates are found by averaging the cooresponding coordinates of the face vertices.
*/
int FindFaceCenterBC(double fourPoints[][3],double center[3])
{
	int i,j;
	double temp;
	for(j=0;j<3;j++)
	{
		temp=0.0;
		for(i=0;i<4;i++)
			temp=temp+fourPoints[i][j];
		center[j]=temp/4.0;
	}
	return 0;
}
