#include "APKCalSurf.h"
/**
This function is used to find the surface area vector of a surface passing through four coplanner points.
 
As input arguments it takes a pointer to an array containing the x, y, z-coordinates of the four points. The surface 
area are calculated as half of the cross product of two cross diagonals passing throught the vertices. 
*/
int CalSurf(double areacof,double *points,double *vector)
{
	int comp;
	double diagonal13[3],diagonal42[3];
	for(comp=0;comp<3;comp++)
	{
		diagonal13[comp]=*(points+6+comp)-*(points+comp);
		diagonal42[comp]=*(points+3+comp)-*(points+9+comp);
	}
	*(vector+0)=0.5*(diagonal13[1]*diagonal42[2]-diagonal13[2]*diagonal42[1]);
	*(vector+1)=0.5*(diagonal13[2]*diagonal42[0]-diagonal13[0]*diagonal42[2]);
	*(vector+2)=0.5*(diagonal13[0]*diagonal42[1]-diagonal13[1]*diagonal42[0]);
	
	return 0;
}

