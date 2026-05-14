#include "APKCheckOutward.h"

/**
This function is used to check whether the calculated surface area vector is an outward normal to that surface or not.

The surface area vectors of the cell faces are used assuming them to be outward. To ensure this, we take the dot 
product of the surface vector and the vector extending from the cell center to the face center. If the dot 
product comes out to be negative, the surface area vector is multiplied by (-1). 
*/
int CheckOutward(double surface[3],double surfPoint[4][3],double centre[3])
{
	int comp,point;
	double faceCentre[3],testVector[3],dotProduct;
	for(comp=0;comp<3;comp++)
	{
		faceCentre[comp]=0.0;
		for(point=0;point<4;point++)
		faceCentre[comp]=faceCentre[comp]+surfPoint[point][comp];
		faceCentre[comp]=faceCentre[comp]/4.0;
		      
		testVector[comp]=faceCentre[comp]-centre[comp];
	}
	dotProduct=0.0;
	for(comp=0;comp<3;comp++)
		dotProduct=dotProduct+testVector[comp]*surface[comp];
	/*Checking if dotproduct is negative*/ 
	if(dotProduct<0.0)
	for(comp=0;comp<3;comp++)
		surface[comp]=(-1.0)*surface[comp];
	return 0;
}
