#include "APKface_center.h"

/**
This function is used to find a specific coordinate (x, y, or z) of the face center at the boundary.

As input arguments it takes the grid coordinate pointer of the vertex with minimum indices of the face, corresponding to 
the coordinate of the center to be found, direction of the face and the vertex size of the Zone.
*/
double face_center(double *crd,int dir,int ni,int nj, int nk)
{
	double cent_crd;

	switch(dir)
	{
		case 1:
		cent_crd =(*(crd)+*(crd+ni)+*(crd+ni*nj)+*(crd+ni*nj+ni))*0.25;
		break;

		case 2:
		cent_crd =(*(crd)+*(crd+ni)+*(crd+ni*nj)+*(crd+ni*nj+ni))/4;
		break;

		case 3:
		cent_crd =(*(crd)+*(crd+1)+*(crd+ni*nj)+*(crd+ni*nj+1))/4;
		break;

		case 4:
		cent_crd =(*(crd)+*(crd+1)+*(crd+ni*nj)+*(crd+ni*nj+1))/4;
		break;

		case 5:
		cent_crd =(*(crd)+*(crd+1)+*(crd+ni)+*(crd+ni+1))/4;
		break;

		case 6:
		cent_crd =(*(crd)+*(crd+1)+*(crd+ni)+*(crd+ni+1))/4;
		break;
	}
	return (cent_crd);
}
