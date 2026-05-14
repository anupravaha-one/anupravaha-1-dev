#include "APKGetSurfacePoints.h"

/**
This function is used to find the x,y and z coordinates of all points making a cell face.
 
The grid-coordinate pointers are supplied to 
the function for vertex coordinates, the face and cell index is identified by the contents of the array temp[][]. 
The values taken by the array temp[2][3] are listed below for different faces-

I_min   
temp[0][0]=NI*NJ,temp[0][1]=NI,temp[0][2]=1,temp[1][0]=k,temp[1][1]=j,temp[1][2]=i;

I_max   
temp[0][0]=NI*NJ,temp[0][1]=NI,temp[0][2]=1,temp[1][0]=k,temp[1][1]=j,temp[1][2]=i+1;  

J_min   
temp[0][0]=NI*NJ,temp[0][1]=1,temp[0][2]=NI,temp[1][0]=k,temp[1][1]=i,temp[1][2]=j;

J_max   
temp[0][0]=NI*NJ,temp[0][1]=1,temp[0][2]=NI,temp[1][0]=k,temp[1][1]=i,temp[1][2]=j+1;

K_min   
temp[0][0]=NI,temp[0][1]=1,temp[0][2]=NI*NJ,temp[1][0]=j,temp[1][1]=i,temp[1][2]=k;

K_max   
temp[0][0]=NI,temp[0][1]=1,temp[0][2]=NI*NJ,temp[1][0]=j,temp[1][1]=i,temp[1][2]=k+1;

The x, y and z coordinates of all the four points of a face are returned in SurfPoints[4][3].
*/
int GetSurfacePoints(double *coordX,double *coordY,double *coordZ,double surfacePoints[4][3],int temp[2][3])
{
	 surfacePoints[0][0]=*(coordX+temp[0][0]* temp[1][0] +  temp[0][1]*temp[1][1]  +  temp[0][2]*temp[1][2]);
	 surfacePoints[0][1]=*(coordY+temp[0][0]* temp[1][0] +  temp[0][1]*temp[1][1]  +  temp[0][2]*temp[1][2]);
	 surfacePoints[0][2]=*(coordZ+temp[0][0]* temp[1][0] +  temp[0][1]*temp[1][1]  +  temp[0][2]*temp[1][2]);
   
	 surfacePoints[1][0]=*(coordX+temp[0][0]* temp[1][0] +  temp[0][1]*(temp[1][1]+1)+temp[0][2]*temp[1][2]);
	 surfacePoints[1][1]=*(coordY+temp[0][0]* temp[1][0] +  temp[0][1]*(temp[1][1]+1)+temp[0][2]*temp[1][2]);
	 surfacePoints[1][2]=*(coordZ+temp[0][0]* temp[1][0] +  temp[0][1]*(temp[1][1]+1)+temp[0][2]*temp[1][2]);

	 surfacePoints[2][0]=*(coordX+temp[0][0]*(temp[1][0]+1)+temp[0][1]*(temp[1][1]+1)+temp[0][2]*temp[1][2]);
	 surfacePoints[2][1]=*(coordY+temp[0][0]*(temp[1][0]+1)+temp[0][1]*(temp[1][1]+1)+temp[0][2]*temp[1][2]);
	 surfacePoints[2][2]=*(coordZ+temp[0][0]*(temp[1][0]+1)+temp[0][1]*(temp[1][1]+1)+temp[0][2]*temp[1][2]);

	 surfacePoints[3][0]=*(coordX+temp[0][0]*(temp[1][0]+1)+temp[0][1]* temp[1][1]  + temp[0][2]*temp[1][2]);
	 surfacePoints[3][1]=*(coordY+temp[0][0]*(temp[1][0]+1)+temp[0][1]* temp[1][1]  + temp[0][2]*temp[1][2]);
	 surfacePoints[3][2]=*(coordZ+temp[0][0]*(temp[1][0]+1)+temp[0][1]* temp[1][1]  + temp[0][2]*temp[1][2]);
	 return 0;
}
