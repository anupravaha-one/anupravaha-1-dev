#include "APKFindBeta.h"

/**
This function is used to find the diffusion lengths called Beta.
 
 he procedure for calculating Beta for a cell 
face is as follows-Calculation of the vector distance between - two cell centers across a cell face and the 
vector distances between the centers of two opposite edges using the coordinates of the face vertices. These 
distances are denoted by dX1[3], dX2[3] and dX3[3] respectively. The direction of these vectors in the direction 
of the increasing cell indices(i,j,k). The magnitude of all the three vectors is found by taking the dot product 
of a vector with itself and taking the square root of the dot product. The unit vectors in the direction of 
these vectors are found by dividing all the components of these vectors by their respective magnitudes.
The surface area vectors of the cell faces are to be represented as the linear combination of these three unit 
vectors, which needs to solve a system of three linear equations to find the coefficients of the linear 
combination denoted by alpha[3]. The system of three linear equations is solved using Cramers' rule which gives 
all the three components of alpha[3]. Beta[3] are defined as alpha[i] divided the magnitude of the ith vector. 
 */
int FindBeta(int face,double centreCell[3],double centreNeighbour[3],double facePoints[][3],double surf[3],double beta[3])
{
	int comp;
	double dX1[3],dX2[3],dX3[3],magdX1,magdX2,magdX3,n1[3],n2[3],n3[3],Det1,Det2,Det3,Det,alpha[3];
	/*finding the vectors */
	switch(face)
	{
		case 0:// West Face
			for(comp=0;comp<3;comp++)
			{
				dX1[comp]=(centreCell[comp]-centreNeighbour[comp]);
				dX2[comp]=((facePoints[1][comp]+facePoints[2][comp])/2.0-(facePoints[0][comp]+facePoints[3][comp])/2.0);
				dX3[comp]=(facePoints[2][comp]+facePoints[3][comp])/2.0-(facePoints[0][comp]+facePoints[1][comp])/2.0;
				/* capital X represents Vector*/
			}
		break;
		 case 1:// East Face
			for(comp=0;comp<3;comp++)
			{
				dX1[comp]=(centreNeighbour[comp]-centreCell[comp]);
				dX2[comp]=((facePoints[1][comp]+facePoints[2][comp])/2.0-(facePoints[0][comp]+facePoints[3][comp])/2.0);
				dX3[comp]=(facePoints[2][comp]+facePoints[3][comp])/2.0-(facePoints[0][comp]+facePoints[1][comp])/2.0;
				/* capital X represents Vector*/
			}
		break;	
	
		case 2:// South Face
			for(comp=0;comp<3;comp++)
			{
				dX1[comp]=((facePoints[1][comp]+facePoints[2][comp])/2.0-(facePoints[0][comp]+facePoints[3][comp])/2.0);
				dX2[comp]=(centreCell[comp]-centreNeighbour[comp]);
				dX3[comp]=(facePoints[2][comp]+facePoints[3][comp])/2.0-(facePoints[0][comp]+facePoints[1][comp])/2.0;
				/*capital X represents Vector*/
			}
		break;
	
		case 3:// North Face
			for(comp=0;comp<3;comp++)
			{
				dX1[comp]=((facePoints[1][comp]+facePoints[2][comp])/2.0-(facePoints[0][comp]+facePoints[3][comp])/2.0);
				dX2[comp]=(centreNeighbour[comp]-centreCell[comp]);
				dX3[comp]=(facePoints[2][comp]+facePoints[3][comp])/2.0-(facePoints[0][comp]+facePoints[1][comp])/2.0;
				/* capital X represents Vector*/
			}
		break;
	
		case 4:// Bottam Face
			for(comp=0;comp<3;comp++)
			{
				dX1[comp]=(facePoints[2][comp]+facePoints[3][comp])/2.0-(facePoints[0][comp]+facePoints[1][comp])/2.0;
				dX2[comp]=(facePoints[0][comp]+facePoints[3][comp])/2.0-(facePoints[1][comp]+facePoints[2][comp])/2.0;
				dX3[comp]=(centreCell[comp]-centreNeighbour[comp]);
				/* capital X represents Vector*/
			}
		break;
		
		case 5:// Top Face
			for(comp=0;comp<3;comp++)
			{
				dX1[comp]= (facePoints[2][comp]+facePoints[3][comp])/2.0-(facePoints[0][comp]+facePoints[1][comp])/2.0;
				dX2[comp]= (facePoints[0][comp]+facePoints[3][comp])/2.0-(facePoints[1][comp]+facePoints[2][comp])/2.0;
				dX3[comp]=(centreNeighbour[comp]-centreCell[comp]);
				/* capital X represents Vector*/
			}
		break;
	}//End of switch case
	
	/*finding the mangitudes*/
	magdX1=0.0;
	magdX2=0.0;
	magdX3=0.0;
	for(comp=0;comp<3;comp++)
	{
		magdX1=magdX1+dX1[comp]*dX1[comp];
		magdX2=magdX2+dX2[comp]*dX2[comp];
		magdX3=magdX3+dX3[comp]*dX3[comp];
	}
   	magdX1=sqrt(magdX1);
   	magdX2=sqrt(magdX2);
   	magdX3=sqrt(magdX3);
	/*finding the unit vectors */
	for(comp=0;comp<3;comp++)
	{  
		if(fabs(magdX1)<1e-10)
			n1[comp]=0.0;
		else    
			n1[comp]=dX1[comp]/magdX1;

		if(fabs(magdX2)<1e-10)
			n2[comp]=0.0;
		else	    
			n2[comp]=dX2[comp]/magdX2;

		if(fabs(magdX3)<1e-10)
			n3[comp]=0.0;
		else
			n3[comp]=dX3[comp]/magdX3;
	}   
	Det=n1[0]*(n2[1]*n3[2]-n2[2]*n3[1])-n1[1]*(n2[0]*n3[2]-n2[2]*n3[0])+n1[2]*(n2[0]*n3[1]-n2[1]*n3[0]);
	Det1=surf[0]*(n2[1]*n3[2]-n2[2]*n3[1])-surf[1]*(n2[0]*n3[2]-n2[2]*n3[0])+surf[2]*(n2[0]*n3[1]-n2[1]*n3[0]);
	Det2=n1[0]*(surf[1]*n3[2]-surf[2]*n3[1])-n1[1]*(surf[0]*n3[2]-surf[2]*n3[0])+n1[2]*(surf[0]*n3[1]-surf[1]*n3[0]);
	Det3=n1[0]*(n2[1]*surf[2]-n2[2]*surf[1])-n1[1]*(n2[0]*surf[2]-n2[2]*surf[0])+n1[2]*(n2[0]*surf[1]-n2[1]*surf[0]);
	/*if(Det<e-10)
	{
		alpha[0]=0.0;
		alpha[1]=0.0;
		alpha[2]=0.0;
	}
	else*/
		if(fabs(Det)<=1.0e-10)
		{
			alpha[0]=0.0;
			alpha[1]=0.0;
			alpha[2]=0.0;
		}
		else
		{
		alpha[0]=Det1/Det;
		alpha[1]=Det2/Det;
		alpha[2]=Det3/Det;
		}
		if(fabs(magdX1)<1e-10)
			beta[0]=0.0;
		else
			beta[0]=alpha[0]/magdX1;
		   
		if(fabs(magdX2)<1e-10)
			beta[1]=0.0;
		else
			beta[1]=alpha[1]/magdX2;
		
		if(fabs(magdX3)<1e-10)
			beta[2]=0.0;
		else
			beta[2]=alpha[2]/magdX3;
	return 0;
}
