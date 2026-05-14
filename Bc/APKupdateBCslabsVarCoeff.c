#include "APKupdateBCslabsVarCoeff.h"

/**
This function is used to calculate the boundaryface values of Variable coeffients like Density,Viscosity etc.
This function calculates the boundary face values of Variable coeffients according to the the function entered
by user in user defined functions described in the file BC_userDefinedVarCoefficients.c.
*/
int updateBCslabsVarCoeff(int variableTag,Zone *zone,ZoneBoundary *bcPtr,int dir,Base *base)
{
	int k,j,i,NI,NJ,NK;
	double *variable;
	int *ptrRange,index,temp[2][3];
	double *coordX,*coordY,*coordZ,surfPoints[4][3],centNeighbour[3];
	
	ptrRange = bcPtr->Point_Range;
	NI = *(zone->Size+0);
	NJ = *(zone->Size+1);
	NK = *(zone->Size+2);
	coordX=(zone->CgnsGrid)->X;
	coordY=(zone->CgnsGrid)->Y; 
	coordZ=(zone->CgnsGrid)->Z;
	
	variable = FindSlabPtrBC(bcPtr,variableTag);
	
	switch(variableTag)
	{
		case 710:
		{
			switch(dir)
			{
				case 1: /*Dir 1 ->Imax,East*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=(NI-2);
						temp[0][0]=NI*NJ;temp[0][1]=NI;temp[0][2]=1;temp[1][0]=k;temp[1][1]=j;temp[1][2]=i+1;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncDensity(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
						//printf("densty=%lf\n",*(variable+index));
						//exit(0);
					}
				break;

				case 2: /*Dir 2 ->Imin,West*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=0;
						temp[0][0]=NI*NJ;temp[0][1]=NI;temp[0][2]=1;temp[1][0]=k;temp[1][1]=j;temp[1][2]=i;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncDensity(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 3:/*Dir 3 ->Jmax,North*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=(NJ-2);
						temp[0][0]=NI*NJ;temp[0][1]=1;temp[0][2]=NI;temp[1][0]=k;temp[1][1]=i;temp[1][2]=j+1;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncDensity(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 4:/*Dir 4 ->Jmin,South*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=0;
						temp[0][0]=NI*NJ;temp[0][1]=1;temp[0][2]=NI;temp[1][0]=k;temp[1][1]=i;temp[1][2]=j;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncDensity(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 5:/*Dir 5 ->Kmax,Top*/
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=(NK-2);
						temp[0][0]=NI;temp[0][1]=1;temp[0][2]=NI*NJ;temp[1][0]=j;temp[1][1]=i;temp[1][2]=k+1;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncDensity(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 6:/*Dir 6 ->Kmin,Bottom*/
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=0;
						temp[0][0]=NI;temp[0][1]=1;temp[0][2]=NI*NJ;temp[1][0]=j;temp[1][1]=i;temp[1][2]=k;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncDensity(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;
				default:
					printf("Case not found,Terminating the Program from VDupdateBcslabVarCoeff\n");
					exit(1);
				break;
			}
// 			printf("The boundary name is %s in zone %s\n\n",bcPtr->BCName,zone->zonename);
// 			printf("The density at boundaries is %lf\n\n",*variable);
		}	
		break;
				
		case 711:
		{	switch(dir)
			{
				case 1: /*Dir 1 ->Imax,East*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=(NI-2);
						temp[0][0]=NI*NJ;temp[0][1]=NI;temp[0][2]=1;temp[1][0]=k;temp[1][1]=j;temp[1][2]=i+1;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncViscosity(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 2: /*Dir 2 ->Imin,West*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=0;
						temp[0][0]=NI*NJ;temp[0][1]=NI;temp[0][2]=1;temp[1][0]=k;temp[1][1]=j;temp[1][2]=i;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncViscosity(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 3:/*Dir 3 ->Jmax,North*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=(NJ-2);
						temp[0][0]=NI*NJ;temp[0][1]=1;temp[0][2]=NI;temp[1][0]=k;temp[1][1]=i;temp[1][2]=j+1;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncViscosity(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 4:/*Dir 4 ->Jmin,South*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=0;
						temp[0][0]=NI*NJ;temp[0][1]=1;temp[0][2]=NI;temp[1][0]=k;temp[1][1]=i;temp[1][2]=j;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncViscosity(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 5:/*Dir 5 ->Kmax,Top*/
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=(NK-2);
						temp[0][0]=NI;temp[0][1]=1;temp[0][2]=NI*NJ;temp[1][0]=j;temp[1][1]=i;temp[1][2]=k+1;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncViscosity(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 6:/*Dir 6 ->Kmin,Bottom*/
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=0;
						temp[0][0]=NI;temp[0][1]=1;temp[0][2]=NI*NJ;temp[1][0]=j;temp[1][1]=i;temp[1][2]=k;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncViscosity(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;
				default:
					printf("Case not found,Terminating the Program from VDupdateBcslabVarCoeff\n");
					exit(1);
				break;
			}
// 			printf("The boundary name is %s in zone %s\n\n",bcPtr->BCName,zone->zonename);
// 			printf("The Viscosity at boundaries is %lf\n\n",*variable);
		}
		break;			
					
		case 712:
		{	switch(dir)
			{
				case 1: /*Dir 1 ->Imax,East*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=(NI-2);
						temp[0][0]=NI*NJ;temp[0][1]=NI;temp[0][2]=1;temp[1][0]=k;temp[1][1]=j;temp[1][2]=i+1;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncThermalConductivity(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 2: /*Dir 2 ->Imin,West*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=0;
						temp[0][0]=NI*NJ;temp[0][1]=NI;temp[0][2]=1;temp[1][0]=k;temp[1][1]=j;temp[1][2]=i;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncThermalConductivity(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 3:/*Dir 3 ->Jmax,North*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=(NJ-2);
						temp[0][0]=NI*NJ;temp[0][1]=1;temp[0][2]=NI;temp[1][0]=k;temp[1][1]=i;temp[1][2]=j+1;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncThermalConductivity(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 4:/*Dir 4 ->Jmin,South*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=0;
						temp[0][0]=NI*NJ;temp[0][1]=1;temp[0][2]=NI;temp[1][0]=k;temp[1][1]=i;temp[1][2]=j;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncThermalConductivity(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 5:/*Dir 5 ->Kmax,Top*/
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=(NK-2);
						temp[0][0]=NI;temp[0][1]=1;temp[0][2]=NI*NJ;temp[1][0]=j;temp[1][1]=i;temp[1][2]=k+1;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncThermalConductivity(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 6:/*Dir 6 ->Kmin,Bottom*/
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=0;
						temp[0][0]=NI;temp[0][1]=1;temp[0][2]=NI*NJ;temp[1][0]=j;temp[1][1]=i;temp[1][2]=k;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncThermalConductivity(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;
				default:
					printf("Case not found,Terminating the Program from VDupdateBcslabVarCoeff\n");
					exit(1);
				break;
			}
		}	
		break;		
					
		case 713:
		{	switch(dir)
			{
				case 1: /*Dir 1 ->Imax,East*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=(NI-2);
						temp[0][0]=NI*NJ;temp[0][1]=NI;temp[0][2]=1;temp[1][0]=k;temp[1][1]=j;temp[1][2]=i+1;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncSpecificHeat(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 2: /*Dir 2 ->Imin,West*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=0;
						temp[0][0]=NI*NJ;temp[0][1]=NI;temp[0][2]=1;temp[1][0]=k;temp[1][1]=j;temp[1][2]=i;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncSpecificHeat(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 3:/*Dir 3 ->Jmax,North*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=(NJ-2);
						temp[0][0]=NI*NJ;temp[0][1]=1;temp[0][2]=NI;temp[1][0]=k;temp[1][1]=i;temp[1][2]=j+1;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncSpecificHeat(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 4:/*Dir 4 ->Jmin,South*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=0;
						temp[0][0]=NI*NJ;temp[0][1]=1;temp[0][2]=NI;temp[1][0]=k;temp[1][1]=i;temp[1][2]=j;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncSpecificHeat(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 5:/*Dir 5 ->Kmax,Top*/
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=(NK-2);
						temp[0][0]=NI;temp[0][1]=1;temp[0][2]=NI*NJ;temp[1][0]=j;temp[1][1]=i;temp[1][2]=k+1;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncSpecificHeat(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 6:/*Dir 6 ->Kmin,Bottom*/
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=0;
						temp[0][0]=NI;temp[0][1]=1;temp[0][2]=NI*NJ;temp[1][0]=j;temp[1][1]=i;temp[1][2]=k;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncSpecificHeat(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;
				default:
					printf("Case not found,Terminating the Program from VDupdateBcslabVarCoeff\n");
					exit(1);
				break;
			}
		}	
		break;
		
		case 717:
		{	switch(dir)
			{
				case 1: /*Dir 1 ->Imax,East*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=(NI-2);
						temp[0][0]=NI*NJ;temp[0][1]=NI;temp[0][2]=1;temp[1][0]=k;temp[1][1]=j;temp[1][2]=i+1;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncChargeDiffusion(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 2: /*Dir 2 ->Imin,West*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=0;
						temp[0][0]=NI*NJ;temp[0][1]=NI;temp[0][2]=1;temp[1][0]=k;temp[1][1]=j;temp[1][2]=i;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncChargeDiffusion(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 3:/*Dir 3 ->Jmax,North*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=(NJ-2);
						temp[0][0]=NI*NJ;temp[0][1]=1;temp[0][2]=NI;temp[1][0]=k;temp[1][1]=i;temp[1][2]=j+1;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncChargeDiffusion(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 4:/*Dir 4 ->Jmin,South*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=0;
						temp[0][0]=NI*NJ;temp[0][1]=1;temp[0][2]=NI;temp[1][0]=k;temp[1][1]=i;temp[1][2]=j;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncChargeDiffusion(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 5:/*Dir 5 ->Kmax,Top*/
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=(NK-2);
						temp[0][0]=NI;temp[0][1]=1;temp[0][2]=NI*NJ;temp[1][0]=j;temp[1][1]=i;temp[1][2]=k+1;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncChargeDiffusion(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 6:/*Dir 6 ->Kmin,Bottom*/
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=0;
						temp[0][0]=NI;temp[0][1]=1;temp[0][2]=NI*NJ;temp[1][0]=j;temp[1][1]=i;temp[1][2]=k;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncChargeDiffusion(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;
				default:
					printf("Case not found,Terminating the Program from VDupdateBcslabVarCoeff\n");
					exit(1);
				break;
			}
		}	
		break;
		
		case 700:
		{	
			if(base->flagTurb == 0 || base->BaseParameters->Turb->numTurbScalars < 1)
			{
				switch(dir)
				{
					case 1: /*Dir 1 ->Imax,East*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						{
							index = k*(NJ-1)+j;
							i=(NI-2);
							temp[0][0]=NI*NJ;temp[0][1]=NI;temp[0][2]=1;temp[1][0]=k;temp[1][1]=j;temp[1][2]=i+1;
							GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
							FindFaceCenterBC(surfPoints,centNeighbour);
							*(variable+index) = BCuserfuncScalCoeff1(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
						}
					break;
	
					case 2: /*Dir 2 ->Imin,West*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						{
							index = k*(NJ-1)+j;
							i=0;
							temp[0][0]=NI*NJ;temp[0][1]=NI;temp[0][2]=1;temp[1][0]=k;temp[1][1]=j;temp[1][2]=i;
							GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
							FindFaceCenterBC(surfPoints,centNeighbour);
							*(variable+index) = BCuserfuncScalCoeff1(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
						}
					break;
	
					case 3:/*Dir 3 ->Jmax,North*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = k*(NI-1)+i;
							j=(NJ-2);
							temp[0][0]=NI*NJ;temp[0][1]=1;temp[0][2]=NI;temp[1][0]=k;temp[1][1]=i;temp[1][2]=j+1;
							GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
							FindFaceCenterBC(surfPoints,centNeighbour);
							*(variable+index) = BCuserfuncScalCoeff1(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
						}
					break;
	
					case 4:/*Dir 4 ->Jmin,South*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = k*(NI-1)+i;
							j=0;
							temp[0][0]=NI*NJ;temp[0][1]=1;temp[0][2]=NI;temp[1][0]=k;temp[1][1]=i;temp[1][2]=j;
							GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
							FindFaceCenterBC(surfPoints,centNeighbour);
							*(variable+index) = BCuserfuncScalCoeff1(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
						}
					break;
	
					case 5:/*Dir 5 ->Kmax,Top*/
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = j*(NI-1)+i;
							k=(NK-2);
							temp[0][0]=NI;temp[0][1]=1;temp[0][2]=NI*NJ;temp[1][0]=j;temp[1][1]=i;temp[1][2]=k+1;
							GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
							FindFaceCenterBC(surfPoints,centNeighbour);
							*(variable+index) = BCuserfuncScalCoeff1(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
						}
					break;
	
					case 6:/*Dir 6 ->Kmin,Bottom*/
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = j*(NI-1)+i;
							k=0;
							temp[0][0]=NI;temp[0][1]=1;temp[0][2]=NI*NJ;temp[1][0]=j;temp[1][1]=i;temp[1][2]=k;
							GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
							FindFaceCenterBC(surfPoints,centNeighbour);
							*(variable+index) = BCuserfuncScalCoeff1(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
						}
					break;
					default:
						printf("Case not found,Terminating the Program from VDupdateBcslabVarCoeff\n");
						exit(1);
					break;
				}
			}
		}	
		break;
		
		case 701:
		{
			if(base->flagTurb == 0 || base->BaseParameters->Turb->numTurbScalars < 2)
			{
				switch(dir)
				{
					case 1: /*Dir 1 ->Imax,East*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						{
							index = k*(NJ-1)+j;
							i=(NI-2);
							temp[0][0]=NI*NJ;temp[0][1]=NI;temp[0][2]=1;temp[1][0]=k;temp[1][1]=j;temp[1][2]=i+1;
							GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
							FindFaceCenterBC(surfPoints,centNeighbour);
							*(variable+index) = BCuserfuncScalCoeff2(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
						}
					break;
	
					case 2: /*Dir 2 ->Imin,West*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						{
							index = k*(NJ-1)+j;
							i=0;
							temp[0][0]=NI*NJ;temp[0][1]=NI;temp[0][2]=1;temp[1][0]=k;temp[1][1]=j;temp[1][2]=i;
							GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
							FindFaceCenterBC(surfPoints,centNeighbour);
							*(variable+index) = BCuserfuncScalCoeff2(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
						}
					break;
	
					case 3:/*Dir 3 ->Jmax,North*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = k*(NI-1)+i;
							j=(NJ-2);
							temp[0][0]=NI*NJ;temp[0][1]=1;temp[0][2]=NI;temp[1][0]=k;temp[1][1]=i;temp[1][2]=j+1;
							GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
							FindFaceCenterBC(surfPoints,centNeighbour);
							*(variable+index) = BCuserfuncScalCoeff2(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
						}
					break;
	
					case 4:/*Dir 4 ->Jmin,South*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = k*(NI-1)+i;
							j=0;
							temp[0][0]=NI*NJ;temp[0][1]=1;temp[0][2]=NI;temp[1][0]=k;temp[1][1]=i;temp[1][2]=j;
							GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
							FindFaceCenterBC(surfPoints,centNeighbour);
							*(variable+index) = BCuserfuncScalCoeff2(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
						}
					break;
	
					case 5:/*Dir 5 ->Kmax,Top*/
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = j*(NI-1)+i;
							k=(NK-2);
							temp[0][0]=NI;temp[0][1]=1;temp[0][2]=NI*NJ;temp[1][0]=j;temp[1][1]=i;temp[1][2]=k+1;
							GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
							FindFaceCenterBC(surfPoints,centNeighbour);
							*(variable+index) = BCuserfuncScalCoeff2(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
						}
					break;
	
					case 6:/*Dir 6 ->Kmin,Bottom*/
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = j*(NI-1)+i;
							k=0;
							temp[0][0]=NI;temp[0][1]=1;temp[0][2]=NI*NJ;temp[1][0]=j;temp[1][1]=i;temp[1][2]=k;
							GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
							FindFaceCenterBC(surfPoints,centNeighbour);
							*(variable+index) = BCuserfuncScalCoeff2(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
						}
					break;
					default:
						printf("Case not found,Terminating the Program from VDupdateBcslabVarCoeff\n");
						exit(1);
					break;
				}		
			}
		}
		break;
		
		case 702:
		{
			if(base->flagTurb == 0 || base->BaseParameters->Turb->numTurbScalars < 3)
			{
				switch(dir)
				{
					case 1: /*Dir 1 ->Imax,East*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						{
							index = k*(NJ-1)+j;
							i=(NI-2);
							temp[0][0]=NI*NJ;temp[0][1]=NI;temp[0][2]=1;temp[1][0]=k;temp[1][1]=j;temp[1][2]=i+1;
							GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
							FindFaceCenterBC(surfPoints,centNeighbour);
							*(variable+index) = BCuserfuncScalCoeff3(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
						}
					break;
	
					case 2: /*Dir 2 ->Imin,West*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						{
							index = k*(NJ-1)+j;
							i=0;
							temp[0][0]=NI*NJ;temp[0][1]=NI;temp[0][2]=1;temp[1][0]=k;temp[1][1]=j;temp[1][2]=i;
							GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
							FindFaceCenterBC(surfPoints,centNeighbour);
							*(variable+index) = BCuserfuncScalCoeff3(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
						}
					break;
	
					case 3:/*Dir 3 ->Jmax,North*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = k*(NI-1)+i;
							j=(NJ-2);
							temp[0][0]=NI*NJ;temp[0][1]=1;temp[0][2]=NI;temp[1][0]=k;temp[1][1]=i;temp[1][2]=j+1;
							GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
							FindFaceCenterBC(surfPoints,centNeighbour);
							*(variable+index) = BCuserfuncScalCoeff3(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
						}
					break;
	
					case 4:/*Dir 4 ->Jmin,South*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = k*(NI-1)+i;
							j=0;
							temp[0][0]=NI*NJ;temp[0][1]=1;temp[0][2]=NI;temp[1][0]=k;temp[1][1]=i;temp[1][2]=j;
							GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
							FindFaceCenterBC(surfPoints,centNeighbour);
							*(variable+index) = BCuserfuncScalCoeff3(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
						}
					break;
	
					case 5:/*Dir 5 ->Kmax,Top*/
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = j*(NI-1)+i;
							k=(NK-2);
							temp[0][0]=NI;temp[0][1]=1;temp[0][2]=NI*NJ;temp[1][0]=j;temp[1][1]=i;temp[1][2]=k+1;
							GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
							FindFaceCenterBC(surfPoints,centNeighbour);
							*(variable+index) = BCuserfuncScalCoeff3(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
						}
					break;
	
					case 6:/*Dir 6 ->Kmin,Bottom*/
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = j*(NI-1)+i;
							k=0;
							temp[0][0]=NI;temp[0][1]=1;temp[0][2]=NI*NJ;temp[1][0]=j;temp[1][1]=i;temp[1][2]=k;
							GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
							FindFaceCenterBC(surfPoints,centNeighbour);
							*(variable+index) = BCuserfuncScalCoeff3(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
						}
					break;
					default:
						printf("Case not found,Terminating the Program from VDupdateBcslabVarCoeff\n");
						exit(1);
					break;
				}		
			}
		}	
		break;
		
		case 703:
		{
			if(base->flagTurb == 0 || base->BaseParameters->Turb->numTurbScalars < 4)
			{
				switch(dir)
				{
					case 1: /*Dir 1 ->Imax,East*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						{
							index = k*(NJ-1)+j;
							i=(NI-2);
							temp[0][0]=NI*NJ;temp[0][1]=NI;temp[0][2]=1;temp[1][0]=k;temp[1][1]=j;temp[1][2]=i+1;
							GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
							FindFaceCenterBC(surfPoints,centNeighbour);
							*(variable+index) = BCuserfuncScalCoeff4(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
						}
					break;
	
					case 2: /*Dir 2 ->Imin,West*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						{
							index = k*(NJ-1)+j;
							i=0;
							temp[0][0]=NI*NJ;temp[0][1]=NI;temp[0][2]=1;temp[1][0]=k;temp[1][1]=j;temp[1][2]=i;
							GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
							FindFaceCenterBC(surfPoints,centNeighbour);
							*(variable+index) = BCuserfuncScalCoeff4(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
						}
					break;
	
					case 3:/*Dir 3 ->Jmax,North*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = k*(NI-1)+i;
							j=(NJ-2);
							temp[0][0]=NI*NJ;temp[0][1]=1;temp[0][2]=NI;temp[1][0]=k;temp[1][1]=i;temp[1][2]=j+1;
							GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
							FindFaceCenterBC(surfPoints,centNeighbour);
							*(variable+index) = BCuserfuncScalCoeff4(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
						}
					break;
	
					case 4:/*Dir 4 ->Jmin,South*/
						for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = k*(NI-1)+i;
							j=0;
							temp[0][0]=NI*NJ;temp[0][1]=1;temp[0][2]=NI;temp[1][0]=k;temp[1][1]=i;temp[1][2]=j;
							GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
							FindFaceCenterBC(surfPoints,centNeighbour);
							*(variable+index) = BCuserfuncScalCoeff4(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
						}
					break;
	
					case 5:/*Dir 5 ->Kmax,Top*/
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = j*(NI-1)+i;
							k=(NK-2);
							temp[0][0]=NI;temp[0][1]=1;temp[0][2]=NI*NJ;temp[1][0]=j;temp[1][1]=i;temp[1][2]=k+1;
							GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
							FindFaceCenterBC(surfPoints,centNeighbour);
							*(variable+index) = BCuserfuncScalCoeff4(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
						}
					break;
	
					case 6:/*Dir 6 ->Kmin,Bottom*/
						for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
						for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
						{
							index = j*(NI-1)+i;
							k=0;
							temp[0][0]=NI;temp[0][1]=1;temp[0][2]=NI*NJ;temp[1][0]=j;temp[1][1]=i;temp[1][2]=k;
							GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
							FindFaceCenterBC(surfPoints,centNeighbour);
							*(variable+index) = BCuserfuncScalCoeff4(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
						}
					break;
					default:
						printf("Case not found,Terminating the Program from VDupdateBcslabVarCoeff\n");
						exit(1);
					break;
				}		
			}
		}	
		break;
		
		case 704:
		{	switch(dir)
			{
				case 1: /*Dir 1 ->Imax,East*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=(NI-2);
						temp[0][0]=NI*NJ;temp[0][1]=NI;temp[0][2]=1;temp[1][0]=k;temp[1][1]=j;temp[1][2]=i+1;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncScalCoeff5(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 2: /*Dir 2 ->Imin,West*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=0;
						temp[0][0]=NI*NJ;temp[0][1]=NI;temp[0][2]=1;temp[1][0]=k;temp[1][1]=j;temp[1][2]=i;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncScalCoeff5(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 3:/*Dir 3 ->Jmax,North*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=(NJ-2);
						temp[0][0]=NI*NJ;temp[0][1]=1;temp[0][2]=NI;temp[1][0]=k;temp[1][1]=i;temp[1][2]=j+1;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncScalCoeff5(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 4:/*Dir 4 ->Jmin,South*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=0;
						temp[0][0]=NI*NJ;temp[0][1]=1;temp[0][2]=NI;temp[1][0]=k;temp[1][1]=i;temp[1][2]=j;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncScalCoeff5(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 5:/*Dir 5 ->Kmax,Top*/
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=(NK-2);
						temp[0][0]=NI;temp[0][1]=1;temp[0][2]=NI*NJ;temp[1][0]=j;temp[1][1]=i;temp[1][2]=k+1;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncScalCoeff5(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 6:/*Dir 6 ->Kmin,Bottom*/
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=0;
						temp[0][0]=NI;temp[0][1]=1;temp[0][2]=NI*NJ;temp[1][0]=j;temp[1][1]=i;temp[1][2]=k;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncScalCoeff5(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;
				default:
					printf("Case not found,Terminating the Program from VDupdateBcslabVarCoeff\n");
					exit(1);
				break;
			}		
		}	
		break;
		
		case 705:
		{	switch(dir)
			{
				case 1: /*Dir 1 ->Imax,East*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=(NI-2);
						temp[0][0]=NI*NJ;temp[0][1]=NI;temp[0][2]=1;temp[1][0]=k;temp[1][1]=j;temp[1][2]=i+1;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncScalCoeff6(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 2: /*Dir 2 ->Imin,West*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=0;
						temp[0][0]=NI*NJ;temp[0][1]=NI;temp[0][2]=1;temp[1][0]=k;temp[1][1]=j;temp[1][2]=i;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncScalCoeff6(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 3:/*Dir 3 ->Jmax,North*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=(NJ-2);
						temp[0][0]=NI*NJ;temp[0][1]=1;temp[0][2]=NI;temp[1][0]=k;temp[1][1]=i;temp[1][2]=j+1;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncScalCoeff6(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 4:/*Dir 4 ->Jmin,South*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=0;
						temp[0][0]=NI*NJ;temp[0][1]=1;temp[0][2]=NI;temp[1][0]=k;temp[1][1]=i;temp[1][2]=j;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncScalCoeff6(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 5:/*Dir 5 ->Kmax,Top*/
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=(NK-2);
						temp[0][0]=NI;temp[0][1]=1;temp[0][2]=NI*NJ;temp[1][0]=j;temp[1][1]=i;temp[1][2]=k+1;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncScalCoeff6(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 6:/*Dir 6 ->Kmin,Bottom*/
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=0;
						temp[0][0]=NI;temp[0][1]=1;temp[0][2]=NI*NJ;temp[1][0]=j;temp[1][1]=i;temp[1][2]=k;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncScalCoeff6(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;
				default:
					printf("Case not found,Terminating the Program from VDupdateBcslabVarCoeff\n");
					exit(1);
				break;
			}		
		}	
		break;
		
		case 706:
		{	switch(dir)
			{
				case 1: /*Dir 1 ->Imax,East*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=(NI-2);
						temp[0][0]=NI*NJ;temp[0][1]=NI;temp[0][2]=1;temp[1][0]=k;temp[1][1]=j;temp[1][2]=i+1;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncScalCoeff7(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 2: /*Dir 2 ->Imin,West*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=0;
						temp[0][0]=NI*NJ;temp[0][1]=NI;temp[0][2]=1;temp[1][0]=k;temp[1][1]=j;temp[1][2]=i;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncScalCoeff7(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 3:/*Dir 3 ->Jmax,North*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=(NJ-2);
						temp[0][0]=NI*NJ;temp[0][1]=1;temp[0][2]=NI;temp[1][0]=k;temp[1][1]=i;temp[1][2]=j+1;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncScalCoeff7(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 4:/*Dir 4 ->Jmin,South*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=0;
						temp[0][0]=NI*NJ;temp[0][1]=1;temp[0][2]=NI;temp[1][0]=k;temp[1][1]=i;temp[1][2]=j;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncScalCoeff7(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 5:/*Dir 5 ->Kmax,Top*/
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=(NK-2);
						temp[0][0]=NI;temp[0][1]=1;temp[0][2]=NI*NJ;temp[1][0]=j;temp[1][1]=i;temp[1][2]=k+1;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncScalCoeff7(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 6:/*Dir 6 ->Kmin,Bottom*/
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=0;
						temp[0][0]=NI;temp[0][1]=1;temp[0][2]=NI*NJ;temp[1][0]=j;temp[1][1]=i;temp[1][2]=k;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncScalCoeff7(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;
				default:
					printf("Case not found,Terminating the Program from VDupdateBcslabVarCoeff\n");
					exit(1);
				break;
			}		
		}	
		break;
		
		case 707:
		{	switch(dir)
			{
				case 1: /*Dir 1 ->Imax,East*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=(NI-2);
						temp[0][0]=NI*NJ;temp[0][1]=NI;temp[0][2]=1;temp[1][0]=k;temp[1][1]=j;temp[1][2]=i+1;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncScalCoeff8(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 2: /*Dir 2 ->Imin,West*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=0;
						temp[0][0]=NI*NJ;temp[0][1]=NI;temp[0][2]=1;temp[1][0]=k;temp[1][1]=j;temp[1][2]=i;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncScalCoeff8(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 3:/*Dir 3 ->Jmax,North*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=(NJ-2);
						temp[0][0]=NI*NJ;temp[0][1]=1;temp[0][2]=NI;temp[1][0]=k;temp[1][1]=i;temp[1][2]=j+1;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncScalCoeff8(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 4:/*Dir 4 ->Jmin,South*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=0;
						temp[0][0]=NI*NJ;temp[0][1]=1;temp[0][2]=NI;temp[1][0]=k;temp[1][1]=i;temp[1][2]=j;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncScalCoeff8(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 5:/*Dir 5 ->Kmax,Top*/
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=(NK-2);
						temp[0][0]=NI;temp[0][1]=1;temp[0][2]=NI*NJ;temp[1][0]=j;temp[1][1]=i;temp[1][2]=k+1;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncScalCoeff8(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 6:/*Dir 6 ->Kmin,Bottom*/
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=0;
						temp[0][0]=NI;temp[0][1]=1;temp[0][2]=NI*NJ;temp[1][0]=j;temp[1][1]=i;temp[1][2]=k;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncScalCoeff8(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;
				default:
					printf("Case not found,Terminating the Program from VDupdateBcslabVarCoeff\n");
					exit(1);
				break;
			}		
		}	
		break;
		
		case 708:
		{	switch(dir)
			{
				case 1: /*Dir 1 ->Imax,East*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=(NI-2);
						temp[0][0]=NI*NJ;temp[0][1]=NI;temp[0][2]=1;temp[1][0]=k;temp[1][1]=j;temp[1][2]=i+1;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncScalCoeff9(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 2: /*Dir 2 ->Imin,West*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=0;
						temp[0][0]=NI*NJ;temp[0][1]=NI;temp[0][2]=1;temp[1][0]=k;temp[1][1]=j;temp[1][2]=i;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncScalCoeff9(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 3:/*Dir 3 ->Jmax,North*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=(NJ-2);
						temp[0][0]=NI*NJ;temp[0][1]=1;temp[0][2]=NI;temp[1][0]=k;temp[1][1]=i;temp[1][2]=j+1;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncScalCoeff9(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 4:/*Dir 4 ->Jmin,South*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=0;
						temp[0][0]=NI*NJ;temp[0][1]=1;temp[0][2]=NI;temp[1][0]=k;temp[1][1]=i;temp[1][2]=j;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncScalCoeff9(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 5:/*Dir 5 ->Kmax,Top*/
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=(NK-2);
						temp[0][0]=NI;temp[0][1]=1;temp[0][2]=NI*NJ;temp[1][0]=j;temp[1][1]=i;temp[1][2]=k+1;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncScalCoeff9(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 6:/*Dir 6 ->Kmin,Bottom*/
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=0;
						temp[0][0]=NI;temp[0][1]=1;temp[0][2]=NI*NJ;temp[1][0]=j;temp[1][1]=i;temp[1][2]=k;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncScalCoeff9(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;
				default:
					printf("Case not found,Terminating the Program from VDupdateBcslabVarCoeff\n");
					exit(1);
				break;
			}		
		}	
		break;
		
		case 709:
		{	switch(dir)
			{
				case 1: /*Dir 1 ->Imax,East*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=(NI-2);
						temp[0][0]=NI*NJ;temp[0][1]=NI;temp[0][2]=1;temp[1][0]=k;temp[1][1]=j;temp[1][2]=i+1;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncScalCoeff10(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 2: /*Dir 2 ->Imin,West*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=0;
						temp[0][0]=NI*NJ;temp[0][1]=NI;temp[0][2]=1;temp[1][0]=k;temp[1][1]=j;temp[1][2]=i;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncScalCoeff10(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 3:/*Dir 3 ->Jmax,North*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=(NJ-2);
						temp[0][0]=NI*NJ;temp[0][1]=1;temp[0][2]=NI;temp[1][0]=k;temp[1][1]=i;temp[1][2]=j+1;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncScalCoeff10(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 4:/*Dir 4 ->Jmin,South*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=0;
						temp[0][0]=NI*NJ;temp[0][1]=1;temp[0][2]=NI;temp[1][0]=k;temp[1][1]=i;temp[1][2]=j;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncScalCoeff10(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 5:/*Dir 5 ->Kmax,Top*/
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=(NK-2);
						temp[0][0]=NI;temp[0][1]=1;temp[0][2]=NI*NJ;temp[1][0]=j;temp[1][1]=i;temp[1][2]=k+1;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncScalCoeff10(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 6:/*Dir 6 ->Kmin,Bottom*/
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=0;
						temp[0][0]=NI;temp[0][1]=1;temp[0][2]=NI*NJ;temp[1][0]=j;temp[1][1]=i;temp[1][2]=k;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncScalCoeff10(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;
				default:
					printf("Case not found,Terminating the Program from VDupdateBcslabVarCoeff\n");
					exit(1);
				break;
			}
		}	
		break;
		
		//for multiphase flow added by sushil
		case 720:
		
		case 721:
		
		case 722:
		
		case 723:
		
		case 724:
		
		break;
		//end
		case 725:
		{	switch(dir)
			{
				case 1: /*Dir 1 ->Imax,East*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=(NI-2);
						temp[0][0]=NI*NJ;temp[0][1]=NI;temp[0][2]=1;temp[1][0]=k;temp[1][1]=j;temp[1][2]=i+1;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncCPdrop(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 2: /*Dir 2 ->Imin,West*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					{
						index = k*(NJ-1)+j;
						i=0;
						temp[0][0]=NI*NJ;temp[0][1]=NI;temp[0][2]=1;temp[1][0]=k;temp[1][1]=j;temp[1][2]=i;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncCPdrop(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 3:/*Dir 3 ->Jmax,North*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=(NJ-2);
						temp[0][0]=NI*NJ;temp[0][1]=1;temp[0][2]=NI;temp[1][0]=k;temp[1][1]=i;temp[1][2]=j+1;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncCPdrop(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 4:/*Dir 4 ->Jmin,South*/
					for(k= *(ptrRange+2)-1;k < *(ptrRange+5)-1;k++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = k*(NI-1)+i;
						j=0;
						temp[0][0]=NI*NJ;temp[0][1]=1;temp[0][2]=NI;temp[1][0]=k;temp[1][1]=i;temp[1][2]=j;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncCPdrop(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 5:/*Dir 5 ->Kmax,Top*/
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=(NK-2);
						temp[0][0]=NI;temp[0][1]=1;temp[0][2]=NI*NJ;temp[1][0]=j;temp[1][1]=i;temp[1][2]=k+1;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncCPdrop(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;

				case 6:/*Dir 6 ->Kmin,Bottom*/
					for(j= *(ptrRange+1)-1;j < *(ptrRange+4)-1;j++)
					for(i= *(ptrRange+0)-1;i < *(ptrRange+3)-1;i++)
					{
						index = j*(NI-1)+i;
						k=0;
						temp[0][0]=NI;temp[0][1]=1;temp[0][2]=NI*NJ;temp[1][0]=j;temp[1][1]=i;temp[1][2]=k;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						FindFaceCenterBC(surfPoints,centNeighbour);
						*(variable+index) = BCuserfuncCPdrop(bcPtr,index,centNeighbour[0],centNeighbour[1],centNeighbour[2],base);
					}
				break;
				default:
					printf("Case not found,Terminating the Program from VDupdateBcslabVarCoeff\n");
					exit(1);
				break;
			}
		}	
		break;
		default:
			printf("Vartag [%d] not found, Terminating the program from updateBcslabsVarCoeff()\n", variableTag);
			exit(1);
		break;
	}
	return 0;
}

