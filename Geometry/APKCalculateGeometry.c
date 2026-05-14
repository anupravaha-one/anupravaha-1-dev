#include "APKCalculateGeometry.h"

/**
CalculateGeometry() is used to calculate the cell face surface areas, cell volumes, and diffusion lenghts for 
a Base.
 
The functions loops over all the Zones, all the cells in each Zone and over all the faces of each cell. 
Inside the face loop, the program is split into six cases for each of the six faces of a cell. In the each of 
the case statements two main tasks are performed - finding the coordinates of vertices which constitute a cell 
face and finding the adjacent cell center. The coordinates of the face vertices are found using the function 
GetSurfacePoints(). For finding the neighboring cell center we need to know whether the neighboring cell lies 
in the present Zone or in the Neighboring Zone across the inerface or that's a boundary. If it lies in the 
neighboring Zone or at boundary, first we need to find whether it lies across an interface or at a boundary 
and if across an interface the interface pointer needs to be found. The interface pointer is found using the 
function FindInterfaceBCPtrMod() along with a character tag which says whether Its a boundary or an interface. 
If the current cell lies near a boundary, the neighboring cell center is essentially the center of the face 
that lies on the boundary. In such a case the face center is found using FindFaceCenterBC(). In case the current 
cell lies near an interface the procedure for finding the neighboring cell center is listed below-Find the 
pointer of the neighboring Zone using FindZonePtr(). Find the indices of the neighboring cell and transform 
the indices into the convention of neighboring Zone using TransformCellIndices1(). Using the transformed indices 
and the neighboring Zone pointer, find the addresses (pointers) of the coordinates of the vertex which has the 
minimum indices(same as the cell indices) in the neighboring cell. Using these pointers find the coordinates of 
the neighboring cell center using FindCellCentre(). 

In case the neighboring cell lies inside the present Zone only the cell center coordinates are found using 
FindCellCentre() by supplying it the indices of the neighboring cell and pointers of the coordinates of the   
vertex with minimum indices in the cell. 

After finding the neighboring cell centers, cell face area vectors are found for each face of a cell inside 
the face loop using CalSurf(). CheckOutward() is used to ensure the direction of the surface area vector is 
outward. If the direction is inward the vector is multiplied by (-1) to make it outward. 

FindDeltaX() is used to find the cell center distances across a cell face in the face loop.

The diffusion lenghts are found using FindBeta() inside the face loop only. 

After termination of the face loop (still inside the cell loop) cell volumes are calculated as follows-
Retrieving the surface area vectors of three adjacent cell faces. 
Finding the vector sum of the surface area vectors found above. 
Finding the diagonal vector which passes through the intersection of the three adjacent faces.
Finding the dot product of the diagonal and the vector sum of the surface vectors. 
The cell volume is found by dividing the dot product by 3.0.
*/

int CalculateGeometry(Base *CgnsBase)
{
// 	FILE *r;
// 	r=fopen("SOLUTIONS/zoneinfo.dat","w");
	
	int i,j,k,comp,face,*ptrSize,*sizePtrNeighbor,index[3];
	int indices1[3],indices2[3],temp[2][3],nt,np,slabIndex,NumslabCells;
 
	double surfPoints[4][3],volume,adjacentFace[6][3],principalDiag[3],areaSum[3],surfVector[3],centroid[3],areacof;
  
	double *ptrPoint1X,*ptrPoint1Y,*ptrPoint1Z,beta[3],delta_X[3];
  
	double *ptrNeighbour1X,*ptrNeighbour1Y,*ptrNeighbour1Z,centNeighbour[3];
	
	double *coordX,*coordY,*coordZ,X[6],Y[6],Z[6],del[9];

	char tag,name[33];
	double *ptrPoint[3];
	Zone *ptrZone;
	void *facePtr;
	int NI,NJ,NK,nZones,numZones,Tag=1; 
	int numBC,nbc,numIface,nIface;
	ZoneBoundary *bcPtr;
	Interface *IfacePtr;
	
	nt=CgnsBase->BaseParameters->ntheta;
	np=CgnsBase->BaseParameters->nphi;

	nZones=((CgnsBase)->Number_of_zones);
// 	fprintf(r,"num of Zones=%d\n",nZones); 
	if(CgnsBase->BaseParameters->AxisymAzimuthangle!=NULL)
	{
		areacof=(8.0*atan(1.0))/(*(CgnsBase->BaseParameters->AxisymAzimuthangle));
	}
	else
	areacof=1.0;
	for(numZones=0;numZones<nZones;numZones++)
	{
// 		fprintf(r,"I am in Zone=%s\n",(CgnsBase->CgnsZone+numZones)->zonename); 
		((CgnsBase->CgnsZone+numZones)->ZGeom) =(ZoneGeometry *)malloc(sizeof(ZoneGeometry));
		ptrSize=((CgnsBase)->CgnsZone+numZones)->Size;
		
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
	
		(((CgnsBase->CgnsZone+numZones)->ZGeom)->SurfaceArea) = (double *)malloc((NI-1)*(NJ-1)*(NK-1)*6*3*sizeof(double));
		(((CgnsBase->CgnsZone+numZones)->ZGeom)->Beta) = (double *)malloc((NI-1)*(NJ-1)*(NK-1)*6*3*sizeof(double));
		(((CgnsBase->CgnsZone+numZones)->ZGeom)->deltaX) = (double *)malloc((NI-1)*(NJ-1)*(NK-1)*6*3*sizeof(double));
		(((CgnsBase->CgnsZone+numZones)->ZGeom)->Volume) = (double *)malloc((NI-1)*(NJ-1)*(NK-1)*sizeof(double));
		
		if(CgnsBase->BaseParameters->flagRadiation==1||CgnsBase->BaseParameters->flagflowRadiation==1)
		{
			(((CgnsBase->CgnsZone+numZones)->ZGeom)->Omega)= (double *)malloc(nt*np*sizeof(double));
			(((CgnsBase->CgnsZone+numZones)->ZGeom)->VectorX)=(double *)malloc(nt*np*sizeof(double));
			(((CgnsBase->CgnsZone+numZones)->ZGeom)->VectorY)=(double *)malloc(nt*np*sizeof(double));
			(((CgnsBase->CgnsZone+numZones)->ZGeom)->VectorZ)=(double *)malloc(nt*np*sizeof(double));
			
			(((CgnsBase->CgnsZone+numZones)->ZGeom)->DmF)= (double *)malloc((NI-1)*(NJ-1)*(NK-1)*6*nt*np*sizeof(double));
			(((CgnsBase->CgnsZone+numZones)->ZGeom)->NormalX)= (double *)malloc((NI-1)*(NJ-1)*(NK-1)*6*nt*np*sizeof(double));
			(((CgnsBase->CgnsZone+numZones)->ZGeom)->NormalY)= (double *)malloc((NI-1)*(NJ-1)*(NK-1)*6*nt*np*sizeof(double));
			(((CgnsBase->CgnsZone+numZones)->ZGeom)->NormalZ)= (double *)malloc((NI-1)*(NJ-1)*(NK-1)*6*nt*np*sizeof(double));
			if((((CgnsBase->CgnsZone+numZones)->ZGeom)->DmF)==NULL)
				printf("\nUnable to allocate memory\n");
		}
		
		if(CgnsBase->BaseParameters->flag2d == 0)//3D Case
			(((CgnsBase->CgnsZone+numZones)->ZGeom)->DeltaInverse) = (double *)malloc((NI-1)*(NJ-1)*(NK-1)*9*sizeof(double));
		else //2Dcase
			(((CgnsBase->CgnsZone+numZones)->ZGeom)->DeltaInverse) = (double *)malloc((NI-1)*(NJ-1)*(NK-1)*4*sizeof(double));
		
		if(CgnsBase->flagTurb == 1)
		{
			numBC =  (CgnsBase->CgnsZone+numZones)->Number_of_BC_Patches;
			for(nbc=0;nbc<numBC;nbc++)
			{
				//printf("Enter boundary tag--The boundary name is %s\n ",((CgnsBase->CgnsZone+numZones)->ZoneBc+nbc)->BCName);
				//scanf("%d",&((CgnsBase->CgnsZone+numZones)->ZoneBc+nbc)->Boundary_Tag);
				printf("boundarytag = %d\n",((CgnsBase->CgnsZone+numZones)->ZoneBc+nbc)->Boundary_Tag);
				if(((CgnsBase->CgnsZone+numZones)->ZoneBc+nbc)->Boundary_Tag == 22 && (CgnsBase->CgnsZone+numZones)->flagSolidFluid == 0)
				{
					((CgnsBase->CgnsZone+numZones)->ZoneBc+nbc)->deltaP = (double *)initializeslab(1,((CgnsBase->CgnsZone+numZones)->ZoneBc+nbc)->Direction,((CgnsBase->CgnsZone+numZones)->ZoneBc+nbc)->Point_Range);
					((CgnsBase->CgnsZone+numZones)->ZoneBc+nbc)->uTau = (double *)initializeslab(1,((CgnsBase->CgnsZone+numZones)->ZoneBc+nbc)->Direction,((CgnsBase->CgnsZone+numZones)->ZoneBc+nbc)->Point_Range);
					if(((CgnsBase->CgnsZone+numZones)->ZoneBc+nbc)->deltaP == NULL)
					{
						printf("Memory Allocation failed for deltaP,exiting from CalculateGeometry\n");
						exit(1);
					}
				}
				else
				{
					printf("**** boundarytag = %d ... deltaP -> NULL, uTau -> NULL ****\n",((CgnsBase->CgnsZone+numZones)->ZoneBc+nbc)->Boundary_Tag);
					((CgnsBase->CgnsZone+numZones)->ZoneBc+nbc)->deltaP = NULL;
					((CgnsBase->CgnsZone+numZones)->ZoneBc+nbc)->uTau = NULL;
				}
			}
			
			numIface =  (CgnsBase->CgnsZone+numZones)->Number_of_Interfaces;
			for(nIface=0;nIface<numIface;nIface++)
			{
				IfacePtr = (CgnsBase->CgnsZone+numZones)->Iface+nIface;
				ptrZone = (Zone *)FindZonePtr(IfacePtr->Adjoining_Zone_Name,CgnsBase);
				
				if(ptrZone->flagSolidFluid == 1 && (CgnsBase->CgnsZone+numZones)->flagSolidFluid == 0)
				{
					IfacePtr->deltaP = (double *)initializeslab(1,IfacePtr->Direction,IfacePtr->Point_Range);
					IfacePtr->uTau = (double *)initializeslab(1,IfacePtr->Direction,IfacePtr->Point_Range);
					if(IfacePtr->deltaP == NULL || IfacePtr->uTau == NULL)
					{
						printf("Memory Allocation failed for deltaP or uTau @ Interface Direction = %d, Zone =%s,exiting from CalculateGeometry\n", IfacePtr->Direction, (CgnsBase->CgnsZone+numZones)->zonename);
						exit(1);
					}
					printf("**** Interface Direction = %d Zone =%s ... Memory is allocated for deltaP, uTau  ****\n",IfacePtr->Direction, (CgnsBase->CgnsZone+numZones)->zonename);
				}
				else
				{
					printf("**** Interface Direction = %d Zone =%s ... deltaP -> NULL, uTau -> NULL ****\n",IfacePtr->Direction, (CgnsBase->CgnsZone+numZones)->zonename);
					IfacePtr->deltaP = NULL;
					IfacePtr->uTau = NULL;
				}
			}
		}
		
		if(CgnsBase->flagWf == 1)
		{
			numBC =  (CgnsBase->CgnsZone+numZones)->Number_of_BC_Patches;
			for(nbc=0;nbc<numBC;nbc++)
			{
				//printf("Enter boundary tag--The boundary name is %s\n ",((CgnsBase->CgnsZone+numZones)->ZoneBc+nbc)->BCName);
				//scanf("%d",&((CgnsBase->CgnsZone+numZones)->ZoneBc+nbc)->Boundary_Tag);
				printf("boundarytag = %d\n",((CgnsBase->CgnsZone+numZones)->ZoneBc+nbc)->Boundary_Tag);
				if(((CgnsBase->CgnsZone+numZones)->ZoneBc+nbc)->Boundary_Tag == 22 && (CgnsBase->CgnsZone+numZones)->flagSolidFluid == 0)
				{
					((CgnsBase->CgnsZone+numZones)->ZoneBc+nbc)->deltaP = (double *)initializeslab(1,((CgnsBase->CgnsZone+numZones)->ZoneBc+nbc)->Direction,((CgnsBase->CgnsZone+numZones)->ZoneBc+nbc)->Point_Range);
					
					if(((CgnsBase->CgnsZone+numZones)->ZoneBc+nbc)->deltaP == NULL)
					{
						printf("Memory Allocation failed for deltaP,exiting from CalculateGeometry\n");
						exit(1);
					}
				}
				else
				{
					printf("**** boundarytag = %d ... deltaP -> NULL ****\n",((CgnsBase->CgnsZone+numZones)->ZoneBc+nbc)->Boundary_Tag);
					((CgnsBase->CgnsZone+numZones)->ZoneBc+nbc)->deltaP = NULL;
					
				}
			}
			
			numIface =  (CgnsBase->CgnsZone+numZones)->Number_of_Interfaces;
			for(nIface=0;nIface<numIface;nIface++)
			{
				IfacePtr = (CgnsBase->CgnsZone+numZones)->Iface+nIface;
				ptrZone = (Zone *)FindZonePtr(IfacePtr->Adjoining_Zone_Name,CgnsBase);
				
				if(ptrZone->flagSolidFluid == 1 && (CgnsBase->CgnsZone+numZones)->flagSolidFluid == 0)
				{
					IfacePtr->deltaP = (double *)initializeslab(1,IfacePtr->Direction,IfacePtr->Point_Range);
				
					if(IfacePtr->deltaP == NULL)
					{
						printf("Memory Allocation failed for deltaP Interface Direction = %d, Zone =%s,exiting from CalculateGeometry\n", IfacePtr->Direction, (CgnsBase->CgnsZone+numZones)->zonename);
						exit(1);
					}
					printf("**** Interface Direction = %d Zone =%s ... Memory is allocated for deltaP  ****\n",IfacePtr->Direction, (CgnsBase->CgnsZone+numZones)->zonename);
				}
				else
				{
					printf("**** Interface Direction = %d Zone =%s ... deltaP -> NULL, ****\n",IfacePtr->Direction, (CgnsBase->CgnsZone+numZones)->zonename);
					IfacePtr->deltaP = NULL;
					
				}
			}
		}
		
		numBC =  (CgnsBase->CgnsZone+numZones)->Number_of_BC_Patches;
		for(nbc=0;nbc<numBC;nbc++)
		{
			bcPtr = (CgnsBase->CgnsZone+numZones)->ZoneBc+nbc;
			
			switch(bcPtr->Direction)
			{
				case 1:
				case 2:
					NumslabCells = (NJ-1)*(NK-1);
					break;
				case 3:
				case 4:
					NumslabCells = (NI-1)*(NK-1);
					break;
				case 5:
				case 6:
					NumslabCells = (NI-1)*(NJ-1);
					break;
			}
			bcPtr->FaceCenter = (double *)malloc(NumslabCells*3*sizeof(double));
		}
		
		coordX=((CgnsBase->CgnsZone+numZones)->CgnsGrid)->X;
		coordY=((CgnsBase->CgnsZone+numZones)->CgnsGrid)->Y; 
		coordZ=((CgnsBase->CgnsZone+numZones)->CgnsGrid)->Z;
		for(k=0;k<(NK-1);k++)
		for(j=0;j<(NJ-1);j++)
		for(i=0;i<(NI-1);i++)
		{
			ptrPoint1X=coordX+NI*NJ*k+NI*j+i;
			ptrPoint1Y=coordY+NI*NJ*k+NI*j+i;
			ptrPoint1Z=coordZ+NI*NJ*k+NI*j+i;/*point1 is the point with minimum indices in a cell i.e (i,j,k)*/
			FindCellCentre(ptrSize,ptrPoint1X,ptrPoint1Y,ptrPoint1Z,centroid);
			for(face=0;face<6;face++)
			{
				tag = 'n';
				switch(face)
				{
					case 0:/*Imin->West Face*/
						
						temp[0][0]=NI*NJ,temp[0][1]=NI,temp[0][2]=1,temp[1][0]=k,temp[1][1]=j,temp[1][2]=i;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						tag = 'a';//a default value
						if(i==0)
						{
							facePtr=(void *)FindInterfaceBCPtrMod((CgnsBase->CgnsZone+numZones),2,name,&tag);
							if(tag=='i')
							{
								facePtr = (Interface *)facePtr;
								ptrZone=(Zone *)FindZonePtr(name,CgnsBase);/* working properly*/
								Tag=strcmp(((CgnsBase->CgnsZone+numZones)->SolidFluidName),(ptrZone->SolidFluidName));
								if(Tag!=0)/*Done By Ashish */
								{
									
									FindFaceCenterBC(surfPoints,centNeighbour);
								}
								else
								{
									sizePtrNeighbor=ptrZone->Size;/*tested working well*/
									indices1[0]=i,indices1[1]=j+1,indices1[2]=k+1;
									TransformCellIndices1(indices1,facePtr,indices2);
									AdjoiningZonePointerXYZ(ptrZone,indices2,ptrPoint);
									ptrNeighbour1X=ptrPoint[0];
									ptrNeighbour1Y=ptrPoint[1];
									ptrNeighbour1Z=ptrPoint[2];
									FindCellCentre(sizePtrNeighbor,ptrNeighbour1X,ptrNeighbour1Y,ptrNeighbour1Z,centNeighbour);
								}
							}
							else
							if(tag=='b')
								FindFaceCenterBC(surfPoints,centNeighbour);
						}
						else
						{
							ptrNeighbour1X=ptrPoint1X-1;
							ptrNeighbour1Y=ptrPoint1Y-1;
							ptrNeighbour1Z=ptrPoint1Z-1;
							FindCellCentre(ptrSize,ptrNeighbour1X,ptrNeighbour1Y,ptrNeighbour1Z,centNeighbour);
						}
						X[0]=centNeighbour[0];
						Y[0]=centNeighbour[1];
						Z[0]=centNeighbour[2];
						break;
					
					case 1:/*Imax->East Face*/
						//printf("i m in case one\n");
						temp[0][0]=NI*NJ,temp[0][1]=NI,temp[0][2]=1,temp[1][0]=k,temp[1][1]=j,temp[1][2]=i+1;  
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						tag = 'a';//a default value
						if(i==NI-2)
						{
							facePtr=(void *)FindInterfaceBCPtrMod((CgnsBase->CgnsZone+numZones),1,name,&tag);
							if(tag=='i')
							{
								facePtr = (Interface *)facePtr;
								ptrZone=(Zone *)FindZonePtr(name,CgnsBase);
								Tag=strcmp(((CgnsBase->CgnsZone+numZones)->SolidFluidName),(ptrZone->SolidFluidName));
								if(Tag!=0 )/*Done By Ashish */
								{
									
									FindFaceCenterBC(surfPoints,centNeighbour);
								}
								else
								{
									sizePtrNeighbor=ptrZone->Size;
									indices1[0]=i+2,indices1[1]=j+1,indices1[2]=k+1;
									TransformCellIndices1(indices1,facePtr,indices2);
									AdjoiningZonePointerXYZ(ptrZone,indices2,ptrPoint);
									ptrNeighbour1X=ptrPoint[0];
									ptrNeighbour1Y=ptrPoint[1];
									ptrNeighbour1Z=ptrPoint[2];
		
									FindCellCentre(sizePtrNeighbor,ptrNeighbour1X,ptrNeighbour1Y,ptrNeighbour1Z,centNeighbour);
								}
							}
							else
							if(tag=='b')
								FindFaceCenterBC(surfPoints,centNeighbour);
						}
						else
						{
							ptrNeighbour1X=ptrPoint1X+1;
							ptrNeighbour1Y=ptrPoint1Y+1;
							ptrNeighbour1Z=ptrPoint1Z+1;
							FindCellCentre(ptrSize,ptrNeighbour1X,ptrNeighbour1Y,ptrNeighbour1Z,centNeighbour);
						}
						X[1]=centNeighbour[0];
						Y[1]=centNeighbour[1];
						Z[1]=centNeighbour[2];
						break;
			
					case 2:/*Jmin->South face*/
						temp[0][0]=NI*NJ,temp[0][1]=1,temp[0][2]=NI,temp[1][0]=k,temp[1][1]=i,temp[1][2]=j;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						tag = 'a';//a default value
						if(j==0)
						{
							facePtr=(void *)FindInterfaceBCPtrMod((CgnsBase->CgnsZone+numZones),4,name,&tag);
							if(tag=='i')
							{
								facePtr = (Interface *)facePtr;
								ptrZone=(Zone *)FindZonePtr(name,CgnsBase);
								Tag=strcmp(((CgnsBase->CgnsZone+numZones)->SolidFluidName),(ptrZone->SolidFluidName));
								if(Tag!=0 )/*Done By Ashish */
								{
									
									FindFaceCenterBC(surfPoints,centNeighbour);
								}
								else
								{
									sizePtrNeighbor=ptrZone->Size;
									indices1[0]=i+1,indices1[1]=j,indices1[2]=k+1;
									TransformCellIndices1(indices1,facePtr,indices2);
									AdjoiningZonePointerXYZ(ptrZone,indices2,ptrPoint);
									ptrNeighbour1X=ptrPoint[0];
									ptrNeighbour1Y=ptrPoint[1];
									ptrNeighbour1Z=ptrPoint[2];
		
									FindCellCentre(sizePtrNeighbor,ptrNeighbour1X,ptrNeighbour1Y,ptrNeighbour1Z,centNeighbour);
								}
							}
							else
							if(tag=='b')
								FindFaceCenterBC(surfPoints,centNeighbour);
						}
						else
						{ 
							ptrNeighbour1X=ptrPoint1X-NI;
							ptrNeighbour1Y=ptrPoint1Y-NI;
							ptrNeighbour1Z=ptrPoint1Z-NI;
							FindCellCentre(ptrSize,ptrNeighbour1X,ptrNeighbour1Y,ptrNeighbour1Z,centNeighbour);
						}
						X[2]=centNeighbour[0];
						Y[2]=centNeighbour[1];
						Z[2]=centNeighbour[2];
						break;
					case 3:/*Jmax->North Face*/
						temp[0][0]=NI*NJ,temp[0][1]=1,temp[0][2]=NI,temp[1][0]=k,temp[1][1]=i,temp[1][2]=j+1;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						tag = 'a';//a default value
						if(j==NJ-2)
						{
							facePtr=(void *)FindInterfaceBCPtrMod((CgnsBase->CgnsZone+numZones),3,name,&tag);
							if(tag=='i')
							{
								facePtr = (Interface *)facePtr;
								ptrZone=(Zone *)FindZonePtr(name,CgnsBase);
								Tag=strcmp(((CgnsBase->CgnsZone+numZones)->SolidFluidName),(ptrZone->SolidFluidName));
								if(Tag!=0 )/*Done By Ashish */
								{
									
									FindFaceCenterBC(surfPoints,centNeighbour);
								}
								else
								{
									sizePtrNeighbor=ptrZone->Size;
									indices1[0]=i+1,indices1[1]=j+2,indices1[2]=k+1;
									TransformCellIndices1(indices1,facePtr,indices2);
									AdjoiningZonePointerXYZ(ptrZone,indices2,ptrPoint);
									ptrNeighbour1X=ptrPoint[0];
									ptrNeighbour1Y=ptrPoint[1];
									ptrNeighbour1Z=ptrPoint[2];
		
									FindCellCentre(sizePtrNeighbor,ptrNeighbour1X,ptrNeighbour1Y,ptrNeighbour1Z,centNeighbour);
								}
							}
							else
							if(tag=='b')
								FindFaceCenterBC(surfPoints,centNeighbour);
						}
						else
						{
							ptrNeighbour1X=ptrPoint1X+NI;
							ptrNeighbour1Y=ptrPoint1Y+NI;
							ptrNeighbour1Z=ptrPoint1Z+NI;
							FindCellCentre(ptrSize,ptrNeighbour1X,ptrNeighbour1Y,ptrNeighbour1Z,centNeighbour);
						}
						X[3]=centNeighbour[0];
						Y[3]=centNeighbour[1];
						Z[3]=centNeighbour[2];
						break;
					case 4:/*Kmin->Bottam Face*/
						temp[0][0]=NI,temp[0][1]=1,temp[0][2]=NI*NJ,temp[1][0]=j,temp[1][1]=i,temp[1][2]=k;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						tag = 'a';//a default value
						if(k==0)
						{
							facePtr=(void *)FindInterfaceBCPtrMod((CgnsBase->CgnsZone+numZones),6,name,&tag);
							if(tag=='i')
							{
								facePtr = (Interface *)facePtr;
								ptrZone=(Zone *)FindZonePtr(name,CgnsBase);
								Tag=strcmp(((CgnsBase->CgnsZone+numZones)->SolidFluidName),(ptrZone->SolidFluidName));
								if(Tag!=0 )/*Done By Ashish */
								{
									
									FindFaceCenterBC(surfPoints,centNeighbour);
								}
								else
								{
									sizePtrNeighbor=ptrZone->Size;
									indices1[0]=i+1,indices1[1]=j+1,indices1[2]=k;
									TransformCellIndices1(indices1,facePtr,indices2);
									AdjoiningZonePointerXYZ(ptrZone,indices2,ptrPoint);
									ptrNeighbour1X=ptrPoint[0];
									ptrNeighbour1Y=ptrPoint[1];
									ptrNeighbour1Z=ptrPoint[2];
		
									FindCellCentre(sizePtrNeighbor,ptrNeighbour1X,ptrNeighbour1Y,ptrNeighbour1Z,centNeighbour);
								}
							}
							else
							if(tag=='b')
								FindFaceCenterBC(surfPoints,centNeighbour);
						}
						else
						{
							ptrNeighbour1X=ptrPoint1X-NI*NJ;
							ptrNeighbour1Y=ptrPoint1Y-NI*NJ;
							ptrNeighbour1Z=ptrPoint1Z-NI*NJ;
							FindCellCentre(ptrSize,ptrNeighbour1X,ptrNeighbour1Y,ptrNeighbour1Z,centNeighbour);
						}
						X[4]=centNeighbour[0];
						Y[4]=centNeighbour[1];
						Z[4]=centNeighbour[2];
						break;
				
					case 5:/*Kmax->Top face*/
						temp[0][0]=NI,temp[0][1]=1,temp[0][2]=NI*NJ,temp[1][0]=j,temp[1][1]=i,temp[1][2]=k+1;
						GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
						tag = 'a';//a default value
						if(k==NK-2)
						{
							facePtr=(void *)FindInterfaceBCPtrMod((CgnsBase->CgnsZone+numZones),5,name,&tag);
							if(tag=='i')
							{
								facePtr = (Interface *)facePtr;
								ptrZone=(Zone *)FindZonePtr(name,CgnsBase);
								Tag=strcmp(((CgnsBase->CgnsZone+numZones)->SolidFluidName),(ptrZone->SolidFluidName));
								if(Tag!=0 )/*Done By Ashish */
								{
									
									FindFaceCenterBC(surfPoints,centNeighbour);
								}
								else
								{
									sizePtrNeighbor=ptrZone->Size;
									indices1[0]=i+1,indices1[1]=j+1,indices1[2]=k+2;
									TransformCellIndices1(indices1,facePtr,indices2);
									AdjoiningZonePointerXYZ(ptrZone,indices2,ptrPoint);
									ptrNeighbour1X=ptrPoint[0];
									ptrNeighbour1Y=ptrPoint[1];
									ptrNeighbour1Z=ptrPoint[2];
		
									FindCellCentre(sizePtrNeighbor,ptrNeighbour1X,ptrNeighbour1Y,ptrNeighbour1Z,centNeighbour);
								}
							}
							else
							if(tag=='b')
								FindFaceCenterBC(surfPoints,centNeighbour);
						}
						else
						{
							ptrNeighbour1X=ptrPoint1X+NI*NJ;
							ptrNeighbour1Y=ptrPoint1Y+NI*NJ;
							ptrNeighbour1Z=ptrPoint1Z+NI*NJ;
							FindCellCentre(ptrSize,ptrNeighbour1X,ptrNeighbour1Y,ptrNeighbour1Z,centNeighbour);
						}
						X[5]=centNeighbour[0];
						Y[5]=centNeighbour[1];
						Z[5]=centNeighbour[2];
						break;
				}/*#######################End of Switch statement#############################33*/
				
				CalSurf(areacof,&surfPoints[0][0],surfVector);
				CheckOutward(surfVector,surfPoints,centroid);
				
				FindDeltaX(centroid,centNeighbour,delta_X);
				
				if(tag=='b' || (tag=='i' && Tag!=0))
				{
					ChangeFaceCenter(face,centroid,centNeighbour,delta_X,surfVector);
					
					switch(face)
					{
						case 0:
						case 1:
							slabIndex = (NJ-1)*k+j;
							break;
						case 2:
						case 3:
							slabIndex = (NI-1)*k+i;
							break;
						case 4:
						case 5:
							slabIndex = (NI-1)*j+i;
							break;
					}
					if(tag == 'b')//To be stored in conjugate interfaces also <<< Narsimha 8th March, 2008
					{
						bcPtr = (ZoneBoundary *)facePtr;
						*(bcPtr->FaceCenter+slabIndex*3+0) = centNeighbour[0];
						*(bcPtr->FaceCenter+slabIndex*3+1) = centNeighbour[1];
						*(bcPtr->FaceCenter+slabIndex*3+2) = centNeighbour[2];
					}
				}
				
				X[face]=centNeighbour[0];
				Y[face]=centNeighbour[1];
				Z[face]=centNeighbour[2];
	
				if(tag=='b' && (CgnsBase->flagTurb==1 || CgnsBase->flagWf==1))
				{
					bcPtr = (ZoneBoundary *) facePtr;
					
					if((CgnsBase->CgnsZone+numZones)->flagSolidFluid == 0 && bcPtr->Boundary_Tag == 22)
					{
						if(bcPtr->Direction == 1 || bcPtr->Direction == 2)
						  {
						*(bcPtr->deltaP+(NJ-1)*k+j) = sqrt(delta_X[0]*delta_X[0] + delta_X[1]*delta_X[1] + delta_X[2]*delta_X[2]);
						printf("direction = %d deltaP = %lf\n",bcPtr->Direction,*(bcPtr->deltaP+(NJ-1)*k+j));
						  }
						
						if(bcPtr->Direction == 3 || bcPtr->Direction == 4)
						  {
						*(bcPtr->deltaP+(NI-1)*k+i) = sqrt(delta_X[0]*delta_X[0] + delta_X[1]*delta_X[1] + delta_X[2]*delta_X[2]);
						//	printf("direction = %d deltaP = %lf\n",bcPtr->Direction,*(bcPtr->deltaP+(NI-1)*k+i));
						  }
						
						if(bcPtr->Direction == 5 || bcPtr->Direction == 6)
						  {
						*(bcPtr->deltaP+(NI-1)*j+i) = sqrt(delta_X[0]*delta_X[0] + delta_X[1]*delta_X[1] + delta_X[2]*delta_X[2]);
						//	printf("direction = %d deltaP = %lf\n",bcPtr->Direction,*(bcPtr->deltaP+(NI-1)*j+i));
						  }
					}
				}
				
				if(tag=='i' && ptrZone->flagSolidFluid == 1 && (CgnsBase->CgnsZone+numZones)->flagSolidFluid == 0 && (CgnsBase->flagTurb==1 || CgnsBase->flagWf==1))
				{
					IfacePtr = (Interface *) facePtr;
					
					if(IfacePtr->Direction == 1 || IfacePtr->Direction == 2)
						*(IfacePtr->deltaP+(NJ-1)*k+j) = sqrt(delta_X[0]*delta_X[0] + delta_X[1]*delta_X[1] + delta_X[2]*delta_X[2]);
					
					if(IfacePtr->Direction == 3 || IfacePtr->Direction == 4)
						*(IfacePtr->deltaP+(NI-1)*k+i) = sqrt(delta_X[0]*delta_X[0] + delta_X[1]*delta_X[1] + delta_X[2]*delta_X[2]);
					
					if(IfacePtr->Direction == 5 || IfacePtr->Direction == 6)
						*(IfacePtr->deltaP+(NI-1)*j+i) = sqrt(delta_X[0]*delta_X[0] + delta_X[1]*delta_X[1] + delta_X[2]*delta_X[2]);
				}
				
				for(comp=0;comp<3;comp++)
				{
					*(((CgnsBase->CgnsZone+numZones)->ZGeom)->deltaX+((NI-1)*(NJ-1)*k+(NI-1)*j+i)*6*3+3*face+comp)=delta_X[comp];
				} 
					
				if(CgnsBase->BaseParameters->flagRadiation==1||CgnsBase->BaseParameters->flagflowRadiation==1)
				{
					index[0]=i;index[1]=j;index[2]=k;
					CalculateDmF(CgnsBase,(CgnsBase->CgnsZone+numZones),index,face,surfVector);
					CalculateNormalVector(CgnsBase,(CgnsBase->CgnsZone+numZones),index,face,surfVector);
				}
				
				FindBeta(face,centroid,centNeighbour,surfPoints,surfVector,beta);
				
				ModifyBeta(face,beta,tag,Tag);
			
				for(comp=0;comp<3;comp++)
				{
					if(isnan(beta[comp])!=0)
					{
						printf("COMP=%d\tFACE=%d\tINDEX i=%d\tj=%d\tk=%d\tZone number=%d\n",comp,face,i,j,k,numZones+1);
						printf("Nan Encountered in Beta Caluclations\nExiting....\n");
						exit(1);
					}
					*(((CgnsBase->CgnsZone+numZones)->ZGeom)->Beta+((NI-1)*(NJ-1)*k+(NI-1)*j+i)*6*3+3*face+comp)=beta[comp]; 
				}
				for(comp=0;comp<3;comp++)
					*(((CgnsBase->CgnsZone+numZones)->ZGeom)->SurfaceArea+((NI-1)*(NJ-1)*k+(NI-1)*j+i)*6*3+3*face+comp)=*(surfVector+comp);
					
// 				fprintf(r,"face=%d\tbeta[0]:%lf\tbeta[1]:%lf\tbeta[2]:%lf\t\n",face,beta[0],beta[1],beta[2]);
// 				fprintf(r,"face=%d\tSarea[0]:%lf\tSarea[1]:%lf\tSarea[2]:%lf\t\n",face,*(surfVector+0),*(surfVector+1),*(surfVector+2));
			}//################End of for looping over faces###################
			
			//Calculating volumes
			principalDiag[0]=*(coordX+NI*NJ*k+NI*(j+1)+i)-*(coordX+NI*NJ*(k+1)+NI*j+(i+1));
			principalDiag[1]=*(coordY+NI*NJ*k+NI*(j+1)+i)-*(coordY+NI*NJ*(k+1)+NI*j+(i+1));
			principalDiag[2]=*(coordZ+NI*NJ*k+NI*(j+1)+i)-*(coordZ+NI*NJ*(k+1)+NI*j+(i+1));
			
			for(comp=0;comp<3;comp++)
			{
				adjacentFace[0][comp]=*(((CgnsBase->CgnsZone+numZones)->ZGeom)->SurfaceArea+((NI-1)*(NJ-1)*k+(NI-1)*j+i)*6*3+3*0+comp);
				adjacentFace[1][comp]=*(((CgnsBase->CgnsZone+numZones)->ZGeom)->SurfaceArea+((NI-1)*(NJ-1)*k+(NI-1)*j+i)*6*3+3*3+comp);
				adjacentFace[2][comp]=*(((CgnsBase->CgnsZone+numZones)->ZGeom)->SurfaceArea+((NI-1)*(NJ-1)*k+(NI-1)*j+i)*6*3+3*4+comp);
			}
			for(comp=0;comp<3;comp++)
			{ 
				areaSum[comp]=adjacentFace[0][comp]+adjacentFace[1][comp]+adjacentFace[2][comp];
// 				fprintf(r,"i=%d\t,j=%d\t,k=%d\t,areaSum[%d]:%lf\t\n",i,j,k,comp,areaSum[comp]);
			}
			
			volume=0.0;
			for(comp=0;comp<3;comp++)
				volume=volume+areaSum[comp]*principalDiag[comp];
			
			volume=fabs(volume/3.0);
// 			fprintf(r,"i=%d\t,j=%d\t,k=%d\t,volume:%12.8e\t\n",i,j,k,volume);
			*(((CgnsBase->CgnsZone+numZones)->ZGeom)->Volume+(NI-1)*(NJ-1)*k+(NI-1)*j+i)=volume;
			
			if(CgnsBase->BaseParameters->flag2d == 0)//3D Case
			{
				del[0]=X[1]-X[0];del[1]=Y[1]-Y[0];del[2]=Z[1]-Z[0];
				del[3]=X[3]-X[2];del[4]=Y[3]-Y[2];del[5]=Z[3]-Z[2];
				del[6]=X[5]-X[4];del[7]=Y[5]-Y[4];del[8]=Z[5]-Z[4];
				FindDeltaInverse3D(del,(((CgnsBase->CgnsZone+numZones)->ZGeom)->DeltaInverse+((NI-1)*(NJ-1)*k+(NI-1)*j+i)*9));
			}
			else
			{
				del[0]=X[1]-X[0];del[1]=Y[1]-Y[0];
				del[2]=X[3]-X[2];del[3]=Y[3]-Y[2];
				FindDeltaInverse2D(del,(((CgnsBase->CgnsZone+numZones)->ZGeom)->DeltaInverse+((NI-1)*(NJ-1)*k+(NI-1)*j+i)*4));
			}
			
		}    /*End of for looping over NI,NJ,NK*/ 
		
		if(CgnsBase->BaseParameters->flagRadiation==1||CgnsBase->BaseParameters->flagflowRadiation==1)
		{
			CalculateOmega(CgnsBase,((CgnsBase->CgnsZone+numZones)->ZGeom)->Omega);
	
			CalculateVector(CgnsBase, ((CgnsBase->CgnsZone+numZones)->ZGeom)->VectorX, ((CgnsBase->CgnsZone+numZones)->ZGeom)->VectorY, ((CgnsBase->CgnsZone+numZones)->ZGeom)->VectorZ);
		} 
		
	
	}/* ################End of for loop for nZones##############################*/
	
// 	fclose(r);
	return 0;
}

int FindDeltaX(double center[3],double centerN[3],double delta_X[3])
{
	delta_X[0]=fabs(centerN[0]-center[0]);
	delta_X[1]=fabs(centerN[1]-center[1]);
	delta_X[2]=fabs(centerN[2]-center[2]);
	return 0;
}
/**
This function FindDeltaInverse2D finds the inverse of a 2D matrix.
*/
int FindDeltaInverse2D(double delta[4],double *addressDeltaInverse)
{
	double dx1,dx2,dy1,dy2;
	double det;
	
	dx1=delta[0];
	dy1=delta[1];
	dx2=delta[2];
	dy2=delta[3];
	
	det=(dx1*dy2-dy1*dx2);
	if(det==0)
		printf("\n..........Singularity Encountered in Matrix inversion..........\n");

	*(addressDeltaInverse+0) = dy2/det;
	*(addressDeltaInverse+1) = -dy1/det; 
	*(addressDeltaInverse+2) = -dx2/det; 
	*(addressDeltaInverse+3) = dx1/det; 
	return 0;
}

/**
This function FindDeltaInverse3D finds the inverse of a 3D matrix.
*/
int FindDeltaInverse3D(double delta[9],double *addressDeltaInverse)
{
	double dx1,dx2,dx3,dy1,dy2,dy3,dz1,dz2,dz3;
	double det;
	
	dx1=delta[0];
	dy1=delta[1];
	dz1=delta[2];
	
	dx2=delta[3];
	dy2=delta[4];
	dz2=delta[5];
	
	dx3=delta[6];	
	dy3=delta[7];
	dz3=delta[8];
	
	det=dx1*(dy2*dz3-dy3*dz2)-dy1*(dx2*dz3-dx3*dz2)+dz1*(dx2*dy3-dx3*dy2);
	
	if(det==0)
		printf("\n..........Singularity Encountered in Matrix inversion..........\n");
	
	*(addressDeltaInverse+0) =  (dy2*dz3-dy3*dz2)/det;
	*(addressDeltaInverse+1) = -(dy1*dz3-dy3*dz1)/det;
	*(addressDeltaInverse+2) =  (dy1*dz2-dy2*dz1)/det;
	
	*(addressDeltaInverse+3) = -(dx2*dz3-dx3*dz2)/det;
	*(addressDeltaInverse+4) =  (dx1*dz3-dx3*dz1)/det;
	*(addressDeltaInverse+5) = -(dx1*dz2-dx2*dz1)/det;
	
	*(addressDeltaInverse+6) =  (dx2*dy3-dx3*dy2)/det;
	*(addressDeltaInverse+7) = -(dx1*dy3-dx3*dy1)/det;
	*(addressDeltaInverse+8) =  (dx1*dy2-dx2*dy1)/det;
	
	return 0;
}

int ChangeFaceCenter(int face, double center[3],double centerN[3],double delta_X[3],double surfVector[3])
{
	int comp;
	double unitOutwardNormal[3],SurfArea,dotProduct;
	
	SurfArea = sqrt(surfVector[0] * surfVector[0] + surfVector[1] * surfVector[1] + surfVector[2] * surfVector[2]);
	if(fabs(SurfArea)!=0.0)
	{
		unitOutwardNormal[0] = surfVector[0] / SurfArea;
		unitOutwardNormal[1] = surfVector[1] / SurfArea;
		unitOutwardNormal[2] = surfVector[2] / SurfArea;
		if((isnan(unitOutwardNormal[0])!=0) || (isnan(unitOutwardNormal[1])!=0) || (isnan(unitOutwardNormal[2])!=0))
		{
			printf("Nan encountered in surface area calculation\n");
			printf("Surface Vector is %e\t%e\t%e\t%e\n",surfVector[0],surfVector[1],surfVector[2],SurfArea);
			//exit(12);
		}
		//Computing "delta_X" which is in the increasing INDEX direction
		switch(face)
		{
			case 0:
			case 2:
			case 4:
				delta_X[0] = (center[0]-centerN[0]);
				delta_X[1] = (center[1]-centerN[1]);
				delta_X[2] = (center[2]-centerN[2]);
				break;
			case 1:
			case 3:
			case 5:
				delta_X[0] = (centerN[0]-center[0]);
				delta_X[1] = (centerN[1]-center[1]);
				delta_X[2] = (centerN[2]-center[2]);
				break;
		}
		
		dotProduct = delta_X[0] * unitOutwardNormal[0] + delta_X[1] * unitOutwardNormal[1] + delta_X[2] * unitOutwardNormal[2];
		
		//Calculating coordinates of the point where the surface normal passing through the cell centroid intersects the face
		for(comp=0; comp <3; comp++)
			centerN[comp] = center[comp] + fabs(dotProduct) * unitOutwardNormal[comp];
		
		//Recomputing new "delta_X" which is in the increasing INDEX direction
		switch(face)
		{
			case 0:
			case 2:
			case 4:
				delta_X[0] = (center[0]-centerN[0]);
				delta_X[1] = (center[1]-centerN[1]);
				delta_X[2] = (center[2]-centerN[2]);
				break;
			case 1:
			case 3:
			case 5:
				delta_X[0] = (centerN[0]-center[0]);
				delta_X[1] = (centerN[1]-center[1]);
				delta_X[2] = (centerN[2]-center[2]);
				break;
		}
	}
	return 0;
}

int ModifyBeta(int face,double beta[3],char tag,int Tag)
{
	if(tag=='b' || (tag=='i' && Tag!=0))
	{
		//printf("Before modification ... face=%d beta[0]=%e beta[1]=%e beta[2]=%e\n", face,beta[0],beta[1],beta[2]);
		switch(face)
		{
			case 0:
			case 1:
				beta[1]=0.0;
				beta[2]=0.0;
				break;
			case 2:
			case 3:
				beta[0]=0.0;
				beta[2]=0.0;
				break;
			case 4:
			case 5:
				beta[0]=0.0;
				beta[1]=0.0;
				break;
		}
		//printf("After modification  ... face=%d beta[0]=%e beta[1]=%e beta[2]=%e\n", face,beta[0],beta[1],beta[2]);
	}
	
	return 0;
}

