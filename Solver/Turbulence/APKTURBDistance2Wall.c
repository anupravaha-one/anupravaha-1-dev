#include "APKTURBDistance2Wall.h"
int Distance2Wall(Base *base)
{
	Zone *zone;
	int i,j,k,n=0,cntZone,numZones,NI,NJ,NK;
	int indexP,*ptrSize,temp[2][3],slabCellIndex;
	int countBC,numBC,numBCCells=0,count=0,countIface,numIface;
	double *ptrPoint1X,*ptrPoint1Y,*ptrPoint1Z,*coordX,*coordY,*coordZ;
	double surfPoints[4][3],BCFaceCenter[3],centroid[3],x,y,z,xb,yb,zb,miny;
	ZoneBoundary *zBc;
	Interface *Iface;
	//FILE *fp;
	
	struct FaceCenter
	{
		double X,Y,Z;
	};
	struct FaceCenter *FCptr;
	//fp = fopen("SOLUTIONS/miny.dat", "w");

	numZones=base->Number_of_zones;

	//Calculating the number of Wall Boundary Cells in the whole domain
	for(cntZone=0; cntZone < numZones; cntZone++)
	{
		zone=(base->CgnsZone+cntZone);

		ptrSize= zone->Size;

		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);

		if(zone->flagSolidFluid == 0)	//If the zone is Fluid
		{
			numBC = zone->Number_of_BC_Patches;
			for(countBC=0; countBC < numBC; countBC++)
			{
				zBc = (zone->ZoneBc+countBC);
				if(zBc->Boundary_Tag == 22)
				{
					switch(zBc->Direction)
					{
						case 1://Imax face
						numBCCells = numBCCells + (NJ-1)*(NK-1);
						break;
						case 2://Imin face
						numBCCells = numBCCells + (NJ-1)*(NK-1);
						break;
						case 3://Jmax face
						numBCCells = numBCCells + (NI-1)*(NK-1);
						break;
						case 4://Imin face
						numBCCells = numBCCells + (NI-1)*(NK-1);
						break;
						case 5://Kmax face
						numBCCells = numBCCells + (NI-1)*(NJ-1);
						break;
						case 6://Imin face
						numBCCells = numBCCells + (NI-1)*(NJ-1);
						break;
						default:
						printf("*** Improper Boundary direction in the function Distance2Wall() ***\n");
						break;
					}
				}
			}
		}
		else	//If the boundary is SOLID ---> Compulsory in case of Conjugate Heat Transfer
		{
			numBCCells = numBCCells + 2 * ((NI-1)*(NJ-1) + (NJ-1)*(NK-1) + (NI-1)*(NK-1));
		}
	}

	printf("---------------- No. of Wall Boundary Cells = %d -----------------\n", numBCCells);
	FCptr = (struct FaceCenter *)malloc(numBCCells * sizeof(struct FaceCenter));

	//Storing the X,Y,Z coordinates of all the Wall Boundary Cells in a series of sturcture pointed by "FCptr"
	for(cntZone=0; cntZone < numZones; cntZone++)
	{
		zone=(base->CgnsZone+cntZone);

		ptrSize= zone->Size;

		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		coordX = zone->CgnsGrid->X;
		coordY = zone->CgnsGrid->Y; 
		coordZ = zone->CgnsGrid->Z;
		
		if(zone->flagSolidFluid == 0)	//If the zone is Fluid
		{
			numBC = zone->Number_of_BC_Patches;
			for(countBC=0; countBC < numBC; countBC++)
			{
				zBc = (zone->ZoneBc+countBC);
				if(zBc->Boundary_Tag == 22)
				{
					switch(zBc->Direction)
					{
						case 1://Imax face
						
							for(k=0; k<NK-1; k++)
							for(j=0; j<NJ-1; j++)
							{
								temp[0][0]=NI*NJ,temp[0][1]=NI,temp[0][2]=1,temp[1][0]=k,temp[1][1]=j,temp[1][2]=NI-1;  
								GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
								FindFaceCenterBC(surfPoints,BCFaceCenter);
								(FCptr+n)->X = BCFaceCenter[0];
								(FCptr+n)->Y = BCFaceCenter[1];
								(FCptr+n)->Z = BCFaceCenter[2];
								n++;
							}
						break;
						case 2://Imin face
							for(k=0; k<NK-1; k++)
							for(j=0; j<NJ-1; j++)
							{
								temp[0][0]=NI*NJ,temp[0][1]=NI,temp[0][2]=1,temp[1][0]=k,temp[1][1]=j,temp[1][2]=0;  
								GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
								FindFaceCenterBC(surfPoints,BCFaceCenter);
								(FCptr+n)->X = BCFaceCenter[0];
								(FCptr+n)->Y = BCFaceCenter[1];
								(FCptr+n)->Z = BCFaceCenter[2];
								n++;
							}
						break;
						case 3://Jmax face
							for(k=0; k<NK-1; k++)
							for(i=0; i<NI-1; i++)
							{
								temp[0][0]=NI*NJ,temp[0][1]=1,temp[0][2]=NI,temp[1][0]=k,temp[1][1]=i,temp[1][2]=NJ-1;
								GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
								FindFaceCenterBC(surfPoints,BCFaceCenter);
								(FCptr+n)->X = BCFaceCenter[0];
								(FCptr+n)->Y = BCFaceCenter[1];
								(FCptr+n)->Z = BCFaceCenter[2];
								n++;
							}
						break;
						case 4://Jmin face
							for(k=0; k<NK-1; k++)
							for(i=0; i<NI-1; i++)
							{
								temp[0][0]=NI*NJ,temp[0][1]=1,temp[0][2]=NI,temp[1][0]=k,temp[1][1]=i,temp[1][2]=0;
								GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
								FindFaceCenterBC(surfPoints,BCFaceCenter);
								(FCptr+n)->X = BCFaceCenter[0];
								(FCptr+n)->Y = BCFaceCenter[1];
								(FCptr+n)->Z = BCFaceCenter[2];
								n++;
							}
						break;
						case 5://Kmax face
							for(j=0; j<NJ-1; j++)
							for(i=0; i<NI-1; i++)
							{
								temp[0][0]=NI,temp[0][1]=1,temp[0][2]=NI*NJ,temp[1][0]=j,temp[1][1]=i,temp[1][2]=NK-1;
								GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
								FindFaceCenterBC(surfPoints,BCFaceCenter);
								(FCptr+n)->X = BCFaceCenter[0];
								(FCptr+n)->Y = BCFaceCenter[1];
								(FCptr+n)->Z = BCFaceCenter[2];
								n++;
							}
						break;
						case 6://Kmin face
							for(j=0; j<NJ-1; j++)
							for(i=0; i<NI-1; i++)
							{
								temp[0][0]=NI,temp[0][1]=1,temp[0][2]=NI*NJ,temp[1][0]=j,temp[1][1]=i,temp[1][2]=0;
								GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
								FindFaceCenterBC(surfPoints,BCFaceCenter);
								(FCptr+n)->X = BCFaceCenter[0];
								(FCptr+n)->Y = BCFaceCenter[1];
								(FCptr+n)->Z = BCFaceCenter[2];
								n++;
							}
						break;
						default:
						printf("*** Improper Boundary direction in the function Distance2Wall() ***\n");
						break;
					}
				}
			}
		}
		else	//If the boundary is SOLID ---> Compulsory in case of Conjugate Heat Transfer
		{
			//Looping over all the six faces of the Zone
			for(countBC=1; countBC <= 6; countBC++)
			{
				switch(countBC)
				{
					case 1://Imax face
						for(k=0; k<NK-1; k++)
						for(j=0; j<NJ-1; j++)
						{
							temp[0][0]=NI*NJ,temp[0][1]=NI,temp[0][2]=1,temp[1][0]=k,temp[1][1]=j,temp[1][2]=NI-1;  
							GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
							FindFaceCenterBC(surfPoints,BCFaceCenter);
							(FCptr+n)->X = BCFaceCenter[0];
							(FCptr+n)->Y = BCFaceCenter[1];
							(FCptr+n)->Z = BCFaceCenter[2];
							n++;
						}
					break;
					case 2://Imin face
						for(k=0; k<NK-1; k++)
						for(j=0; j<NJ-1; j++)
						{
							temp[0][0]=NI*NJ,temp[0][1]=NI,temp[0][2]=1,temp[1][0]=k,temp[1][1]=j,temp[1][2]=0;  
							GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
							FindFaceCenterBC(surfPoints,BCFaceCenter);
							(FCptr+n)->X = BCFaceCenter[0];
							(FCptr+n)->Y = BCFaceCenter[1];
							(FCptr+n)->Z = BCFaceCenter[2];
							n++;
						}
					break;
					case 3://Jmax face
						for(k=0; k<NK-1; k++)
						for(i=0; i<NI-1; i++)
						{
							temp[0][0]=NI*NJ,temp[0][1]=1,temp[0][2]=NI,temp[1][0]=k,temp[1][1]=i,temp[1][2]=NJ-1;
							GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
							FindFaceCenterBC(surfPoints,BCFaceCenter);
							(FCptr+n)->X = BCFaceCenter[0];
							(FCptr+n)->Y = BCFaceCenter[1];
							(FCptr+n)->Z = BCFaceCenter[2];
							n++;
						}
					break;
					case 4://Jmin face
						for(k=0; k<NK-1; k++)
						for(i=0; i<NI-1; i++)
						{
							temp[0][0]=NI*NJ,temp[0][1]=1,temp[0][2]=NI,temp[1][0]=k,temp[1][1]=i,temp[1][2]=0;
							GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
							FindFaceCenterBC(surfPoints,BCFaceCenter);
							(FCptr+n)->X = BCFaceCenter[0];
							(FCptr+n)->Y = BCFaceCenter[1];
							(FCptr+n)->Z = BCFaceCenter[2];
							n++;
						}
					break;
					case 5://Kmax face
						for(j=0; j<NJ-1; j++)
						for(i=0; i<NI-1; i++)
						{
							temp[0][0]=NI,temp[0][1]=1,temp[0][2]=NI*NJ,temp[1][0]=j,temp[1][1]=i,temp[1][2]=NK-1;
							GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
							FindFaceCenterBC(surfPoints,BCFaceCenter);
							(FCptr+n)->X = BCFaceCenter[0];
							(FCptr+n)->Y = BCFaceCenter[1];
							(FCptr+n)->Z = BCFaceCenter[2];
							n++;
						}
					break;
					case 6://Kmin face
						for(j=0; j<NJ-1; j++)
						for(i=0; i<NI-1; i++)
						{
							temp[0][0]=NI,temp[0][1]=1,temp[0][2]=NI*NJ,temp[1][0]=j,temp[1][1]=i,temp[1][2]=0;
							GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
							FindFaceCenterBC(surfPoints,BCFaceCenter);
							(FCptr+n)->X = BCFaceCenter[0];
							(FCptr+n)->Y = BCFaceCenter[1];
							(FCptr+n)->Z = BCFaceCenter[2];
							n++;
						}
					break;
					default:
					printf("*** Improper Boundary direction in the function Distance2Wall() ***\n");
					break;
				}
			}
		}
	}
	
	//Computing the minimum distance from each cell to any arbitrary noslip solid wall
	if(numBCCells != 0)
	{
		for(cntZone=0; cntZone < numZones; cntZone++)
		{
			zone=(base->CgnsZone+cntZone);
	
			ptrSize= zone->Size;
	
			NI=*(ptrSize+0);
			NJ=*(ptrSize+1);
			NK=*(ptrSize+2);
			coordX = zone->CgnsGrid->X;
			coordY = zone->CgnsGrid->Y; 
			coordZ = zone->CgnsGrid->Z;
			
			if(zone->flagSolidFluid == 0)	//If the zone is Fluid
			{
				for(k=0;k<(NK-1);k++)
				for(j=0;j<(NJ-1);j++)
				{
					for(i=0;i<(NI-1);i++)
					{
						indexP = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
						ptrPoint1X=coordX+NI*NJ*k+NI*j+i;
						ptrPoint1Y=coordY+NI*NJ*k+NI*j+i;
						ptrPoint1Z=coordZ+NI*NJ*k+NI*j+i;//point1 is the point with minimum indices in a cell i.e (i,j,k)
						FindCellCentre(ptrSize,ptrPoint1X,ptrPoint1Y,ptrPoint1Z,centroid);
						x = centroid[0];
						y = centroid[1];
						z = centroid[2];
						
						xb = (FCptr+0)->X;
						yb = (FCptr+0)->Y;
						zb = (FCptr+0)->Z;
						miny = sqrt((x-xb)*(x-xb) + (y-yb)*(y-yb) + (z-zb)*(z-zb));
						
						for(count=1; count < numBCCells; count++)
						{
							xb = (FCptr+count)->X;
							yb = (FCptr+count)->Y;
							zb = (FCptr+count)->Z;
							
							miny = min(miny , sqrt((x-xb)*(x-xb) + (y-yb)*(y-yb) + (z-zb)*(z-zb)));
						}
						
						*(zone->ZGeom->minWallDistance+indexP) = miny;
						//fprintf(fp,"%e\t", *(zone->ZGeom->minWallDistance+indexP));
					}
					//fprintf(fp,"\n");
				}
				//fprintf(fp,"\n\n");
				
				//Calculating the minimum distance to any Solid Wall for the Boundary Slabs
				numBC = zone->Number_of_BC_Patches;
				for(countBC=0; countBC < numBC; countBC++)
				{
					zBc = (zone->ZoneBc+countBC);
					if(zBc->ZSlab->minWallDistance != NULL)
					{
						switch(zBc->Direction)
						{
							case 1://Imax face
							
								for(k=0; k<NK-1; k++)
								for(j=0; j<NJ-1; j++)
								{
									temp[0][0]=NI*NJ,temp[0][1]=NI,temp[0][2]=1,temp[1][0]=k,temp[1][1]=j,temp[1][2]=NI-1;  
									GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
									FindFaceCenterBC(surfPoints,BCFaceCenter);
									x = BCFaceCenter[0];
									y = BCFaceCenter[1];
									z = BCFaceCenter[2];
	
									xb = (FCptr+0)->X;
									yb = (FCptr+0)->Y;
									zb = (FCptr+0)->Z;
									miny = sqrt((x-xb)*(x-xb) + (y-yb)*(y-yb) + (z-zb)*(z-zb));
									
									for(count=1; count < numBCCells; count++)
									{
										xb = (FCptr+count)->X;
										yb = (FCptr+count)->Y;
										zb = (FCptr+count)->Z;
										
										miny = min(miny , sqrt((x-xb)*(x-xb) + (y-yb)*(y-yb) + (z-zb)*(z-zb)));
									}
									
									slabCellIndex = (NJ-1)*k+j;
									*(zBc->ZSlab->minWallDistance+slabCellIndex) = miny;
								}
							break;
							case 2://Imin face
								for(k=0; k<NK-1; k++)
								for(j=0; j<NJ-1; j++)
								{
									temp[0][0]=NI*NJ,temp[0][1]=NI,temp[0][2]=1,temp[1][0]=k,temp[1][1]=j,temp[1][2]=0;  
									GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
									FindFaceCenterBC(surfPoints,BCFaceCenter);
									x = BCFaceCenter[0];
									y = BCFaceCenter[1];
									z = BCFaceCenter[2];
									
									xb = (FCptr+0)->X;
									yb = (FCptr+0)->Y;
									zb = (FCptr+0)->Z;
									miny = sqrt((x-xb)*(x-xb) + (y-yb)*(y-yb) + (z-zb)*(z-zb));
									
									for(count=1; count < numBCCells; count++)
									{
										xb = (FCptr+count)->X;
										yb = (FCptr+count)->Y;
										zb = (FCptr+count)->Z;
										
										miny = min(miny , sqrt((x-xb)*(x-xb) + (y-yb)*(y-yb) + (z-zb)*(z-zb)));
									}
									
									slabCellIndex = (NJ-1)*k+j;
									*(zBc->ZSlab->minWallDistance+slabCellIndex) = miny;
								}
							break;
							case 3://Jmax face
								for(k=0; k<NK-1; k++)
								for(i=0; i<NI-1; i++)
								{
									temp[0][0]=NI*NJ,temp[0][1]=1,temp[0][2]=NI,temp[1][0]=k,temp[1][1]=i,temp[1][2]=NJ-1;
									GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
									FindFaceCenterBC(surfPoints,BCFaceCenter);
									x = BCFaceCenter[0];
									y = BCFaceCenter[1];
									z = BCFaceCenter[2];
									
									xb = (FCptr+0)->X;
									yb = (FCptr+0)->Y;
									zb = (FCptr+0)->Z;
									miny = sqrt((x-xb)*(x-xb) + (y-yb)*(y-yb) + (z-zb)*(z-zb));
									
									for(count=1; count < numBCCells; count++)
									{
										xb = (FCptr+count)->X;
										yb = (FCptr+count)->Y;
										zb = (FCptr+count)->Z;
										
										miny = min(miny , sqrt((x-xb)*(x-xb) + (y-yb)*(y-yb) + (z-zb)*(z-zb)));
									}
									
									slabCellIndex = (NI-1)*k+i;
									*(zBc->ZSlab->minWallDistance+slabCellIndex) = miny;
								}
							break;
							case 4://Jmin face
								for(k=0; k<NK-1; k++)
								for(i=0; i<NI-1; i++)
								{
									temp[0][0]=NI*NJ,temp[0][1]=1,temp[0][2]=NI,temp[1][0]=k,temp[1][1]=i,temp[1][2]=0;
									GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
									FindFaceCenterBC(surfPoints,BCFaceCenter);
									x = BCFaceCenter[0];
									y = BCFaceCenter[1];
									z = BCFaceCenter[2];
									
									xb = (FCptr+0)->X;
									yb = (FCptr+0)->Y;
									zb = (FCptr+0)->Z;
									miny = sqrt((x-xb)*(x-xb) + (y-yb)*(y-yb) + (z-zb)*(z-zb));
									
									for(count=1; count < numBCCells; count++)
									{
										xb = (FCptr+count)->X;
										yb = (FCptr+count)->Y;
										zb = (FCptr+count)->Z;
										
										miny = min(miny , sqrt((x-xb)*(x-xb) + (y-yb)*(y-yb) + (z-zb)*(z-zb)));
									}
									
									slabCellIndex = (NI-1)*k+i;
									*(zBc->ZSlab->minWallDistance+slabCellIndex) = miny;
								}
							break;
							case 5://Kmax face
								for(j=0; j<NJ-1; j++)
								for(i=0; i<NI-1; i++)
								{
									temp[0][0]=NI,temp[0][1]=1,temp[0][2]=NI*NJ,temp[1][0]=j,temp[1][1]=i,temp[1][2]=NK-1;
									GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
									FindFaceCenterBC(surfPoints,BCFaceCenter);
									x = BCFaceCenter[0];
									y = BCFaceCenter[1];
									z = BCFaceCenter[2];
									
									xb = (FCptr+0)->X;
									yb = (FCptr+0)->Y;
									zb = (FCptr+0)->Z;
									miny = sqrt((x-xb)*(x-xb) + (y-yb)*(y-yb) + (z-zb)*(z-zb));
									
									for(count=1; count < numBCCells; count++)
									{
										xb = (FCptr+count)->X;
										yb = (FCptr+count)->Y;
										zb = (FCptr+count)->Z;
										
										miny = min(miny , sqrt((x-xb)*(x-xb) + (y-yb)*(y-yb) + (z-zb)*(z-zb)));
									}
									
									slabCellIndex = (NI-1)*j+i;
									*(zBc->ZSlab->minWallDistance+slabCellIndex) = miny;
								}
							break;
							case 6://Kmin face
								for(j=0; j<NJ-1; j++)
								for(i=0; i<NI-1; i++)
								{
									temp[0][0]=NI,temp[0][1]=1,temp[0][2]=NI*NJ,temp[1][0]=j,temp[1][1]=i,temp[1][2]=0;
									GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
									FindFaceCenterBC(surfPoints,BCFaceCenter);
									x = BCFaceCenter[0];
									y = BCFaceCenter[1];
									z = BCFaceCenter[2];
									
									xb = (FCptr+0)->X;
									yb = (FCptr+0)->Y;
									zb = (FCptr+0)->Z;
									miny = sqrt((x-xb)*(x-xb) + (y-yb)*(y-yb) + (z-zb)*(z-zb));
									
									for(count=1; count < numBCCells; count++)
									{
										xb = (FCptr+count)->X;
										yb = (FCptr+count)->Y;
										zb = (FCptr+count)->Z;
										
										miny = min(miny , sqrt((x-xb)*(x-xb) + (y-yb)*(y-yb) + (z-zb)*(z-zb)));
									}
									
									slabCellIndex = (NI-1)*j+i;
									*(zBc->ZSlab->minWallDistance+slabCellIndex) = miny;
								}
							break;
							default:
							printf("*** Improper Boundary direction in the function Distance2Wall() ***\n");
							break;
						}
					}
				}
				
				//Calculating the minimum distance to any Solid Wall for the Interface Slabs
				numIface = zone->Number_of_Interfaces;
				for(countIface=0; countIface < numIface; countIface++)
				{
					Iface = (zone->Iface+countIface);
					if(Iface->ZSlab->minWallDistance != NULL)
					{
						switch(Iface->Direction)
						{
							case 1://Imax face
							
								for(k=0; k<NK-1; k++)
								for(j=0; j<NJ-1; j++)
								{
									temp[0][0]=NI*NJ,temp[0][1]=NI,temp[0][2]=1,temp[1][0]=k,temp[1][1]=j,temp[1][2]=NI-1;  
									GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
									FindFaceCenterBC(surfPoints,BCFaceCenter);
									x = BCFaceCenter[0];
									y = BCFaceCenter[1];
									z = BCFaceCenter[2];
									
									xb = (FCptr+0)->X;
									yb = (FCptr+0)->Y;
									zb = (FCptr+0)->Z;
									miny = sqrt((x-xb)*(x-xb) + (y-yb)*(y-yb) + (z-zb)*(z-zb));
									
									for(count=1; count < numBCCells; count++)
									{
										xb = (FCptr+count)->X;
										yb = (FCptr+count)->Y;
										zb = (FCptr+count)->Z;
										
										miny = min(miny , sqrt((x-xb)*(x-xb) + (y-yb)*(y-yb) + (z-zb)*(z-zb)));
									}
									
									slabCellIndex = (NJ-1)*k+j;
									*(Iface->ZSlab->minWallDistance+slabCellIndex) = miny;
								}
							break;
							case 2://Imin face
								for(k=0; k<NK-1; k++)
								for(j=0; j<NJ-1; j++)
								{
									temp[0][0]=NI*NJ,temp[0][1]=NI,temp[0][2]=1,temp[1][0]=k,temp[1][1]=j,temp[1][2]=0;  
									GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
									FindFaceCenterBC(surfPoints,BCFaceCenter);
									x = BCFaceCenter[0];
									y = BCFaceCenter[1];
									z = BCFaceCenter[2];
									
									xb = (FCptr+0)->X;
									yb = (FCptr+0)->Y;
									zb = (FCptr+0)->Z;
									miny = sqrt((x-xb)*(x-xb) + (y-yb)*(y-yb) + (z-zb)*(z-zb));
									
									for(count=1; count < numBCCells; count++)
									{
										xb = (FCptr+count)->X;
										yb = (FCptr+count)->Y;
										zb = (FCptr+count)->Z;
										
										miny = min(miny , sqrt((x-xb)*(x-xb) + (y-yb)*(y-yb) + (z-zb)*(z-zb)));
									}
									
									slabCellIndex = (NJ-1)*k+j;
									*(Iface->ZSlab->minWallDistance+slabCellIndex) = miny;
								}
							break;
							case 3://Jmax face
								for(k=0; k<NK-1; k++)
								for(i=0; i<NI-1; i++)
								{
									temp[0][0]=NI*NJ,temp[0][1]=1,temp[0][2]=NI,temp[1][0]=k,temp[1][1]=i,temp[1][2]=NJ-1;
									GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
									FindFaceCenterBC(surfPoints,BCFaceCenter);
									x = BCFaceCenter[0];
									y = BCFaceCenter[1];
									z = BCFaceCenter[2];
									
									xb = (FCptr+0)->X;
									yb = (FCptr+0)->Y;
									zb = (FCptr+0)->Z;
									miny = sqrt((x-xb)*(x-xb) + (y-yb)*(y-yb) + (z-zb)*(z-zb));
									
									for(count=1; count < numBCCells; count++)
									{
										xb = (FCptr+count)->X;
										yb = (FCptr+count)->Y;
										zb = (FCptr+count)->Z;
										
										miny = min(miny , sqrt((x-xb)*(x-xb) + (y-yb)*(y-yb) + (z-zb)*(z-zb)));
									}
									
									slabCellIndex = (NI-1)*k+i;
									*(Iface->ZSlab->minWallDistance+slabCellIndex) = miny;
								}
							break;
							case 4://Jmin face
								for(k=0; k<NK-1; k++)
								for(i=0; i<NI-1; i++)
								{
									temp[0][0]=NI*NJ,temp[0][1]=1,temp[0][2]=NI,temp[1][0]=k,temp[1][1]=i,temp[1][2]=0;
									GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
									FindFaceCenterBC(surfPoints,BCFaceCenter);
									x = BCFaceCenter[0];
									y = BCFaceCenter[1];
									z = BCFaceCenter[2];
									
									xb = (FCptr+0)->X;
									yb = (FCptr+0)->Y;
									zb = (FCptr+0)->Z;
									miny = sqrt((x-xb)*(x-xb) + (y-yb)*(y-yb) + (z-zb)*(z-zb));
									
									for(count=1; count < numBCCells; count++)
									{
										xb = (FCptr+count)->X;
										yb = (FCptr+count)->Y;
										zb = (FCptr+count)->Z;
										
										miny = min(miny , sqrt((x-xb)*(x-xb) + (y-yb)*(y-yb) + (z-zb)*(z-zb)));
									}
									
									slabCellIndex = (NI-1)*k+i;
									*(Iface->ZSlab->minWallDistance+slabCellIndex) = miny;
								}
							break;
							case 5://Kmax face
								for(j=0; j<NJ-1; j++)
								for(i=0; i<NI-1; i++)
								{
									temp[0][0]=NI,temp[0][1]=1,temp[0][2]=NI*NJ,temp[1][0]=j,temp[1][1]=i,temp[1][2]=NK-1;
									GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
									FindFaceCenterBC(surfPoints,BCFaceCenter);
									x = BCFaceCenter[0];
									y = BCFaceCenter[1];
									z = BCFaceCenter[2];
									
									xb = (FCptr+0)->X;
									yb = (FCptr+0)->Y;
									zb = (FCptr+0)->Z;
									miny = sqrt((x-xb)*(x-xb) + (y-yb)*(y-yb) + (z-zb)*(z-zb));
									
									for(count=1; count < numBCCells; count++)
									{
										xb = (FCptr+count)->X;
										yb = (FCptr+count)->Y;
										zb = (FCptr+count)->Z;
										
										miny = min(miny , sqrt((x-xb)*(x-xb) + (y-yb)*(y-yb) + (z-zb)*(z-zb)));
									}
									
									slabCellIndex = (NI-1)*j+i;
									*(Iface->ZSlab->minWallDistance+slabCellIndex) = miny;
								}
							break;
							case 6://Kmin face
								for(j=0; j<NJ-1; j++)
								for(i=0; i<NI-1; i++)
								{
									temp[0][0]=NI,temp[0][1]=1,temp[0][2]=NI*NJ,temp[1][0]=j,temp[1][1]=i,temp[1][2]=0;
									GetSurfacePoints(coordX,coordY,coordZ,surfPoints,temp);
									FindFaceCenterBC(surfPoints,BCFaceCenter);
									x = BCFaceCenter[0];
									y = BCFaceCenter[1];
									z = BCFaceCenter[2];
									
									xb = (FCptr+0)->X;
									yb = (FCptr+0)->Y;
									zb = (FCptr+0)->Z;
									miny = sqrt((x-xb)*(x-xb) + (y-yb)*(y-yb) + (z-zb)*(z-zb));
									
									for(count=1; count < numBCCells; count++)
									{
										xb = (FCptr+count)->X;
										yb = (FCptr+count)->Y;
										zb = (FCptr+count)->Z;
										
										miny = min(miny , sqrt((x-xb)*(x-xb) + (y-yb)*(y-yb) + (z-zb)*(z-zb)));
									}
									
									slabCellIndex = (NI-1)*j+i;
									*(Iface->ZSlab->minWallDistance+slabCellIndex) = miny;
								}
							break;
							default:
							printf("*** Improper Boundary direction in the function Distance2Wall() ***\n");
							break;
						}
					}
				}
			}
		}
	}
	else
	{
		printf("** NOTE: There is no SOLID wall in the geometry\n");
		printf(" ... Ideally making the minimum distance to wall for all the points in the domain equal to 9999999999.0\n");
		for(cntZone=0; cntZone < numZones; cntZone++)
		{
			zone=(base->CgnsZone+cntZone);
	
			ptrSize= zone->Size;
	
			NI=*(ptrSize+0);
			NJ=*(ptrSize+1);
			NK=*(ptrSize+2);
			
			miny = 9999999999.0;
			if(zone->flagSolidFluid == 0)	//If the zone is Fluid
			{
				//Looping over all the interior cells
				for(k=0;k<(NK-1);k++)
				for(j=0;j<(NJ-1);j++)
				for(i=0;i<(NI-1);i++)
				{
					indexP = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
					*(zone->ZGeom->minWallDistance+indexP) = miny;
				}
				
				//Looping over all the Boundary Slabs
				numBC = zone->Number_of_BC_Patches;
				for(countBC=0; countBC < numBC; countBC++)
				{
					zBc = (zone->ZoneBc+countBC);
					if(zBc->ZSlab->minWallDistance != NULL)
					{
						switch(zBc->Direction)
						{
							case 1://Imax face
							case 2://Imin face
								for(k=0; k<NK-1; k++)
								for(j=0; j<NJ-1; j++)
								{
									slabCellIndex = (NJ-1)*k+j;
									*(zBc->ZSlab->minWallDistance+slabCellIndex) = miny;
								}
							break;
							case 3://Jmax face
							case 4://Jmin face
								for(k=0; k<NK-1; k++)
								for(i=0; i<NI-1; i++)
								{
									slabCellIndex = (NI-1)*k+i;
									*(zBc->ZSlab->minWallDistance+slabCellIndex) = miny;
								}
							break;
							case 5://Kmax face
							case 6://Kmin face
								for(j=0; j<NJ-1; j++)
								for(i=0; i<NI-1; i++)
								{
									slabCellIndex = (NI-1)*j+i;
									*(zBc->ZSlab->minWallDistance+slabCellIndex) = miny;
								}
							break;
							default:
							printf("*** Improper Boundary direction in the function Distance2Wall() ***\n");
							break;
						}
					}
				}
				
				//Looping over all the Interface Slabs
				numIface = zone->Number_of_Interfaces;
				for(countIface=0; countIface < numIface; countIface++)
				{
					Iface = (zone->Iface+countIface);
					if(Iface->ZSlab->minWallDistance != NULL)
					{
						switch(Iface->Direction)
						{
							case 1://Imax face
							case 2://Imin face
								for(k=0; k<NK-1; k++)
								for(j=0; j<NJ-1; j++)
								{
									slabCellIndex = (NJ-1)*k+j;
									*(Iface->ZSlab->minWallDistance+slabCellIndex) = miny;
								}
							break;
							case 3://Jmax face
							case 4://Jmin face
								for(k=0; k<NK-1; k++)
								for(i=0; i<NI-1; i++)
								{
									slabCellIndex = (NI-1)*k+i;
									*(Iface->ZSlab->minWallDistance+slabCellIndex) = miny;
								}
							break;
							case 5://Kmax face
							case 6://Kmin face
								for(j=0; j<NJ-1; j++)
								for(i=0; i<NI-1; i++)
								{
									slabCellIndex = (NI-1)*j+i;
									*(Iface->ZSlab->minWallDistance+slabCellIndex) = miny;
								}
							break;
							default:
							printf("*** Improper Boundary direction in the function Distance2Wall() ***\n");
							break;
						}
					}
				}
			}
		}
	}

	free(FCptr);
	//fclose(fp);
	return 0;
}

