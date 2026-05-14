#include "APKQuickParticleTracking.h" 

int QuickParticleTracking(Base *base,double XP,double YP,double ZP)
{
	int i,j,k,test=1;
	//double *ptrPoint1X,*ptrPoint1Y,*ptrPoint1Z;
	double ptrPointEastNorth[3],ptrPointWestSouth[3];//ptrPointEastNorthY,ptrPointEastNorthZ;
	double *coordX,*coordY,*coordZ;
	double surfVector[3],testVector[3];
	int NI,NJ,NK,indices1[3],indices2[3];
	double checkDotProduct;
	char tag,name[33];
	Zone *neighbourZone,*desiredZone;
	void *facePtr;
//	int izone,nZones,Tag=1,count=0; 

/***************** SECONDARY SEARCH ALGORITHM, which locates the exact zone and the i,j,k cell index ********/
	
	desiredZone=((base->BaseParameters)->ParticlePosition->zone);
	i=((base->BaseParameters)->ParticlePosition->i);
	j=((base->BaseParameters)->ParticlePosition->j);
	k=((base->BaseParameters)->ParticlePosition->k);
	
	do
	{
		coordX=(desiredZone->CgnsGrid)->X;
		coordY=(desiredZone->CgnsGrid)->Y; 
		coordZ=(desiredZone->CgnsGrid)->Z;
		
		NI=*((desiredZone)->Size+0);
		NJ=*((desiredZone)->Size+1);
		NK=*((desiredZone)->Size+2);
		
		ptrPointWestSouth[0] = *(coordX+NI*NJ*k+NI*j+i); // Test point for West and South
		ptrPointWestSouth[1] = *(coordY+NI*NJ*k+NI*j+i);
		ptrPointWestSouth[2] = *(coordZ+NI*NJ*k+NI*j+i);
		
		ptrPointEastNorth[0] = *(coordX+NI*NJ*k+NI*(j+1)+(i+1)); // Test point for East and North
		ptrPointEastNorth[1] = *(coordY+NI*NJ*k+NI*(j+1)+(i+1));
		ptrPointEastNorth[2] = *(coordZ+NI*NJ*k+NI*(j+1)+(i+1));
		
		testVector[0] = XP - ptrPointEastNorth[0];
		testVector[1] = YP - ptrPointEastNorth[1];
		testVector[2] = ZP - ptrPointEastNorth[2];
		
		// East side Surface Area Vectors //
		surfVector[0] = *(((desiredZone->ZGeom)->SurfaceArea)+((NI-1)*(NJ-1)*k+(NI-1)*j+i)*6*3+3*1+0);
		surfVector[1] = *(((desiredZone->ZGeom)->SurfaceArea)+((NI-1)*(NJ-1)*k+(NI-1)*j+i)*6*3+3*1+1);
		surfVector[2] = *(((desiredZone->ZGeom)->SurfaceArea)+((NI-1)*(NJ-1)*k+(NI-1)*j+i)*6*3+3*1+2);
		
		checkDotProduct = testVector[0]*surfVector[0] + testVector[1]*surfVector[1] + testVector[2]*surfVector[2];
		
		//printf("testVectors are %12.8lf %12.8lf %12.8lf\t SurfVectors are %12.8lf %12.8lf %12.8lf\n",testVector[0],testVector[1],testVector[2],surfVector[0],surfVector[1],surfVector[3]);
		//printf("checkDotProduct1 = %12.8lf\n",checkDotProduct);
		//exit(1);
		if(checkDotProduct>0.0) // Implies East is satisfied
		{
			// North side Surface Area Vectors //
			surfVector[0] = *(((desiredZone->ZGeom)->SurfaceArea)+((NI-1)*(NJ-1)*k+(NI-1)*j+i)*6*3+3*3+0);
			surfVector[1] = *(((desiredZone->ZGeom)->SurfaceArea)+((NI-1)*(NJ-1)*k+(NI-1)*j+i)*6*3+3*3+1);
			surfVector[2] = *(((desiredZone->ZGeom)->SurfaceArea)+((NI-1)*(NJ-1)*k+(NI-1)*j+i)*6*3+3*3+2);
			
			checkDotProduct = testVector[0]*surfVector[0] + testVector[1]*surfVector[1] + testVector[2]*surfVector[2];
			
			//printf("checkDotProduct2 = %12.8lf\n",checkDotProduct);
			if(checkDotProduct>0.0) // Implies North is satisfied
			{
				if(i==NI-2 && j==NJ-2)
				{
					facePtr=(void *)FindInterfaceBCPtrMod(desiredZone,1,name,&tag);
					if(tag=='i')
					{	// Implies moving only East side
						facePtr = (Interface *)facePtr;
						neighbourZone=(Zone *)FindZonePtr(name,base);
						
						indices1[0]=i+2,indices1[1]=j+1,indices1[2]=k+1;
						TransformCellIndices1(indices1,facePtr,indices2);

						i = indices2[0]-1; j = indices2[1]-1; k = indices2[2]-1;
						desiredZone = neighbourZone;
					}
					else if(tag=='b')
					{
						facePtr=(void *)FindInterfaceBCPtrMod(desiredZone,3,name,&tag);
						if(tag=='i')
						{	// Implies moving only North side
							facePtr = (Interface *)facePtr;
							neighbourZone=(Zone *)FindZonePtr(name,base);

							indices1[0]=i+1,indices1[1]=j+2,indices1[2]=k+1;
							TransformCellIndices1(indices1,facePtr,indices2);

							i = indices2[0]-1; j = indices2[1]-1; k = indices2[2]-1;
							desiredZone = neighbourZone;
						}
						else if(tag=='b')
						{
							printf("Situation where the corner cell has both sides boundaries\n");
							printf("Algorithm not designed to handle such situations\n");
							printf("Exiting from ParticleTracking.c\n");
							exit(1);
						}
					}
				}
				else if(i==NI-2 && j!=NJ-2)
				{
					facePtr=(void *)FindInterfaceBCPtrMod(desiredZone,1,name,&tag);
					if(tag=='i')
					{	// Implies moving only NorthEast side
						facePtr = (Interface *)facePtr;
						neighbourZone=(Zone *)FindZonePtr(name,base);
						
						indices1[0]=i+2,indices1[1]=j+2,indices1[2]=k+1;
						TransformCellIndices1(indices1,facePtr,indices2);

						i = indices2[0]-1; j = indices2[1]-1; k = indices2[2]-1;
						desiredZone = neighbourZone;
					}
					else if(tag=='b')
					{	// Implies moving only North side
						i = i;
						j = j+1;
					}
				}
				else if(i!=NI-2 && j==NJ-2)
				{
					facePtr=(void *)FindInterfaceBCPtrMod(desiredZone,3,name,&tag);
					if(tag=='i')
					{	// Implies moving only NorthEast side
						facePtr = (Interface *)facePtr;
						neighbourZone=(Zone *)FindZonePtr(name,base);
						
						indices1[0]=i+2,indices1[1]=j+2,indices1[2]=k+1;
						TransformCellIndices1(indices1,facePtr,indices2);

						i = indices2[0]-1; j = indices2[1]-1; k = indices2[2]-1;
						desiredZone = neighbourZone;
					}
					else if(tag=='b')
					{	// Implies moving only East side
						i = i+1;
						j = j;
					}
				}
				else
				{	// Implies moving only NorthEast side
					i = i+1;
					j = j+1;
				}
			}
			else// Mostly South
			{
				testVector[0] = XP - ptrPointWestSouth[0];
				testVector[1] = YP - ptrPointWestSouth[1];
				testVector[2] = ZP - ptrPointWestSouth[2];
				
				// South side Surface Area Vectors //
				surfVector[0] = *(((desiredZone->ZGeom)->SurfaceArea)+((NI-1)*(NJ-1)*k+(NI-1)*j+i)*6*3+3*2+0);
				surfVector[1] = *(((desiredZone->ZGeom)->SurfaceArea)+((NI-1)*(NJ-1)*k+(NI-1)*j+i)*6*3+3*2+1);
				surfVector[2] = *(((desiredZone->ZGeom)->SurfaceArea)+((NI-1)*(NJ-1)*k+(NI-1)*j+i)*6*3+3*2+2);
				
				checkDotProduct = testVector[0]*surfVector[0] + testVector[1]*surfVector[1] + testVector[2]*surfVector[2];
				//printf("checkDotProduct3 = %12.8lf\n",checkDotProduct);
				if(checkDotProduct>0.0) // Implies South Satisfied.
				{
					if(i==NI-2 && j==0)
					{
						facePtr=(void *)FindInterfaceBCPtrMod(desiredZone,1,name,&tag);
						if(tag=='i')
						{	// Implies moving only East side
							facePtr = (Interface *)facePtr;
							neighbourZone=(Zone *)FindZonePtr(name,base);

							indices1[0]=i+2,indices1[1]=j+1,indices1[2]=k+1;
							TransformCellIndices1(indices1,facePtr,indices2);

							i = indices2[0]-1; j = indices2[1]-1; k = indices2[2]-1;
							desiredZone = neighbourZone;
						}
						else if(tag=='b')
						{
							facePtr=(void *)FindInterfaceBCPtrMod(desiredZone,4,name,&tag);
							if(tag=='i')
							{	// Implies moving only South side
								facePtr = (Interface *)facePtr;
								neighbourZone=(Zone *)FindZonePtr(name,base);

								indices1[0]=i+1,indices1[1]=j,indices1[2]=k+1;
								TransformCellIndices1(indices1,facePtr,indices2);

								i = indices2[0]-1; j = indices2[1]-1; k = indices2[2]-1;
								desiredZone = neighbourZone;
							}
							else if(tag=='b')
							{
								printf("Situation where the corner cell has both sides boundaries\n");
								printf("Algorithm not designed to handle such situations\n");
								printf("Exiting from ParticleTracking.c\n");
								exit(1);
							}
						}
					}
					else if(i==NI-2 && j!=0)
					{
						facePtr=(void *)FindInterfaceBCPtrMod(desiredZone,1,name,&tag);
						if(tag=='i')
						{	// Implies moving only SouthEast side
							facePtr = (Interface *)facePtr;
							neighbourZone=(Zone *)FindZonePtr(name,base);

							indices1[0]=i+2,indices1[1]=j,indices1[2]=k+1;
							TransformCellIndices1(indices1,facePtr,indices2);

							i = indices2[0]-1; j = indices2[1]-1; k = indices2[2]-1;
							desiredZone = neighbourZone;
						}
						else if(tag=='b')
						{	// Implies moving only South side
							i = i;
							j = j-1;
						}
					}
					else if(i!=NI-2 && j==0)
					{
						facePtr=(void *)FindInterfaceBCPtrMod(desiredZone,4,name,&tag);
						if(tag=='i')
						{	// Implies moving only SouthEast side
							facePtr = (Interface *)facePtr;
							neighbourZone=(Zone *)FindZonePtr(name,base);

							indices1[0]=i+2,indices1[1]=j,indices1[2]=k+1;
							TransformCellIndices1(indices1,facePtr,indices2);

							i = indices2[0]-1; j = indices2[1]-1; k = indices2[2]-1;
							desiredZone = neighbourZone;
						}
						else if(tag=='b')
						{	// Implies moving only East side
							i = i+1;
							j = j;
						}
					}
					else
					{
						i = i+1;
						j = j-1;
					}
				}
				else // Only East
				{
					if(i==NI-2)
					{
						facePtr=(void *)FindInterfaceBCPtrMod(desiredZone,1,name,&tag);
						if(tag=='i')
						{	// Implies moving only East side
							facePtr = (Interface *)facePtr;
							neighbourZone=(Zone *)FindZonePtr(name,base);

							indices1[0]=i+2,indices1[1]=j+1,indices1[2]=k+1;
							TransformCellIndices1(indices1,facePtr,indices2);

							i = indices2[0]-1; j = indices2[1]-1; k = indices2[2]-1;
							desiredZone = neighbourZone;
						}
						else if(tag=='b')
						{
							printf("Only East case and East is a boundary. Can't move anywhere hence exiting from ParticleTracking\n");
							exit(1);
						}
					}
					else
					{
						i = i+1;
						j = j;
					}
				}
			}
		}
		else // Implies main if loop for East is not Satisfied going for West side.
		{
			testVector[0] = XP - ptrPointWestSouth[0];
			testVector[1] = YP - ptrPointWestSouth[1];
			testVector[2] = ZP - ptrPointWestSouth[2];

			// West side Surface Area Vectors //
			surfVector[0] = *(((desiredZone->ZGeom)->SurfaceArea)+((NI-1)*(NJ-1)*k+(NI-1)*j+i)*6*3+3*0+0);
			surfVector[1] = *(((desiredZone->ZGeom)->SurfaceArea)+((NI-1)*(NJ-1)*k+(NI-1)*j+i)*6*3+3*0+1);
			surfVector[2] = *(((desiredZone->ZGeom)->SurfaceArea)+((NI-1)*(NJ-1)*k+(NI-1)*j+i)*6*3+3*0+2);

			checkDotProduct = testVector[0]*surfVector[0] + testVector[1]*surfVector[1] + testVector[2]*surfVector[2];
			//printf("checkDotProduct4 = %12.8lf\n",checkDotProduct);
			if(checkDotProduct>0.0) // West side is satisfied
			{
				// North side Surface Area Vectors //
				surfVector[0] = *(((desiredZone->ZGeom)->SurfaceArea)+((NI-1)*(NJ-1)*k+(NI-1)*j+i)*6*3+3*3+0);
				surfVector[1] = *(((desiredZone->ZGeom)->SurfaceArea)+((NI-1)*(NJ-1)*k+(NI-1)*j+i)*6*3+3*3+1);
				surfVector[2] = *(((desiredZone->ZGeom)->SurfaceArea)+((NI-1)*(NJ-1)*k+(NI-1)*j+i)*6*3+3*3+2);

				checkDotProduct = testVector[0]*surfVector[0] + testVector[1]*surfVector[1] + testVector[2]*surfVector[2];
				//printf("checkDotProduct5 = %12.8lf\n",checkDotProduct);
				if(checkDotProduct>0.0) // North is Satisfied
				{
					if(i==0 && j==NJ-2)
					{
						facePtr=(void *)FindInterfaceBCPtrMod(desiredZone,2,name,&tag);
						if(tag=='i')
						{	// Implies moving only West side
							facePtr = (Interface *)facePtr;
							neighbourZone=(Zone *)FindZonePtr(name,base);

							indices1[0]=i,indices1[1]=j+1,indices1[2]=k+1;
							TransformCellIndices1(indices1,facePtr,indices2);

							i = indices2[0]-1; j = indices2[1]-1; k = indices2[2]-1;
							desiredZone = neighbourZone;
						}
						else if(tag=='b')
						{
							facePtr=(void *)FindInterfaceBCPtrMod(desiredZone,3,name,&tag);
							if(tag=='i')
							{	// Implies moving only North side
								facePtr = (Interface *)facePtr;
								neighbourZone=(Zone *)FindZonePtr(name,base);

								indices1[0]=i+1,indices1[1]=j+2,indices1[2]=k+1;
								TransformCellIndices1(indices1,facePtr,indices2);

								i = indices2[0]-1; j = indices2[1]-1; k = indices2[2]-1;
								desiredZone = neighbourZone;
							}
							else if(tag=='b')
							{
								printf("Situation where the corner cell has both sides boundaries\n");
								printf("Algorithm not designed to handle such situations\n");
								printf("Exiting from ParticleTracking.c\n");
								exit(1);
							}
						}
					}
					else if(i==0 && j!=NJ-2)
					{
						facePtr=(void *)FindInterfaceBCPtrMod(desiredZone,2,name,&tag);
						if(tag=='i')
						{	// Implies moving only NorthWest side
							facePtr = (Interface *)facePtr;
							neighbourZone=(Zone *)FindZonePtr(name,base);

							indices1[0]=i,indices1[1]=j+2,indices1[2]=k+1;
							TransformCellIndices1(indices1,facePtr,indices2);

							i = indices2[0]-1; j = indices2[1]-1; k = indices2[2]-1;
							desiredZone = neighbourZone;
						}
						else if(tag=='b')
						{	// Implies moving only North side
							i = i;
							j = j+1;
						}
					}
					else if(i!=0 && j==NJ-2)
					{
						facePtr=(void *)FindInterfaceBCPtrMod(desiredZone,3,name,&tag);
						if(tag=='i')
						{	// Implies moving only NorthWest side
							facePtr = (Interface *)facePtr;
							neighbourZone=(Zone *)FindZonePtr(name,base);

							indices1[0]=i,indices1[1]=j+2,indices1[2]=k+1;
							TransformCellIndices1(indices1,facePtr,indices2);

							i = indices2[0]-1; j = indices2[1]-1; k = indices2[2]-1;
							desiredZone = neighbourZone;
						}
						else if(tag=='b')
						{	// Implies moving only West side
							i = i-1;
							j = j;
						}
					}
					else
					{	// Implies moving only NorthWest side
						i = i-1;
						j = j+1;
					}
				}
				else	// South Side
				{
					testVector[0] = XP - ptrPointWestSouth[0];
					testVector[1] = YP - ptrPointWestSouth[1];
					testVector[2] = ZP - ptrPointWestSouth[2];

					// South side Surface Area Vectors //
					surfVector[0] = *(((desiredZone->ZGeom)->SurfaceArea)+((NI-1)*(NJ-1)*k+(NI-1)*j+i)*6*3+3*2+0);
					surfVector[1] = *(((desiredZone->ZGeom)->SurfaceArea)+((NI-1)*(NJ-1)*k+(NI-1)*j+i)*6*3+3*2+1);
					surfVector[2] = *(((desiredZone->ZGeom)->SurfaceArea)+((NI-1)*(NJ-1)*k+(NI-1)*j+i)*6*3+3*2+2);

					checkDotProduct = testVector[0]*surfVector[0] + testVector[1]*surfVector[1] + testVector[2]*surfVector[2];
					//printf("checkDotProduct6 = %12.8lf\n",checkDotProduct);
					if(checkDotProduct>0.0) // Implies South Satisfied.
					{
						if(i==0 && j==0)
						{
							facePtr=(void *)FindInterfaceBCPtrMod(desiredZone,2,name,&tag);
							if(tag=='i')
							{	// Implies moving only West side
								facePtr = (Interface *)facePtr;
								neighbourZone=(Zone *)FindZonePtr(name,base);

								indices1[0]=i,indices1[1]=j+1,indices1[2]=k+1;
								TransformCellIndices1(indices1,facePtr,indices2);

								i = indices2[0]-1; j = indices2[1]-1; k = indices2[2]-1;
								desiredZone = neighbourZone;
							}
							else if(tag=='b')
							{
								facePtr=(void *)FindInterfaceBCPtrMod(desiredZone,4,name,&tag);
								if(tag=='i')
								{	// Implies moving only South side
									facePtr = (Interface *)facePtr;
									neighbourZone=(Zone *)FindZonePtr(name,base);

									indices1[0]=i+1,indices1[1]=j,indices1[2]=k+1;
									TransformCellIndices1(indices1,facePtr,indices2);

									i = indices2[0]-1; j = indices2[1]-1; k = indices2[2]-1;
									desiredZone = neighbourZone;
								}
								else if(tag=='b')
								{
									printf("Situation where the corner cell has both sides boundaries\n");
									printf("Algorithm not designed to handle such situations\n");
									printf("Exiting from ParticleTracking.c\n");
									exit(1);
								}
							}
						}
						else if(i==0 && j!=0)
						{
							facePtr=(void *)FindInterfaceBCPtrMod(desiredZone,2,name,&tag);
							if(tag=='i')
							{	// Implies moving only SouthWest side
								facePtr = (Interface *)facePtr;
								neighbourZone=(Zone *)FindZonePtr(name,base);

								indices1[0]=i,indices1[1]=j,indices1[2]=k+1;
								TransformCellIndices1(indices1,facePtr,indices2);

								i = indices2[0]-1; j = indices2[1]-1; k = indices2[2]-1;
								desiredZone = neighbourZone;
							}
							else if(tag=='b')
							{	// Implies moving only South side
								i = i;
								j = j-1;
							}
						}
						else if(i!=0 && j==0)
						{
							facePtr=(void *)FindInterfaceBCPtrMod(desiredZone,4,name,&tag);
							if(tag=='i')
							{	// Implies moving only SouthWest side
								facePtr = (Interface *)facePtr;
								neighbourZone=(Zone *)FindZonePtr(name,base);

								indices1[0]=i,indices1[1]=j,indices1[2]=k+1;
								TransformCellIndices1(indices1,facePtr,indices2);

								i = indices2[0]-1; j = indices2[1]-1; k = indices2[2]-1;
								desiredZone = neighbourZone;
							}
							else if(tag=='b')
							{	// Implies moving only West side
								i = i-1;
								j = j;
							}
						}
						else
						{
							i = i-1;
							j = j-1;
						}
					}
					else // Only West
					{
						if(i==0)
						{
							facePtr=(void *)FindInterfaceBCPtrMod(desiredZone,2,name,&tag);
							if(tag=='i')
							{	// Implies moving only East side
								facePtr = (Interface *)facePtr;
								neighbourZone=(Zone *)FindZonePtr(name,base);

								indices1[0]=i,indices1[1]=j+1,indices1[2]=k+1;
								TransformCellIndices1(indices1,facePtr,indices2);

								i = indices2[0]-1; j = indices2[1]-1; k = indices2[2]-1;
								desiredZone = neighbourZone;
							}
							else if(tag=='b')
							{
								printf("Only West case and West is a boundary. Can't move anywhere hence exiting from ParticleTracking\n");
								exit(1);
							}
						}
						else
						{
							i = i-1;
							j = j;
						}
					}
				}
			} // End of If loop for West
			else // Implies Neither East nor West, but may be Only North or Only South
			{
				// Checking for North
				testVector[0] = XP - ptrPointEastNorth[0];
				testVector[1] = YP - ptrPointEastNorth[1];
				testVector[2] = ZP - ptrPointEastNorth[2];
				
				// North side Surface Area Vectors //
				surfVector[0] = *(((desiredZone->ZGeom)->SurfaceArea)+((NI-1)*(NJ-1)*k+(NI-1)*j+i)*6*3+3*3+0);
				surfVector[1] = *(((desiredZone->ZGeom)->SurfaceArea)+((NI-1)*(NJ-1)*k+(NI-1)*j+i)*6*3+3*3+1);
				surfVector[2] = *(((desiredZone->ZGeom)->SurfaceArea)+((NI-1)*(NJ-1)*k+(NI-1)*j+i)*6*3+3*3+2);

				checkDotProduct = testVector[0]*surfVector[0] + testVector[1]*surfVector[1] + testVector[2]*surfVector[2];
				if(checkDotProduct>0.0)
				{
					if(j==NJ-2)
					{
						facePtr=(void *)FindInterfaceBCPtrMod(desiredZone,3,name,&tag);
						if(tag=='i')
						{	// Implies moving only North side
							facePtr = (Interface *)facePtr;
							neighbourZone=(Zone *)FindZonePtr(name,base);

							indices1[0]=i+1,indices1[1]=j+2,indices1[2]=k+1;
							TransformCellIndices1(indices1,facePtr,indices2);

							i = indices2[0]-1; j = indices2[1]-1; k = indices2[2]-1;
							desiredZone = neighbourZone;
						}
						else if(tag=='b')
						{
							printf("Move Only North Situation but reached a north bundary\n");
							printf("Algorithm not designed to handle such situations\n");
							printf("Exiting from ParticleTracking.c\n");
							exit(1);
						}
					}
					else
					{
						i = i; j=j+1;
					}
				}
				else // Checking for South
				{
					testVector[0] = XP - ptrPointWestSouth[0];
					testVector[1] = YP - ptrPointWestSouth[1];
					testVector[2] = ZP - ptrPointWestSouth[2];
					
					surfVector[0] = *(((desiredZone->ZGeom)->SurfaceArea)+((NI-1)*(NJ-1)*k+(NI-1)*j+i)*6*3+3*2+0);
					surfVector[1] = *(((desiredZone->ZGeom)->SurfaceArea)+((NI-1)*(NJ-1)*k+(NI-1)*j+i)*6*3+3*2+1);
					surfVector[2] = *(((desiredZone->ZGeom)->SurfaceArea)+((NI-1)*(NJ-1)*k+(NI-1)*j+i)*6*3+3*2+2);

					checkDotProduct = testVector[0]*surfVector[0] + testVector[1]*surfVector[1] + testVector[2]*surfVector[2];
					
					if(checkDotProduct>0.0) //Implies South satisfied
					{
						if(j==0)
						{
							facePtr=(void *)FindInterfaceBCPtrMod(desiredZone,4,name,&tag);
							if(tag=='i')
							{	// Implies moving only South side
								facePtr = (Interface *)facePtr;
								neighbourZone=(Zone *)FindZonePtr(name,base);

								indices1[0]=i+1,indices1[1]=j,indices1[2]=k+1;
								TransformCellIndices1(indices1,facePtr,indices2);

								i = indices2[0]-1; j = indices2[1]-1; k = indices2[2]-1;
								desiredZone = neighbourZone;
							}
							else if(tag=='b')
							{
								printf("Move Only North Situation but reached a north bundary\n");
								printf("Algorithm not designed to handle such situations\n");
								printf("Exiting from ParticleTracking.c\n");
								exit(1);
							}
						}
						else
						{
							i = i; j = j-1;
						}
					}
					else
					{
						test = 0;
//						printf("\n********************* Required Cell Data ******************\n");
//						printf("Required i=%d j=%d zoneName=%s\n",i,j,desiredZone->zonename);
//						printf("\n************************************************************\n");
					}
				}
			}
		}
//	printf("Required i=%d j=%d k=%d currentZoneName=%s\n",i,j,k,desiredZone->zonename);
	//fprintf(fp1,"\n%s\n%12.8lf\t%12.8lf\t%12.8lf\t%d\t%d\t%d\n",zone->zonename,XP,YP,ZP,i,j,k);
	//printf("Test = %d\n",test);
	//Index[0]=i;Index[1]=j;Index[2]=k;
	}while(test==1);
	
	((base->BaseParameters)->ParticlePosition->i)=i;
	((base->BaseParameters)->ParticlePosition->j)=j;
	((base->BaseParameters)->ParticlePosition->k)=k;
	((base->BaseParameters)->ParticlePosition->zone)=desiredZone;
	
	//printf("\n%s\t%d\t%d\t%d\n",((base->BaseParameters)->ParticlePosition->zone->zonename),((base->BaseParameters)->ParticlePosition->i),((base->BaseParameters)->ParticlePosition->j),((base->BaseParameters)->ParticlePosition->k));
	//printf("\n%12.8lf\t%12.8lf\t%12.8lf\n",XP,YP,ZP);
	
return 0;
}
