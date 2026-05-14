#include "APKMultigridExtractVariables.h"

/**
The function ExtractVariables updates all the interface slab values in the current zone for the input 
variable, the variables here are Velocities, temperature, Scalars and so on. 
It takes variable tag, Base pointer and Zone pointer as input. It loops over all the faces of a Zone and finds 
Interface pointers, and then with the help of function FillerInterfaceVariables() it updates the interface 
slabs. It also checks whether the adjoining zone is also of the same type as the current zone, i.e; whether the adjoining
zone is  a solid or fluid if the current zone is a fluid etc, and calls FillerConjugateInterface() to update and implement
interface conditions if it is a conjugate heat transfer case. When it encounters boundary in the Zone then it calls
the function  update_boundary() or update_boundaryTurb(),according to whether it is a laminar or turbulent flow problem,
which implement the boundary conditions.
*/
int MultigridExtractVariables(int variableTag,Base *base,Zone *zone,int gdi)
{
	int iface,Tag=1;
	char tag,name[33];
	void *ptr;
	Interface *facePtr;
	ZoneBoundary *bcPtr;
	Zone *ptrZoneNeighbor;
	int i,j,k,NI,NJ,NK,c,*ptrSize;
	
	ptrSize = (zone->grid+gdi)->GSize;
	
	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
	
	for(iface=0;iface<6;iface++)
	{
		ptr = (void *)FindInterfaceBCPtrMod(zone,iface+1,name,&tag);
		if(tag=='i')
		{
			if(gdi == 0)
			{
				facePtr=(Interface *)ptr;
				
				ptrZoneNeighbor = (Zone *)FindZonePtr(name,base);
				
				Tag=strcmp(((zone)->SolidFluidName),(ptrZoneNeighbor->SolidFluidName));
				
				if(Tag == 0)
				{
					FillerInterfaceVariables(variableTag,facePtr,ptrZoneNeighbor,iface);
				}
			}
			else
			{
				switch(iface)
				{
					case 0:
						for(k=0;k<NK-1;k++)
						for(j=0;j<NJ-1;j++)
						{
							*((zone->grid+gdi)->SlabV[0]+(NJ-1)*k+j) = *((zone->grid+gdi)->DummySlabV[0]+(NJ-1)*k+j);
						}
					break;
					case 1:
						for(k=0;k<NK-1;k++)
						for(j=0;j<NJ-1;j++)
						{
							*((zone->grid+gdi)->SlabV[1]+(NJ-1)*k+j) = *((zone->grid+gdi)->DummySlabV[1]+(NJ-1)*k+j);
						}
					break;
					case 2:
						for(k=0;k<NK-1;k++)
						for(i=0;i<NI-1;i++)
						{
							*((zone->grid+gdi)->SlabV[2]+(NI-1)*k+i) = *((zone->grid+gdi)->DummySlabV[2]+(NI-1)*k+i);
						}
					break;
					case 3:
						for(k=0;k<NK-1;k++)
						for(i=0;i<NI-1;i++)
						{
							*((zone->grid+gdi)->SlabV[3]+(NI-1)*k+i) = *((zone->grid+gdi)->DummySlabV[3]+(NI-1)*k+i);
						}
					break;
					case 4:
						for(j=0;j<NJ-1;j++)
						for(i=0;i<NI-1;i++)
						{
							*((zone->grid+gdi)->SlabV[4]+(NI-1)*j+i) = *((zone->grid+gdi)->DummySlabV[4]+(NI-1)*j+i);
						}
					break;
					case 5:
						for(j=0;j<NJ-1;j++)
						for(i=0;i<NI-1;i++)
						{
							*((zone->grid+gdi)->SlabV[5]+(NI-1)*j+i) = *((zone->grid+gdi)->DummySlabV[5]+(NI-1)*j+i);
						}
					break;
				}
			}
		}
		else if(tag=='b')
		{
			bcPtr=(ZoneBoundary *)ptr;
			if(gdi == 0)
			{
				if(zone->Flowsol->flagTurb == 0)
					update_boundary(base,variableTag,zone,bcPtr);
			}
			else
			{
				c = (bcPtr->Data_Set+3)->Bc_Type;
				
				if(base->BaseParameters->flag2d == 1)
				{
					switch(c)
					{
						case 0://Homogeneous Neuman
							switch(iface)
							{
								case 0:
									for(k=0;k<NK-1;k++)
									for(j=0;j<NJ-1;j++)
									{
										*((zone->grid+gdi)->SlabV[0]+(NJ-1)*k+j) = 0.0;//*((zone->grid+gdi)->V+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2);
									}
								break;
								case 1:
									for(k=0;k<NK-1;k++)
									for(j=0;j<NJ-1;j++)
									{
										*((zone->grid+gdi)->SlabV[1]+(NJ-1)*k+j) = 0.0;//*((zone->grid+gdi)->V+(NI-1)*(NJ-1)*k+(NI-1)*j+0);
									}
								break;
								case 2:
									for(k=0;k<NK-1;k++)
									for(i=0;i<NI-1;i++)
									{
										*((zone->grid+gdi)->SlabV[2]+(NI-1)*k+i) = 0.0;//*((zone->grid+gdi)->V+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i);
									}
								break;
								case 3:
									for(k=0;k<NK-1;k++)
									for(i=0;i<NI-1;i++)
									{
										*((zone->grid+gdi)->SlabV[3]+(NI-1)*k+i) = 0.0;//*((zone->grid+gdi)->V+(NI-1)*(NJ-1)*k+(NI-1)*0+i);
									}
								break;
								case 4:
									for(j=0;j<NJ-1;j++)
									for(i=0;i<NI-1;i++)
									{
										*((zone->grid+gdi)->SlabV[4]+(NI-1)*j+i) = *((zone->grid+gdi)->V+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i);
									}
								break;
								case 5:
									for(j=0;j<NJ-1;j++)
									for(i=0;i<NI-1;i++)
									{
										*((zone->grid+gdi)->SlabV[5]+(NI-1)*j+i) = *((zone->grid+gdi)->V+(NI-1)*(NJ-1)*0+(NI-1)*j+i);
									}
								break;
							}
						break;
						case 2://Non-Homogeneous Dirichlet
							switch(iface)
							{
								case 0:
									for(k=0;k<NK-1;k++)
									for(j=0;j<NJ-1;j++)
									{
										*((zone->grid+gdi)->SlabV[0]+(NJ-1)*k+j) = 0.0;
									}
								break;
								case 1:
									for(k=0;k<NK-1;k++)
									for(j=0;j<NJ-1;j++)
									{
										*((zone->grid+gdi)->SlabV[1]+(NJ-1)*k+j) = 0.0;
									}
								break;
								case 2:
									for(k=0;k<NK-1;k++)
									for(i=0;i<NI-1;i++)
									{
										*((zone->grid+gdi)->SlabV[2]+(NI-1)*k+i) = 0.0;
									}
								break;
								case 3:
									for(k=0;k<NK-1;k++)
									for(i=0;i<NI-1;i++)
									{
										*((zone->grid+gdi)->SlabV[3]+(NI-1)*k+i) = 0.0;
									}
								break;
								case 4:
									for(j=0;j<NJ-1;j++)
									for(i=0;i<NI-1;i++)
									{
										*((zone->grid+gdi)->SlabV[4]+(NI-1)*j+i) = *((zone->grid+gdi)->V+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i);
									}
								break;
								case 5:
									for(j=0;j<NJ-1;j++)
									for(i=0;i<NI-1;i++)
									{
										*((zone->grid+gdi)->SlabV[5]+(NI-1)*j+i) = *((zone->grid+gdi)->V+(NI-1)*(NJ-1)*0+(NI-1)*j+i);
									}
								break;
							}
						break;
						case 3://Homogeneous Dirichlet
							switch(iface)
							{
								case 0:
									for(k=0;k<NK-1;k++)
									for(j=0;j<NJ-1;j++)
									{
										*((zone->grid+gdi)->SlabV[0]+(NJ-1)*k+j) = 0.0;
									}
								break;
								case 1:
									for(k=0;k<NK-1;k++)
									for(j=0;j<NJ-1;j++)
									{
										*((zone->grid+gdi)->SlabV[1]+(NJ-1)*k+j) = 0.0;
									}
								break;
								case 2:
									for(k=0;k<NK-1;k++)
									for(i=0;i<NI-1;i++)
									{
										*((zone->grid+gdi)->SlabV[2]+(NI-1)*k+i) = 0.0;
									}
								break;
								case 3:
									for(k=0;k<NK-1;k++)
									for(i=0;i<NI-1;i++)
									{
										*((zone->grid+gdi)->SlabV[3]+(NI-1)*k+i) = 0.0;
									}
								break;
								case 4:
									for(j=0;j<NJ-1;j++)
									for(i=0;i<NI-1;i++)
									{
										*((zone->grid+gdi)->SlabV[4]+(NI-1)*j+i) = *((zone->grid+gdi)->V+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i);
									}
								break;
								case 5:
									for(j=0;j<NJ-1;j++)
									for(i=0;i<NI-1;i++)
									{
										*((zone->grid+gdi)->SlabV[5]+(NI-1)*j+i) = *((zone->grid+gdi)->V+(NI-1)*(NJ-1)*0+(NI-1)*j+i);
									}
								break;
							}
						break;
					}
				}
				else
				{
					switch(c)
					{
						case 0://Homogeneous Neuman
							switch(iface)
							{
								case 0:
									for(k=0;k<NK-1;k++)
									for(j=0;j<NJ-1;j++)
									{
										*((zone->grid+gdi)->SlabV[0]+(NJ-1)*k+j) = *((zone->grid+gdi)->V+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2);
									}
								break;
								case 1:
									for(k=0;k<NK-1;k++)
									for(j=0;j<NJ-1;j++)
									{
										*((zone->grid+gdi)->SlabV[1]+(NJ-1)*k+j) = *((zone->grid+gdi)->V+(NI-1)*(NJ-1)*k+(NI-1)*j+0);
									}
								break;
								case 2:
									for(k=0;k<NK-1;k++)
									for(i=0;i<NI-1;i++)
									{
										*((zone->grid+gdi)->SlabV[2]+(NI-1)*k+i) = *((zone->grid+gdi)->V+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i);
									}
								break;
								case 3:
									for(k=0;k<NK-1;k++)
									for(i=0;i<NI-1;i++)
									{
										*((zone->grid+gdi)->SlabV[3]+(NI-1)*k+i) = *((zone->grid+gdi)->V+(NI-1)*(NJ-1)*k+(NI-1)*0+i);
									}
								break;
								case 4:
									for(j=0;j<NJ-1;j++)
									for(i=0;i<NI-1;i++)
									{
										*((zone->grid+gdi)->SlabV[4]+(NI-1)*j+i) = *((zone->grid+gdi)->V+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i);
									}
								break;
								case 5:
									for(j=0;j<NJ-1;j++)
									for(i=0;i<NI-1;i++)
									{
										*((zone->grid+gdi)->SlabV[5]+(NI-1)*j+i) = *((zone->grid+gdi)->V+(NI-1)*(NJ-1)*0+(NI-1)*j+i);
									}
								break;
							}
						break;
						case 2://Non-//Homogeneous Dirichlet
							switch(iface)
							{
								case 0:
									for(k=0;k<NK-1;k++)
									for(j=0;j<NJ-1;j++)
									{
										*((zone->grid+gdi)->SlabV[0]+(NJ-1)*k+j) = 0.0;
									}
								break;
								case 1:
									for(k=0;k<NK-1;k++)
									for(j=0;j<NJ-1;j++)
									{
										*((zone->grid+gdi)->SlabV[1]+(NJ-1)*k+j) = 0.0;
									}
								break;
								case 2:
									for(k=0;k<NK-1;k++)
									for(i=0;i<NI-1;i++)
									{
										*((zone->grid+gdi)->SlabV[2]+(NI-1)*k+i) = 0.0;
									}
								break;
								case 3:
									for(k=0;k<NK-1;k++)
									for(i=0;i<NI-1;i++)
									{
										*((zone->grid+gdi)->SlabV[3]+(NI-1)*k+i) = 0.0;
									}
								break;
								case 4:
									for(j=0;j<NJ-1;j++)
									for(i=0;i<NI-1;i++)
									{
										*((zone->grid+gdi)->SlabV[4]+(NI-1)*j+i) = 0.0;
									}
								break;
								case 5:
									for(j=0;j<NJ-1;j++)
									for(i=0;i<NI-1;i++)
									{
										*((zone->grid+gdi)->SlabV[5]+(NI-1)*j+i) = 0.0;
									}
								break;
							}
						break;
						case 3://Homogeneous Dirichlet
							switch(iface)
							{
								case 0:
									for(k=0;k<NK-1;k++)
									for(j=0;j<NJ-1;j++)
									{
										*((zone->grid+gdi)->SlabV[0]+(NJ-1)*k+j) = 0.0;
									}
								break;
								case 1:
									for(k=0;k<NK-1;k++)
									for(j=0;j<NJ-1;j++)
									{
										*((zone->grid+gdi)->SlabV[1]+(NJ-1)*k+j) = 0.0;
									}
								break;
								case 2:
									for(k=0;k<NK-1;k++)
									for(i=0;i<NI-1;i++)
									{
										*((zone->grid+gdi)->SlabV[2]+(NI-1)*k+i) = 0.0;
									}
								break;
								case 3:
									for(k=0;k<NK-1;k++)
									for(i=0;i<NI-1;i++)
									{
										*((zone->grid+gdi)->SlabV[3]+(NI-1)*k+i) = 0.0;
									}
								break;
								case 4:
									for(j=0;j<NJ-1;j++)
									for(i=0;i<NI-1;i++)
									{
										*((zone->grid+gdi)->SlabV[4]+(NI-1)*j+i) = 0.0;
									}
								break;
								case 5:
									for(j=0;j<NJ-1;j++)
									for(i=0;i<NI-1;i++)
									{
										*((zone->grid+gdi)->SlabV[5]+(NI-1)*j+i) = 0.0;
									}
								break;
							}
						break;
					}
				}
				
			}
		}
	}
	return 0;
}
