#include "APKDummyVariableRestrict.h"

int DummyVariableRestrict(int variableTag,Base *base,Zone *zone,int gdi)
{
	int iface;
	char tag,name[33];
	void *ptr;
	int i,j,k,ic,jc,kc,NI,NJ,NK,NIC,NJC,NKC;
	
	NI = *((zone->grid+gdi)->GSize+0);
	NJ = *((zone->grid+gdi)->GSize+1);
	NK = *((zone->grid+gdi)->GSize+2);
	
	NIC = *((zone->grid+gdi+1)->GSize+0);
	NJC = *((zone->grid+gdi+1)->GSize+1);
	NKC = *((zone->grid+gdi+1)->GSize+2);
	
	for(iface=0;iface<6;iface++)
	{
		ptr = (void *)FindInterfaceBCPtrMod(zone,iface+1,name,&tag);
		if(tag=='i')
		{
			switch(iface)
			{
				case 0:
					if(((NJ-1)%2 == 0)&&((NK-1)%2 == 0))
					{
						for(kc=0;kc<NKC-1;kc++)
						{
							k = 2*kc;
							for(jc=0;jc<NJC-1;jc++)
							{
								j = 2*jc;
								*((zone->grid+gdi+1)->DummySlabV[1]+(NJC-1)*kc+jc) = 0.25*(*((zone->grid+gdi)->DummySlabV[1]+(NJ-1)*k+j) + *((zone->grid+gdi)->DummySlabV[1]+(NJ-1)*k+j+1)
															 + *((zone->grid+gdi)->DummySlabV[1]+(NJ-1)*(k+1)+j) + *((zone->grid+gdi)->DummySlabV[1]+(NJ-1)*(k+1)+j+1));
							}
						}
					}
					else if(((NJ-1)%2 != 0)&&((NK-1)%2 == 0))
					{
						jc = NJC-2;j = NJ-2;
						for(kc=0;kc<NKC-1;kc++)
						{
							k = 2*kc;
							*((zone->grid+gdi+1)->DummySlabV[1]+(NJC-1)*kc+jc) = 0.5*(*((zone->grid+gdi)->DummySlabV[1]+(NJ-1)*k+j) + *((zone->grid+gdi)->DummySlabV[1]+(NJ-1)*(k+1)+j));
						}
						for(kc=0;kc<NKC-1;kc++)
						{
							k = 2*kc;
							for(jc=0;jc<NJC-2;jc++)
							{
								j = 2*jc;
								*((zone->grid+gdi+1)->DummySlabV[1]+(NJC-1)*kc+jc) = 0.25*(*((zone->grid+gdi)->DummySlabV[1]+(NJ-1)*k+j) + *((zone->grid+gdi)->DummySlabV[1]+(NJ-1)*k+j+1)
															 + *((zone->grid+gdi)->DummySlabV[1]+(NJ-1)*(k+1)+j) + *((zone->grid+gdi)->DummySlabV[1]+(NJ-1)*(k+1)+j+1));
							}
						}
					}
					else if(((NJ-1)%2 == 0)&&((NK-1)%2 != 0))
					{
						kc = NKC-2;k = NK-2;
						for(jc=0;jc<NJC-1;jc++)
						{
							j = 2*jc;
							*((zone->grid+gdi+1)->DummySlabV[1]+(NJC-1)*kc+jc) = 0.5*(*((zone->grid+gdi)->DummySlabV[1]+(NJ-1)*k+j) + *((zone->grid+gdi)->DummySlabV[1]+(NJ-1)*k+j+1));
						}
						for(kc=0;kc<NKC-2;kc++)
						{
							k = 2*kc;
							for(jc=0;jc<NJC-1;jc++)
							{
								j = 2*jc;
								*((zone->grid+gdi+1)->DummySlabV[1]+(NJC-1)*kc+jc) = 0.25*(*((zone->grid+gdi)->DummySlabV[1]+(NJ-1)*k+j) + *((zone->grid+gdi)->DummySlabV[1]+(NJ-1)*k+j+1)
															 + *((zone->grid+gdi)->DummySlabV[1]+(NJ-1)*(k+1)+j) + *((zone->grid+gdi)->DummySlabV[1]+(NJ-1)*(k+1)+j+1));
							}
						}
					}
					else if(((NJ-1)%2 != 0)&&((NK-1)%2 != 0))
					{
						jc = NJC-2;j = NJ-2;
						kc = NKC-2;k = NK-2;
						*((zone->grid+gdi+1)->DummySlabV[1]+(NJC-1)*kc+jc) = *((zone->grid+gdi)->DummySlabV[1]+(NJ-1)*k+j);
						jc = NJC-2;j = NJ-2;
						for(kc=0;kc<NKC-2;kc++)
						{
							k = 2*kc;
							*((zone->grid+gdi+1)->DummySlabV[1]+(NJC-1)*kc+jc) = 0.5*(*((zone->grid+gdi)->DummySlabV[1]+(NJ-1)*k+j) + *((zone->grid+gdi)->DummySlabV[1]+(NJ-1)*(k+1)+j));
						}
						kc = NKC-2;k = NK-2;
						for(jc=0;jc<NJC-2;jc++)
						{
							j = 2*jc;
							*((zone->grid+gdi+1)->DummySlabV[1]+(NJC-1)*kc+jc) = 0.5*(*((zone->grid+gdi)->DummySlabV[1]+(NJ-1)*k+j) + *((zone->grid+gdi)->DummySlabV[1]+(NJ-1)*k+j+1));
						}
						for(kc=0;kc<NKC-2;kc++)
						{
							k = 2*kc;
							for(jc=0;jc<NJC-2;jc++)
							{
								j = 2*jc;
								*((zone->grid+gdi+1)->DummySlabV[1]+(NJC-1)*kc+jc) = 0.25*(*((zone->grid+gdi)->DummySlabV[1]+(NJ-1)*k+j) + *((zone->grid+gdi)->DummySlabV[1]+(NJ-1)*k+j+1)
															 + *((zone->grid+gdi)->DummySlabV[1]+(NJ-1)*(k+1)+j) + *((zone->grid+gdi)->DummySlabV[1]+(NJ-1)*(k+1)+j+1));
							}
						}
					}
				break;
				case 1:
					if(((NJ-1)%2 == 0)&&((NK-1)%2 == 0))
					{
						for(kc=0;kc<NKC-1;kc++)
						{
							k = 2*kc;
							for(jc=0;jc<NJC-1;jc++)
							{
								j = 2*jc;
								*((zone->grid+gdi+1)->DummySlabV[0]+(NJC-1)*kc+jc) = 0.25*(*((zone->grid+gdi)->DummySlabV[0]+(NJ-1)*k+j) + *((zone->grid+gdi)->DummySlabV[0]+(NJ-1)*k+j+1)
															 + *((zone->grid+gdi)->DummySlabV[0]+(NJ-1)*(k+1)+j) + *((zone->grid+gdi)->DummySlabV[0]+(NJ-1)*(k+1)+j+1));
							}
						}
					}
					else if(((NJ-1)%2 != 0)&&((NK-1)%2 == 0))
					{
						jc = NJC-2;j = NJ-2;
						for(kc=0;kc<NKC-1;kc++)
						{
							k = 2*kc;
							*((zone->grid+gdi+1)->DummySlabV[0]+(NJC-1)*kc+jc) = 0.5*(*((zone->grid+gdi)->DummySlabV[0]+(NJ-1)*k+j) + *((zone->grid+gdi)->DummySlabV[0]+(NJ-1)*(k+1)+j));
						}
						for(kc=0;kc<NKC-1;kc++)
						{
							k = 2*kc;
							for(jc=0;jc<NJC-2;jc++)
							{
								j = 2*jc;
								*((zone->grid+gdi+1)->DummySlabV[0]+(NJC-1)*kc+jc) = 0.25*(*((zone->grid+gdi)->DummySlabV[0]+(NJ-1)*k+j) + *((zone->grid+gdi)->DummySlabV[0]+(NJ-1)*k+j+1)
															 + *((zone->grid+gdi)->DummySlabV[0]+(NJ-1)*(k+1)+j) + *((zone->grid+gdi)->DummySlabV[0]+(NJ-1)*(k+1)+j+1));
							}
						}
					}
					else if(((NJ-1)%2 == 0)&&((NK-1)%2 != 0))
					{
						kc = NKC-2;k = NK-2;
						for(jc=0;jc<NJC-1;jc++)
						{
							j = 2*jc;
							*((zone->grid+gdi+1)->DummySlabV[0]+(NJC-1)*kc+jc) = 0.5*(*((zone->grid+gdi)->DummySlabV[0]+(NJ-1)*k+j) + *((zone->grid+gdi)->DummySlabV[0]+(NJ-1)*k+j+1));
						}
						for(kc=0;kc<NKC-2;kc++)
						{
							k = 2*kc;
							for(jc=0;jc<NJC-1;jc++)
							{
								j = 2*jc;
								*((zone->grid+gdi+1)->DummySlabV[0]+(NJC-1)*kc+jc) = 0.25*(*((zone->grid+gdi)->DummySlabV[0]+(NJ-1)*k+j) + *((zone->grid+gdi)->DummySlabV[0]+(NJ-1)*k+j+1)
															 + *((zone->grid+gdi)->DummySlabV[0]+(NJ-1)*(k+1)+j) + *((zone->grid+gdi)->DummySlabV[0]+(NJ-1)*(k+1)+j+1));
							}
						}
					}
					else if(((NJ-1)%2 != 0)&&((NK-1)%2 != 0))
					{
						jc = NJC-2;j = NJ-2;
						kc = NKC-2;k = NK-2;
						*((zone->grid+gdi+1)->DummySlabV[0]+(NJC-1)*kc+jc) = *((zone->grid+gdi)->DummySlabV[0]+(NJ-1)*k+j);
						jc = NJC-2;j = NJ-2;
						for(kc=0;kc<NKC-2;kc++)
						{
							k = 2*kc;
							*((zone->grid+gdi+1)->DummySlabV[0]+(NJC-1)*kc+jc) = 0.5*(*((zone->grid+gdi)->DummySlabV[0]+(NJ-1)*k+j) + *((zone->grid+gdi)->DummySlabV[0]+(NJ-1)*(k+1)+j));
						}
						kc = NKC-2;k = NK-2;
						for(jc=0;jc<NJC-2;jc++)
						{
							j = 2*jc;
							*((zone->grid+gdi+1)->DummySlabV[0]+(NJC-1)*kc+jc) = 0.5*(*((zone->grid+gdi)->DummySlabV[0]+(NJ-1)*k+j) + *((zone->grid+gdi)->DummySlabV[0]+(NJ-1)*k+j+1));
						}
						for(kc=0;kc<NKC-2;kc++)
						{
							k = 2*kc;
							for(jc=0;jc<NJC-2;jc++)
							{
								j = 2*jc;
								*((zone->grid+gdi+1)->DummySlabV[0]+(NJC-1)*kc+jc) = 0.25*(*((zone->grid+gdi)->DummySlabV[0]+(NJ-1)*k+j) + *((zone->grid+gdi)->DummySlabV[0]+(NJ-1)*k+j+1)
															 + *((zone->grid+gdi)->DummySlabV[0]+(NJ-1)*(k+1)+j) + *((zone->grid+gdi)->DummySlabV[0]+(NJ-1)*(k+1)+j+1));
							}
						}
					}
				break;
				case 2:
					if(((NI-1)%2 == 0)&&((NK-1)%2 == 0))
					{
						for(kc=0;kc<NKC-1;kc++)
						{
							k = 2*kc;
							for(ic=0;ic<NIC-1;ic++)
							{
								i = 2*ic;
								*((zone->grid+gdi+1)->DummySlabV[3]+(NIC-1)*kc+ic) = 0.25*(*((zone->grid+gdi)->DummySlabV[3]+(NI-1)*k+i)+*((zone->grid+gdi)->DummySlabV[3]+(NI-1)*k+i+1)
															+*((zone->grid+gdi)->DummySlabV[3]+(NI-1)*(k+1)+i)+*((zone->grid+gdi)->DummySlabV[3]+(NI-1)*(k+1)+i+1));
							}
						}
					}
					else if(((NI-1)%2 != 0)&&((NK-1)%2 == 0))
					{
						ic = NIC-2;i = NI-2;
						for(kc=0;kc<NKC-1;kc++)
						{
							k = 2*kc;
							*((zone->grid+gdi+1)->DummySlabV[3]+(NIC-1)*kc+ic) = 0.5*(*((zone->grid+gdi)->DummySlabV[3]+(NI-1)*k+i)+*((zone->grid+gdi)->DummySlabV[3]+(NI-1)*(k+1)+i));
						}
						for(kc=0;kc<NKC-1;kc++)
						{
							k = 2*kc;
							for(ic=0;ic<NIC-2;ic++)
							{
								i = 2*ic;
								*((zone->grid+gdi+1)->DummySlabV[3]+(NIC-1)*kc+ic) = 0.25*(*((zone->grid+gdi)->DummySlabV[3]+(NI-1)*k+i)+*((zone->grid+gdi)->DummySlabV[3]+(NI-1)*k+i+1)
															+*((zone->grid+gdi)->DummySlabV[3]+(NI-1)*(k+1)+i)+*((zone->grid+gdi)->DummySlabV[3]+(NI-1)*(k+1)+i+1));
							}
						}
					}
					else if(((NI-1)%2 == 0)&&((NK-1)%2 != 0))
					{
						kc = NKC-2;k = NK-2;
						for(ic=0;ic<NIC-1;ic++)
						{
							i = 2*ic;
							*((zone->grid+gdi+1)->DummySlabV[3]+(NIC-1)*kc+ic) = 0.5*(*((zone->grid+gdi)->DummySlabV[3]+(NI-1)*k+i)+*((zone->grid+gdi)->DummySlabV[3]+(NI-1)*k+i+1));
						}
						for(kc=0;kc<NKC-2;kc++)
						{
							k = 2*kc;
							for(ic=0;ic<NIC-1;ic++)
							{
								i = 2*ic;
								*((zone->grid+gdi+1)->DummySlabV[3]+(NIC-1)*kc+ic) = 0.25*(*((zone->grid+gdi)->DummySlabV[3]+(NI-1)*k+i)+*((zone->grid+gdi)->DummySlabV[3]+(NI-1)*k+i+1)
															+*((zone->grid+gdi)->DummySlabV[3]+(NI-1)*(k+1)+i)+*((zone->grid+gdi)->DummySlabV[3]+(NI-1)*(k+1)+i+1));
							}
						}
					}
					else if(((NI-1)%2 != 0)&&((NK-1)%2 != 0))
					{
						ic = NIC-2;i = NI-2;
						kc = NKC-2;k = NK-2;
						*((zone->grid+gdi+1)->DummySlabV[3]+(NIC-1)*kc+ic) = *((zone->grid+gdi)->DummySlabV[3]+(NI-1)*k+i);
						ic = NIC-2;i = NI-2;
						for(kc=0;kc<NKC-2;kc++)
						{
							k = 2*kc;
							*((zone->grid+gdi+1)->DummySlabV[3]+(NIC-1)*kc+ic) = 0.5*(*((zone->grid+gdi)->DummySlabV[3]+(NI-1)*k+i)+*((zone->grid+gdi)->DummySlabV[3]+(NI-1)*(k+1)+i));
						}
						kc = NKC-2;k = NK-2;
						for(ic=0;ic<NIC-2;ic++)
						{
							i = 2*ic;
							*((zone->grid+gdi+1)->DummySlabV[3]+(NIC-1)*kc+ic) = 0.5*(*((zone->grid+gdi)->DummySlabV[3]+(NI-1)*k+i)+*((zone->grid+gdi)->DummySlabV[3]+(NI-1)*k+i+1));
						}
						for(kc=0;kc<NKC-2;kc++)
						{
							k = 2*kc;
							for(ic=0;ic<NIC-2;ic++)
							{
								i = 2*ic;
								*((zone->grid+gdi+1)->DummySlabV[3]+(NIC-1)*kc+ic) = 0.25*(*((zone->grid+gdi)->DummySlabV[3]+(NI-1)*k+i)+*((zone->grid+gdi)->DummySlabV[3]+(NI-1)*k+i+1)
															+*((zone->grid+gdi)->DummySlabV[3]+(NI-1)*(k+1)+i)+*((zone->grid+gdi)->DummySlabV[3]+(NI-1)*(k+1)+i+1));
							}
						}
					}
				break;
				case 3:
					if(((NI-1)%2 == 0)&&((NK-1)%2 == 0))
					{
						for(kc=0;kc<NKC-1;kc++)
						{
							k = 2*kc;
							for(ic=0;ic<NIC-1;ic++)
							{
								i = 2*ic;
								*((zone->grid+gdi+1)->DummySlabV[2]+(NIC-1)*kc+ic) = 0.25*(*((zone->grid+gdi)->DummySlabV[2]+(NI-1)*k+i)+*((zone->grid+gdi)->DummySlabV[2]+(NI-1)*k+i+1)
															+*((zone->grid+gdi)->DummySlabV[2]+(NI-1)*(k+1)+i)+*((zone->grid+gdi)->DummySlabV[2]+(NI-1)*(k+1)+i+1));
							}
						}
					}
					else if(((NI-1)%2 != 0)&&((NK-1)%2 == 0))
					{
						ic = NIC-2;i = NI-2;
						for(kc=0;kc<NKC-1;kc++)
						{
							k = 2*kc;
							*((zone->grid+gdi+1)->DummySlabV[2]+(NIC-1)*kc+ic) = 0.5*(*((zone->grid+gdi)->DummySlabV[2]+(NI-1)*k+i)+*((zone->grid+gdi)->DummySlabV[2]+(NI-1)*(k+1)+i));
						}
						for(kc=0;kc<NKC-1;kc++)
						{
							k = 2*kc;
							for(ic=0;ic<NIC-2;ic++)
							{
								i = 2*ic;
								*((zone->grid+gdi+1)->DummySlabV[2]+(NIC-1)*kc+ic) = 0.25*(*((zone->grid+gdi)->DummySlabV[2]+(NI-1)*k+i)+*((zone->grid+gdi)->DummySlabV[2]+(NI-1)*k+i+1)
															+*((zone->grid+gdi)->DummySlabV[2]+(NI-1)*(k+1)+i)+*((zone->grid+gdi)->DummySlabV[2]+(NI-1)*(k+1)+i+1));
							}
						}
					}
					else if(((NI-1)%2 == 0)&&((NK-1)%2 != 0))
					{
						kc = NKC-2;k = NK-2;
						for(ic=0;ic<NIC-1;ic++)
						{
							i = 2*ic;
							*((zone->grid+gdi+1)->DummySlabV[2]+(NIC-1)*kc+ic) = 0.5*(*((zone->grid+gdi)->DummySlabV[2]+(NI-1)*k+i)+*((zone->grid+gdi)->DummySlabV[2]+(NI-1)*k+i+1));
						}
						for(kc=0;kc<NKC-2;kc++)
						{
							k = 2*kc;
							for(ic=0;ic<NIC-1;ic++)
							{
								i = 2*ic;
								*((zone->grid+gdi+1)->DummySlabV[2]+(NIC-1)*kc+ic) = 0.25*(*((zone->grid+gdi)->DummySlabV[2]+(NI-1)*k+i)+*((zone->grid+gdi)->DummySlabV[2]+(NI-1)*k+i+1)
															+*((zone->grid+gdi)->DummySlabV[2]+(NI-1)*(k+1)+i)+*((zone->grid+gdi)->DummySlabV[2]+(NI-1)*(k+1)+i+1));
							}
						}
					}
					else if(((NI-1)%2 != 0)&&((NK-1)%2 != 0))
					{
						ic = NIC-2;i = NI-2;
						kc = NKC-2;k = NK-2;
						*((zone->grid+gdi+1)->DummySlabV[2]+(NIC-1)*kc+ic) = *((zone->grid+gdi)->DummySlabV[2]+(NI-1)*k+i);
						ic = NIC-2;i = NI-2;
						for(kc=0;kc<NKC-2;kc++)
						{
							k = 2*kc;
							*((zone->grid+gdi+1)->DummySlabV[2]+(NIC-1)*kc+ic) = 0.5*(*((zone->grid+gdi)->DummySlabV[2]+(NI-1)*k+i)+*((zone->grid+gdi)->DummySlabV[2]+(NI-1)*(k+1)+i));
						}
						kc = NKC-2;k = NK-2;
						for(ic=0;ic<NIC-2;ic++)
						{
							i = 2*ic;
							*((zone->grid+gdi+1)->DummySlabV[2]+(NIC-1)*kc+ic) = 0.5*(*((zone->grid+gdi)->DummySlabV[2]+(NI-1)*k+i)+*((zone->grid+gdi)->DummySlabV[2]+(NI-1)*k+i+1));
						}
						for(kc=0;kc<NKC-2;kc++)
						{
							k = 2*kc;
							for(ic=0;ic<NIC-2;ic++)
							{
								i = 2*ic;
								*((zone->grid+gdi+1)->DummySlabV[2]+(NIC-1)*kc+ic) = 0.25*(*((zone->grid+gdi)->DummySlabV[2]+(NI-1)*k+i)+*((zone->grid+gdi)->DummySlabV[2]+(NI-1)*k+i+1)
															+*((zone->grid+gdi)->DummySlabV[2]+(NI-1)*(k+1)+i)+*((zone->grid+gdi)->DummySlabV[2]+(NI-1)*(k+1)+i+1));
							}
						}
					}
				break;
				case 4:
					if(((NI-1)%2 == 0)&&((NJ-1)%2 == 0))
					{
						for(jc=0;jc<NJC-1;jc++)
						{
							j = 2*jc;
							for(ic=0;ic<NIC-1;ic++)
							{
								i = 2*ic;
								*((zone->grid+gdi+1)->DummySlabV[5]+(NIC-1)*jc+ic) = 0.25*(*((zone->grid+gdi)->DummySlabV[5]+(NI-1)*j+i)+*((zone->grid+gdi)->DummySlabV[5]+(NI-1)*j+i+1)
															+*((zone->grid+gdi)->DummySlabV[5]+(NI-1)*(j+1)+i)+*((zone->grid+gdi)->DummySlabV[5]+(NI-1)*(j+1)+i+1));
							}
						}
					}
					else if(((NI-1)%2 != 0)&&((NJ-1)%2 == 0))
					{
						ic = NIC-2;i = NI-2;
						for(jc=0;jc<NJC-1;jc++)
						{
							j = 2*jc;
							*((zone->grid+gdi+1)->DummySlabV[5]+(NIC-1)*jc+ic) = 0.5*(*((zone->grid+gdi)->DummySlabV[5]+(NI-1)*j+i)+*((zone->grid+gdi)->DummySlabV[5]+(NI-1)*(j+1)+i));
						}
						for(jc=0;jc<NJC-1;jc++)
						{
							j = 2*jc;
							for(ic=0;ic<NIC-2;ic++)
							{
								i = 2*ic;
								*((zone->grid+gdi+1)->DummySlabV[5]+(NIC-1)*jc+ic) = 0.25*(*((zone->grid+gdi)->DummySlabV[5]+(NI-1)*j+i)+*((zone->grid+gdi)->DummySlabV[5]+(NI-1)*j+i+1)
															+*((zone->grid+gdi)->DummySlabV[5]+(NI-1)*(j+1)+i)+*((zone->grid+gdi)->DummySlabV[5]+(NI-1)*(j+1)+i+1));
							}
						}
					}
					else if(((NI-1)%2 == 0)&&((NJ-1)%2 != 0))
					{
						jc = NJC-2;j = NJ-2;
						for(ic=0;ic<NIC-1;ic++)
						{
							i = 2*ic;
							*((zone->grid+gdi+1)->DummySlabV[5]+(NIC-1)*jc+ic) = 0.5*(*((zone->grid+gdi)->DummySlabV[5]+(NI-1)*j+i)+*((zone->grid+gdi)->DummySlabV[5]+(NI-1)*j+i+1));
						}
						for(jc=0;jc<NJC-2;jc++)
						{
							j = 2*jc;
							for(ic=0;ic<NIC-1;ic++)
							{
								i = 2*ic;
								*((zone->grid+gdi+1)->DummySlabV[5]+(NIC-1)*jc+ic) = 0.25*(*((zone->grid+gdi)->DummySlabV[5]+(NI-1)*j+i)+*((zone->grid+gdi)->DummySlabV[5]+(NI-1)*j+i+1)
															+*((zone->grid+gdi)->DummySlabV[5]+(NI-1)*(j+1)+i)+*((zone->grid+gdi)->DummySlabV[5]+(NI-1)*(j+1)+i+1));
							}
						}
					}
					else if(((NI-1)%2 != 0)&&((NJ-1)%2 != 0))
					{
						jc = NJC-2;j = NJ-2;
						ic = NIC-2;i = NI-2;
						*((zone->grid+gdi+1)->DummySlabV[5]+(NIC-1)*jc+ic) = *((zone->grid+gdi)->DummySlabV[5]+(NI-1)*j+i);
						ic = NIC-2;i = NI-2;
						for(jc=0;jc<NJC-2;jc++)
						{
							j = 2*jc;
							*((zone->grid+gdi+1)->DummySlabV[5]+(NIC-1)*jc+ic) = 0.5*(*((zone->grid+gdi)->DummySlabV[5]+(NI-1)*j+i)+*((zone->grid+gdi)->DummySlabV[5]+(NI-1)*(j+1)+i));
						}
						jc = NJC-2;j = NJ-2;
						for(ic=0;ic<NIC-2;ic++)
						{
							i = 2*ic;
							*((zone->grid+gdi+1)->DummySlabV[5]+(NIC-1)*jc+ic) = 0.5*(*((zone->grid+gdi)->DummySlabV[5]+(NI-1)*j+i)+*((zone->grid+gdi)->DummySlabV[5]+(NI-1)*j+i+1));
						}
						for(jc=0;jc<NJC-2;jc++)
						{
							j = 2*jc;
							for(ic=0;ic<NIC-2;ic++)
							{
								i = 2*ic;
								*((zone->grid+gdi+1)->DummySlabV[5]+(NIC-1)*jc+ic) = 0.25*(*((zone->grid+gdi)->DummySlabV[5]+(NI-1)*j+i)+*((zone->grid+gdi)->DummySlabV[5]+(NI-1)*j+i+1)
															+*((zone->grid+gdi)->DummySlabV[5]+(NI-1)*(j+1)+i)+*((zone->grid+gdi)->DummySlabV[5]+(NI-1)*(j+1)+i+1));
							}
						}
					}
				break;
				case 5:
					if(((NI-1)%2 == 0)&&((NJ-1)%2 == 0))
					{
						for(jc=0;jc<NJC-1;jc++)
						{
							j = 2*jc;
							for(ic=0;ic<NIC-1;ic++)
							{
								i = 2*ic;
								*((zone->grid+gdi+1)->DummySlabV[4]+(NIC-1)*jc+ic) = 0.25*(*((zone->grid+gdi)->DummySlabV[4]+(NI-1)*j+i)+*((zone->grid+gdi)->DummySlabV[4]+(NI-1)*j+i+1)
															+*((zone->grid+gdi)->DummySlabV[4]+(NI-1)*(j+1)+i)+*((zone->grid+gdi)->DummySlabV[4]+(NI-1)*(j+1)+i+1));
							}
						}
					}
					else if(((NI-1)%2 != 0)&&((NJ-1)%2 == 0))
					{
						ic = NIC-2;i = NI-2;
						for(jc=0;jc<NJC-1;jc++)
						{
							j = 2*jc;
							*((zone->grid+gdi+1)->DummySlabV[4]+(NIC-1)*jc+ic) = 0.5*(*((zone->grid+gdi)->DummySlabV[4]+(NI-1)*j+i)+*((zone->grid+gdi)->DummySlabV[4]+(NI-1)*(j+1)+i));
						}
						for(jc=0;jc<NJC-1;jc++)
						{
							j = 2*jc;
							for(ic=0;ic<NIC-2;ic++)
							{
								i = 2*ic;
								*((zone->grid+gdi+1)->DummySlabV[4]+(NIC-1)*jc+ic) = 0.25*(*((zone->grid+gdi)->DummySlabV[4]+(NI-1)*j+i)+*((zone->grid+gdi)->DummySlabV[4]+(NI-1)*j+i+1)
															+*((zone->grid+gdi)->DummySlabV[4]+(NI-1)*(j+1)+i)+*((zone->grid+gdi)->DummySlabV[4]+(NI-1)*(j+1)+i+1));
							}
						}
					}
					else if(((NI-1)%2 == 0)&&((NJ-1)%2 != 0))
					{
						jc = NJC-2;j = NJ-2;
						for(ic=0;ic<NIC-1;ic++)
						{
							i = 2*ic;
							*((zone->grid+gdi+1)->DummySlabV[4]+(NIC-1)*jc+ic) = 0.5*(*((zone->grid+gdi)->DummySlabV[4]+(NI-1)*j+i)+*((zone->grid+gdi)->DummySlabV[4]+(NI-1)*j+i+1));
						}
						for(jc=0;jc<NJC-2;jc++)
						{
							j = 2*jc;
							for(ic=0;ic<NIC-1;ic++)
							{
								i = 2*ic;
								*((zone->grid+gdi+1)->DummySlabV[4]+(NIC-1)*jc+ic) = 0.25*(*((zone->grid+gdi)->DummySlabV[4]+(NI-1)*j+i)+*((zone->grid+gdi)->DummySlabV[4]+(NI-1)*j+i+1)
															+*((zone->grid+gdi)->DummySlabV[4]+(NI-1)*(j+1)+i)+*((zone->grid+gdi)->DummySlabV[4]+(NI-1)*(j+1)+i+1));
							}
						}
					}
					else if(((NI-1)%2 != 0)&&((NJ-1)%2 != 0))
					{
						jc = NJC-2;j = NJ-2;
						ic = NIC-2;i = NI-2;
						*((zone->grid+gdi+1)->DummySlabV[4]+(NIC-1)*jc+ic) = *((zone->grid+gdi)->DummySlabV[4]+(NI-1)*j+i);
						ic = NIC-2;i = NI-2;
						for(jc=0;jc<NJC-2;jc++)
						{
							j = 2*jc;
							*((zone->grid+gdi+1)->DummySlabV[4]+(NIC-1)*jc+ic) = 0.5*(*((zone->grid+gdi)->DummySlabV[4]+(NI-1)*j+i)+*((zone->grid+gdi)->DummySlabV[4]+(NI-1)*(j+1)+i));
						}
						jc = NJC-2;j = NJ-2;
						for(ic=0;ic<NIC-2;ic++)
						{
							i = 2*ic;
							*((zone->grid+gdi+1)->DummySlabV[4]+(NIC-1)*jc+ic) = 0.5*(*((zone->grid+gdi)->DummySlabV[4]+(NI-1)*j+i)+*((zone->grid+gdi)->DummySlabV[4]+(NI-1)*j+i+1));
						}
						for(jc=0;jc<NJC-2;jc++)
						{
							j = 2*jc;
							for(ic=0;ic<NIC-2;ic++)
							{
								i = 2*ic;
								*((zone->grid+gdi+1)->DummySlabV[4]+(NIC-1)*jc+ic) = 0.25*(*((zone->grid+gdi)->DummySlabV[4]+(NI-1)*j+i)+*((zone->grid+gdi)->DummySlabV[4]+(NI-1)*j+i+1)
															+*((zone->grid+gdi)->DummySlabV[4]+(NI-1)*(j+1)+i)+*((zone->grid+gdi)->DummySlabV[4]+(NI-1)*(j+1)+i+1));
							}
						}
					}
				break;
			}
		}
	}
	return 0;
}
