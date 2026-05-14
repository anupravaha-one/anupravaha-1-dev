#include "APKDummyVariableRestrict2d.h"

int DummyVariableRestrict2d(int variableTag,Base *base,Zone *zone,int gdi)
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
	
	for(iface=0;iface<4;iface++)
	{
		ptr = (void *)FindInterfaceBCPtrMod(zone,iface+1,name,&tag);
		if(tag=='i')
		{
			switch(iface)
			{
				case 0:
					k = 0;kc = 0;
					if((NJ-1)%2 == 0)
					{
						for(jc=0;jc<NJC-1;jc++)
						{
							j = 2*jc;
							*((zone->grid+gdi+1)->DummySlabV[1]+(NJC-1)*kc+jc) = 0.5*(*((zone->grid+gdi)->DummySlabV[1]+(NJ-1)*k+j) + *((zone->grid+gdi)->DummySlabV[1]+(NJ-1)*k+j+1));
						}
					}
					else
					{
						jc = NJC-2;j = NJ-2;
						*((zone->grid+gdi+1)->DummySlabV[1]+(NJC-1)*kc+jc) = *((zone->grid+gdi)->DummySlabV[1]+(NJ-1)*k+j);
						for(jc=0;jc<NJC-2;jc++)
						{
							j = 2*jc;
							*((zone->grid+gdi+1)->DummySlabV[1]+(NJC-1)*kc+jc) = 0.5*(*((zone->grid+gdi)->DummySlabV[1]+(NJ-1)*k+j) + *((zone->grid+gdi)->DummySlabV[1]+(NJ-1)*k+j+1));
						}
					}
				break;
				case 1:
					k = 0;kc = 0;
					if((NJ-1)%2 == 0)
					{
						for(jc=0;jc<NJC-1;jc++)
						{
							j = 2*jc;
							*((zone->grid+gdi+1)->DummySlabV[0]+(NJC-1)*kc+jc) = 0.5*(*((zone->grid+gdi)->DummySlabV[0]+(NJ-1)*k+j) + *((zone->grid+gdi)->DummySlabV[0]+(NJ-1)*k+j+1));
						}
					}
					else
					{
						jc = NJC-2;j = NJ-2;
						*((zone->grid+gdi+1)->DummySlabV[0]+(NJC-1)*kc+jc) = *((zone->grid+gdi)->DummySlabV[0]+(NJ-1)*k+j);
						for(jc=0;jc<NJC-2;jc++)
						{
							j = 2*jc;
							*((zone->grid+gdi+1)->DummySlabV[0]+(NJC-1)*kc+jc) = 0.5*(*((zone->grid+gdi)->DummySlabV[0]+(NJ-1)*k+j) + *((zone->grid+gdi)->DummySlabV[0]+(NJ-1)*k+j+1));
						}
					}
				break;
				case 2:
					k = 0;kc = 0;
					if((NI-1)%2 == 0)
					{
						for(ic=0;ic<NIC-1;ic++)
						{
							i = 2*ic;
							*((zone->grid+gdi+1)->DummySlabV[3]+(NIC-1)*kc+ic) = 0.5*(*((zone->grid+gdi)->DummySlabV[3]+(NI-1)*k+i)+*((zone->grid+gdi)->DummySlabV[3]+(NI-1)*k+i+1));
						}
					}
					else
					{
						ic = NIC-2;i = NI-2;
						*((zone->grid+gdi+1)->DummySlabV[3]+(NIC-1)*kc+ic) = *((zone->grid+gdi)->DummySlabV[3]+(NI-1)*k+i);
						for(ic=0;ic<NIC-2;ic++)
						{
							i = 2*ic;
							*((zone->grid+gdi+1)->DummySlabV[3]+(NIC-1)*kc+ic) = 0.5*(*((zone->grid+gdi)->DummySlabV[3]+(NI-1)*k+i)+*((zone->grid+gdi)->DummySlabV[3]+(NI-1)*k+i+1));
						}
					}
				break;
				case 3:
					k = 0;kc = 0;
					if((NI-1)%2 == 0)
					{
						for(ic=0;ic<NIC-1;ic++)
						{
							i = 2*ic;
							*((zone->grid+gdi+1)->DummySlabV[2]+(NIC-1)*kc+ic) = 0.5*(*((zone->grid+gdi)->DummySlabV[2]+(NI-1)*k+i)+*((zone->grid+gdi)->DummySlabV[2]+(NI-1)*k+i+1));
						}
					}
					else
					{
						ic = NIC-2;i = NI-2;
						*((zone->grid+gdi+1)->DummySlabV[2]+(NIC-1)*kc+ic) = *((zone->grid+gdi)->DummySlabV[2]+(NI-1)*k+i);
						for(ic=0;ic<NIC-2;ic++)
						{
							i = 2*ic;
							*((zone->grid+gdi+1)->DummySlabV[2]+(NIC-1)*kc+ic) = 0.5*(*((zone->grid+gdi)->DummySlabV[2]+(NI-1)*k+i)+*((zone->grid+gdi)->DummySlabV[2]+(NI-1)*k+i+1));
						}
					}
				break;
			}
		}
	}
	return 0;
}
