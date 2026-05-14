#include "APKDummyVariableProlongate2d.h"

int DummyVariableProlongate2d(int variableTag,Base *base,Zone *zone,int gdi)
{
	int iface;
	char tag,name[33];
	void *ptr;
	int i,j,k,ic,jc,kc,NI,NJ,NK,NIC,NJC,NKC;
	
	NI = *((zone->grid+gdi-1)->GSize+0);
	NJ = *((zone->grid+gdi-1)->GSize+1);
	NK = *((zone->grid+gdi-1)->GSize+2);
	
	NIC = *((zone->grid+gdi)->GSize+0);
	NJC = *((zone->grid+gdi)->GSize+1);
	NKC = *((zone->grid+gdi)->GSize+2);
	
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
						jc = 0;j = 0;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2) = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+NI-2) = 0.5*(*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2)+*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+NIC-2));
						jc = NJC-2;j = NJ-3;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2) = 0.5*(*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2)+*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+NIC-2));
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+NI-2) = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2);
						jc = 1;
						for(j=2;j<NJ-3;j=j+2)
						{
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2) = 0.5*(*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2)+*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+NIC-2));
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+NI-2) = 0.5*(*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2)+*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+NIC-2));
							jc++;
						}
					}
					else
					{
						jc = 0;j = 0;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2) = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+NI-2) = 0.5*(*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2)+*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+NIC-2));
						jc = NJC-2;j = NJ-2;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2) = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2);
						jc = 1;
						for(j=2;j<NJ-2;j=j+2)
						{
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2) = 0.5*(*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2)+*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+NIC-2));
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+NI-2) = 0.5*(*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+NIC-2)+*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+NIC-2));
							jc++;
						}
					}
				break;
				case 1:
					k = 0;kc = 0;
					if((NJ-1)%2 == 0)
					{
						jc = 0;j = 0;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+0) = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+0) = 0.5*(*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0)+*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+0));
						jc = NJC-2;j = NJ-3;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+0) = 0.5*(*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0)+*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+0));
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+0) = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0);
						jc = 1;
						for(j=2;j<NJ-3;j=j+2)
						{
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+0) = 0.5*(*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0)+*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+0));
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+0) = 0.5*(*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0)+*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+0));
							jc++;
						}
					}
					else
					{
						jc = 0;j = 0;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+0) = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+0) = 0.5*(*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0)+*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+0));
						jc = NJC-2;j = NJ-2;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+0) = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0);
						jc = 1;
						for(j=2;j<NJ-2;j=j+2)
						{
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*j+0) = 0.5*(*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0)+*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc-1)+0));
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+0) = 0.5*(*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*jc+0)+*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(jc+1)+0));
							jc++;
						}
					}
				break;
				case 2:
					k = 0;kc = 0;
					if((NI-1)%2 == 0)
					{
						ic = 0;i = 0;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i) = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i+1) = 0.5*(*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic)+*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic+1));
						ic = NIC-2;i = NI-3;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i) = 0.5*(*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic)+*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic-1));
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i+1) = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic);
						ic = 1;
						for(i=2;i<NI-3;i=i+2)
						{
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i) = 0.5*(*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic)+*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic-1));
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i+1) = 0.5*(*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic)+*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic+1));
							ic++;
						}
					}
					else
					{
						ic = 0;i = 0;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i) = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i+1) = 0.5*(*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic)+*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic+1));
						ic = NIC-2;i = NI-2;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i) = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic);
						ic = 1;
						for(i=2;i<NI-2;i=i+2)
						{
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i) = 0.5*(*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic)+*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic-1));
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i+1) = 0.5*(*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic)+*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*(NJC-2)+ic+1));
							ic++;
						}
					}
				break;
				case 3:
					k = 0;kc = 0;
					if((NI-1)%2 == 0)
					{
						ic = 0;i = 0;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i) = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i+1) = 0.5*(*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic)+*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic+1));
						ic = NIC-2;i = NI-3;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i) = 0.5*(*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic)+*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic-1));
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i+1) = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic);
						ic = 1;
						for(i=2;i<NI-3;i=i+2)
						{
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i) = 0.5*(*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic)+*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic-1));
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i+1) = 0.5*(*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic)+*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic+1));
							ic++;
						}
					}
					else
					{
						ic = 0;i = 0;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i) = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic);
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i+1) = 0.5*(*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic)+*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic+1));
						ic = NIC-2;i = NI-2;
						*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i) = *((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic);
						ic = 1;
						for(i=2;i<NI-2;i=i+2)
						{
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i) = 0.5*(*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic)+*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic-1));
							*((zone->grid+gdi-1)->DummyV+(NI-1)*(NJ-1)*k+(NI-1)*0+i+1) = 0.5*(*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic)+*((zone->grid+gdi)->DummyV+(NIC-1)*(NJC-1)*kc+(NIC-1)*0+ic+1));
							ic++;
						}
					}
				break;
			}
		}
	}
	return 0;
}
