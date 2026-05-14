#include "APKFillVariable2Dummy2d.h"

int FillVariable2Dummy2d(int variableTag,Base *base,Zone *zone,int gdi)
{
	int iface;
	double *Variable,*DV;
	char tag,name[33];
	void *ptr;
	int i,j,k,NI,NJ,NK;
	
	Variable = (zone->grid+gdi)->V;
	DV = (zone->grid+gdi)->DummyV;
	
	NI = *((zone->grid+gdi)->GSize+0);
	NJ = *((zone->grid+gdi)->GSize+1);
	NK = *((zone->grid+gdi)->GSize+2);
	
	for(iface=0;iface<4;iface++)
	{
		ptr = (void *)FindInterfaceBCPtrMod(zone,iface+1,name,&tag);
		if(tag=='i')
		{
			switch(iface)
			{
				case 0:
					k = 0;
					for(j=0;j<NJ-1;j++)
					{
						*(DV+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2) = *(Variable+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2);
					}
				break;
				case 1:
					k = 0;
					for(j=0;j<NJ-1;j++)
					{
						*(DV+(NI-1)*(NJ-1)*k+(NI-1)*j+0) = *(Variable+(NI-1)*(NJ-1)*k+(NI-1)*j+0);
					}
				break;
				case 2:
					k = 0;
					for(i=0;i<NI-1;i++)
					{
						*(DV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i) = *(Variable+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i);
					}
				break;
				case 3:
					k = 0;
					for(i=0;i<NI-1;i++)
					{
						*(DV+(NI-1)*(NJ-1)*k+(NI-1)*0+i) = *(Variable+(NI-1)*(NJ-1)*k+(NI-1)*0+i);
					}
				break;
			}
		}
	}
	return 0;
}
