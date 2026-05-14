#include "APKFillVariable2Dummy.h"

int FillVariable2Dummy(int variableTag,Base *base,Zone *zone,int gdi)
{
	int iface;
	char tag,name[33];
	void *ptr;
	int i,j,k,NI,NJ,NK;
	double *Variable,*DV;
	
	Variable = (zone->grid+gdi)->V;
	DV = (zone->grid+gdi)->DummyV;
	
	NI = *((zone->grid+gdi)->GSize+0);
	NJ = *((zone->grid+gdi)->GSize+1);
	NK = *((zone->grid+gdi)->GSize+2);
	
	for(iface=0;iface<6;iface++)
	{
		ptr = (void *)FindInterfaceBCPtrMod(zone,iface+1,name,&tag);
		if(tag=='i')
		{
			switch(iface)
			{
				case 0:
					for(k=0;k<NK-1;k++)
					for(j=0;j<NJ-1;j++)
					{
						*(DV+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2) = *(Variable+(NI-1)*(NJ-1)*k+(NI-1)*j+NI-2);
					}
				break;
				case 1:
					for(k=0;k<NK-1;k++)
					for(j=0;j<NJ-1;j++)
					{
						*(DV+(NI-1)*(NJ-1)*k+(NI-1)*j+0) = *(Variable+(NI-1)*(NJ-1)*k+(NI-1)*j+0);
					}
				break;
				case 2:
					for(k=0;k<NK-1;k++)
					for(i=0;i<NI-1;i++)
					{
						*(DV+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i) = *(Variable+(NI-1)*(NJ-1)*k+(NI-1)*(NJ-2)+i);
					}
				break;
				case 3:
					for(k=0;k<NK-1;k++)
					for(i=0;i<NI-1;i++)
					{
						*(DV+(NI-1)*(NJ-1)*k+(NI-1)*0+i) = *(Variable+(NI-1)*(NJ-1)*k+(NI-1)*0+i);
					}
				break;
				case 4:
					for(j=0;j<NJ-1;j++)
					for(i=0;i<NI-1;i++)
					{
						*(DV+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i) = *(Variable+(NI-1)*(NJ-1)*(NK-2)+(NI-1)*j+i);
					}
				break;
				case 5:
					for(j=0;j<NJ-1;j++)
					for(i=0;i<NI-1;i++)
					{
						*(DV+(NI-1)*(NJ-1)*0+(NI-1)*j+i) = *(Variable+(NI-1)*(NJ-1)*0+(NI-1)*j+i);
					}
				break;
			}
		}
	}
	return 0;
}
