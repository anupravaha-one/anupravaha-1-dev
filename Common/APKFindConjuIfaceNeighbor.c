#include "APKFindConjuIfaceNeighbor.h"

int FindConjuIfaceNeighbor(int variableTag,Zone *zone,Zone *zoneNeighbor)
{
	char tag,name[33];
	int iface,flag;
	void *ptr;
	Interface *facePtr;
	
	for(iface=0;iface < 6;iface++)
	{
		ptr = (void *)FindInterfaceBCPtrMod(zoneNeighbor,iface+1,name,&tag);
		if(tag=='i')
		{
			facePtr=(Interface *)ptr;
			//printf("Zname=%s\tname=%s\n",((zone)->zonename),name);
			flag=strcmp(((zone)->zonename),name);
			if(flag==0 )
			{
				break;
			}
		}
	}
	return (iface);
}
