#include "APKFindInterfaceBCPtrMod.h"
/**
This function, FindInterfaceBCPtrMod() is used to find interface and zoneboundary pointer and the related tag i.e either 'i' or 'b'.
This also finds adjacent zonename if zone face is the interface. These are required to get the slab values.
*/

void *FindInterfaceBCPtrMod(Zone *zone,int dir,char *zoneName,char *tag)
{
	void *facePtr;
	*tag='n';

	int cntBC,cntFace;
	

	for(cntFace=0;cntFace<(zone->Number_of_Interfaces);cntFace++)
	{		
		if(dir==((zone->Iface+cntFace)->Direction))
		{	
			*tag='i';
			strcpy(zoneName,((zone->Iface+cntFace)->Adjoining_Zone_Name));
		
			facePtr=(zone->Iface+cntFace);
			//printf("name of the interface= %s\n",facePtr->Interface_Name);
			//printf("dir =%d, face pointer = %u\n",dir,facePtr);
			break;
		}
	}
	if(cntFace==(zone->Number_of_Interfaces))
	{
		for(cntBC=0;cntBC < (zone->Number_of_BC_Patches); cntBC++)
		{
			if(dir==((zone->ZoneBc+cntBC)->Direction))
			{
				//printf("i ve found a boundary\n");
				*tag='b';
				zoneName="boundary";
				facePtr=(void *)(zone->ZoneBc+cntBC);
				//printf("#######################\n");
				break;
			}
		}
	}
	
	//printf("%c\n",*tag);
	//printf("facePtr= %u\n",facePtr);
	if(*tag=='n');
	//	printf("Direction din't match with either interface or boundary\n");
	return(facePtr);
}
