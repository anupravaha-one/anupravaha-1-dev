#include "APKFindInterfaceBCPtr.h"
/**
See the function 'FindInterfaceBcptrMod() for description. Currently that function is implemented
*/
Interface *FindInterfaceBCPtr(Zone *zone,int dir,char *zoneName,char *tag)
{
	Interface *facePtr;
	*tag='n';
	//printf("name of current zone= %s \n",zone->zonename);
	int cntBC,cntFace;
	//printf("number of interfaces = %d\n",	zone->Number_of_Interfaces);
	//printf("number of Boundaries = %d\n",	zone->Number_of_BC_Patches);
	for(cntFace=0;cntFace<(zone->Number_of_Interfaces);cntFace++)
	{
		if(dir==((zone->Iface+cntFace)->Direction))
		{
			*tag='i';
			strcpy(zoneName,((zone->Iface+cntFace)->Adjoining_Zone_Name));
			//printf("Zone name inside %s \n",zoneName);
			facePtr=(zone->Iface+cntFace);
			//printf("name of the interface= %s\n",facePtr->Interface_Name);
			//printf("dir =%d, face pointer = %u\n",dir,facePtr);
			break;
		}
	}
	if(cntFace==(zone->Number_of_Interfaces))
	{
		for(cntBC=0;cntBC < (zone->Number_of_BC_Patches); cntBC++)
			if(dir==((zone->ZoneBc+cntBC)->Direction))
			{
				//printf("i ve found a boundary\n");
				*tag='b';
				zoneName="boundary";
				/*facePtr=(zone->ZoneBc+cntBC);*/
				//printf("#######################\n");
				break;
			}
	}
	//printf("%c\n",*tag);
	//printf("facePtr= %u\n",facePtr);
	if(*tag=='n')
		printf("Direction din't match with either interface or boundary\n");
	return(facePtr);
}
