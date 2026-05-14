#include "APKfind_zone_face_info.h"

/**
In Data structure the following variables are added to make the interface data transfer easier.
int BC_IFACE_TAG[6];
void *ZONE_FACES[6]; 
int ZONE_NEIGHBOUR\_NUM[6];

The BC_IFACE_TAG indicates whether a particular face of a zone is at interface or boundary.  Tag=1 means boundary and Tag=2 means it's an interface.  ZONE_FACES stores the void pointers of the faces of the zone (for both BC and IFACE).
ZONE_NEIGHBOUR_NUM indicates the value of the neighboring zone number(index starting from zero).  Stores -1 for boundaries.  This information can be retrieved wherever necessary.  We can readily use this stored information instead of finding out every time.

APKfind_zone_face_info.c is added in Common folder.  This function calculates the zone face information (above added variables) for each and every zone and stores it in Data Structure.  In function Main.c this added function is called.

This function is added by BHANS.

Refer Profiling and Optimization manual for further details.
**/

int find_zone_face_info(Base *base)
{

	int cntBC,cntFace,numZones,tag,cntZone,dir;
	Zone *zone;	
	
	for(numZones=0;numZones<base->Number_of_zones;numZones++)
	{
		zone = base->CgnsZone+numZones;
		for(dir=1;dir<=6;dir++)
		{
			for(cntFace=0;cntFace<(zone->Number_of_Interfaces);cntFace++)
			{
				if(dir==((zone->Iface+cntFace)->Direction))
					break;
			}
			if(cntFace==(zone->Number_of_Interfaces))
			{
				for(cntBC=0;cntBC < (zone->Number_of_BC_Patches); cntBC++)
				{		if(dir==((zone->ZoneBc+cntBC)->Direction))
							break;
				}
				if(cntBC==zone->Number_of_BC_Patches)
				{
					printf("Tagging is not proper in ZONE Interface information\n");
					printf("Exiting in APKfind_zone_face_info.c\n");
					exit(123);
				}
				else
					zone->BC_IFACE_TAG[dir-1]=1;
			}
			else
				zone->BC_IFACE_TAG[dir-1]=2;
	
			if(zone->BC_IFACE_TAG[dir-1]==2)//In case of Interface
			{	
				for(cntZone=0;cntZone < (base->Number_of_zones);cntZone++)
				{
					tag=strcmp(((zone->Iface+cntFace)->Adjoining_Zone_Name),((base->CgnsZone+cntZone)->zonename));
					if(tag==0)
					break;
				}
				zone->ZONE_NEIGHBOUR_NUM[dir-1]=cntZone;
				zone->ZONE_FACES[dir-1]=(void *)(zone->Iface+cntFace);
				
			}
			else if (zone->BC_IFACE_TAG[dir-1]==1)
			{
				zone->ZONE_NEIGHBOUR_NUM[dir-1]=-1;
				zone->ZONE_FACES[dir-1]=(void *)(zone->ZoneBc+cntBC);
			}
			else
			{
				printf("Wrong tag stored in APKfind_zone_face_info.c\n");
				printf("Exiting.....\n\n");
				exit(123);
			}
		}
	}
// 	int dum=0;
// 	for(numZones=0;numZones<base->Number_of_zones;numZones++)
// 	{
// 		zone = base->CgnsZone+numZones;
// 		for(dir=0;dir<6;dir++)
// 		printf("This directions tag is %d\tNeibour zone number is %d\tvoid ptr=%d %d\n",zone->BC_IFACE_TAG[dir],zone->ZONE_NEIGHBOUR_NUM[dir],zone->ZONE_FACES[dir],dum++);
// 	}
// 	exit(12);

	return 0;
}


