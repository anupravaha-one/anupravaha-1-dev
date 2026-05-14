#include "APKFindZonePtr.h"
/**
This function is used to find the pointer of a Zone, given a zone name.
 
The name of the 
Zone and the Base pointer are the input arguments. The function loops over all the Zone and compares 
the their names with the input name. The pointer of the Zone, whose name matches with the input name 
is return as output.
*/
Zone *FindZonePtr(char nameofZone[33],Base *base)/* checked working properly*/
{
	//printf("name recieved= %s and base recieved is %u\n",nameofZone,base);
	int cntZone,tag=1;
	Zone *ptr;
	for(cntZone=0;cntZone < (base->Number_of_zones);cntZone++)
	{
		tag=strcmp(nameofZone,((base->CgnsZone+cntZone)->zonename));
		if(tag==0)
		break;
	}
	ptr=base->CgnsZone+cntZone;
	//printf("Name of the Zone = %s \n",ptr->zonename);
	/*printf("Name of the Zone = %s \n\n",nameofZone);*/
	//printf("Zone ptr = %u\n",ptr);
	return(ptr);
}
