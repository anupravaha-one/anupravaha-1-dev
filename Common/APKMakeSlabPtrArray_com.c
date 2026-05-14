#include "APKMakeSlabPtrArray_com.h"

/**
This function is originally MakeSlabPtrArray.com of Version 3.0.  During profiling it has been observed that this original function is consuming 40% of the total run time (including the other function calls from the function).  This original function finds out the slab pointers of all the variables interfaces/boundaries.  This search for slabs was consuming most of the run time.  In the original function, this search is replaced by directly pointing the pointers to an pointer array which is stored in the data structure.  But at the starting of the program, these array are to be filled with correct addresses.  

Refer Profiling and Optimization manual for further details.

This function is added by BHANS.
**/


int MakeSlabPtrArray_com(Zone *zone,double *ptrArray[6],int vartag)
{
	int cntDir;
	void *ptr;
	//char zoneName[33];
	//char tag;
	int tag;
	Interface *facePtr;
	ZoneBoundary *BCPtr;
	
	for(cntDir=0;cntDir<6;cntDir++)
	{		
		//ptr=(void *)FindInterfaceBCPtrMod(zone,cntDir+1,zoneName,&tag);
		ptr = zone->ZONE_FACES[cntDir];
		tag = zone->BC_IFACE_TAG[cntDir];
		//if(tag=='i')
		if(tag==2)
		{
			facePtr=(Interface *)ptr;
			ptrArray[cntDir]=FindSlabPtrIface(facePtr,vartag);
		}
		else if(tag==1)
		//if(tag=='b')
		{	
			BCPtr=(ZoneBoundary *)ptr;
			ptrArray[cntDir]= FindSlabPtrBC(BCPtr,vartag);
		}
		else
		{
			printf("Proper tag is not present in APKMakeSlabPtrArray_com.c\n\n");
			exit(1);
		}
	}
	return 0;
}
