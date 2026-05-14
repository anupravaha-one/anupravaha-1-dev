#include "APKMakeDiffusionCoeffSlabPtrArray.h"
/**
This function is used to make an array of pointers of diffusion coefficients present in all the boundaries and interfaces
corresponding to a particular zone and corrsponding to a particular variable tag.
This functions takes the zone pointer and variable tag, and loops over all the zone boundaries and interfaces, and calls the
functions FindDiffusionCoeffSlabPtrBC() and FindDiffusionCoeffSlabPtrIface() respectively, which return the slab pointers
of diffusion coefficients corresponding to that particular variable tag. All these slab pointers are stored in an array.
*/

//NOTE: BHANESH: REVIEW this function & UPLOAD it for use as this is efficient than the existing function 
//<<< Narsimha - 22nd May, 2008
// int MakeDiffusionCoeffSlabPtrArray(Zone *zone,double *ptrArray[6],int vartag)
// {
// 	int index;
// 	Interface *facePtr;
// 	ZoneBoundary *BCPtr;
// 	
// 	for(index=0; index < zone->Number_of_Interfaces; index++)
// 	{
// 		facePtr = zone->Iface+index;
// 		ptrArray[facePtr->Direction-1] = FindDiffusionCoeffSlabPtrIface(zone,facePtr,vartag);
// 	}
// 	
// 	for(index=0; index < zone->Number_of_BC_Patches; index++)
// 	{
// 		BCPtr = zone->ZoneBc+index;
// 		ptrArray[BCPtr->Direction-1] = FindDiffusionCoeffSlabPtrBC(zone,BCPtr,vartag);
// 	}
// 	
// 	return 0;
// }

int MakeDiffusionCoeffSlabPtrArray(Zone *zone,double *ptrArray[6],int vartag)
{
	int cntDir;
	void *ptr;
	char tag,zoneName[33];
	Interface *facePtr;
	ZoneBoundary *BCPtr;
	
	for(cntDir=0;cntDir<6;cntDir++)
	{
		ptr=(void *)FindInterfaceBCPtrMod(zone,cntDir+1,zoneName,&tag);
		if(tag=='i')
		{
			facePtr=(Interface *)ptr;
			ptrArray[cntDir]=FindDiffusionCoeffSlabPtrIface(zone,facePtr,vartag);
		}
		else
		if(tag=='b')
		{
			BCPtr=(ZoneBoundary *)ptr;
			ptrArray[cntDir]= FindDiffusionCoeffSlabPtrBC(zone,BCPtr,vartag);
		}
		else
		{
			printf("Exiting from MakeDiffusionCoeffSlabPtrArray.c\n");
			exit(1);
		}
	}
	return 0;
}

