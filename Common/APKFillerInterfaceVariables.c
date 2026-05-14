#include "APKFillerInterfaceVariables.h"

int FillerInterfaceVariables(int variableTag,Interface *facePtr,Zone *ptrZone,int iface)
{
	double *slabPointer,*domainPointer;
	int i,totalCells;
	
	slabPointer=FindSlabPtrIface(facePtr,variableTag);
	domainPointer=MakeDomainArray(ptrZone,variableTag);
	
	totalCells=facePtr->ZSlab->Total_Cells;
	
	for(i=0;i<totalCells;i++)
	{
		*(slabPointer+i)= *(domainPointer+(*(facePtr->ZSlab->Target_Indices+i)));
	}
	return 0;
}

