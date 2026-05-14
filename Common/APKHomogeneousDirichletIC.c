#include "APKHomogeneousDirichletIC.h"
/**This function is used to apply Homogeneous Dirichlet condition at the SOLID-FLUID zone interface.
In general this function can be called for any variable which needs Homogeneous Dirichlet condition to be applied at the
SOLID-FLUID interface.
 */

int HomogeneousDirichletIC(int variableTag,Interface *facePtr)
{
	int targetIndices,iface,totalCells,NIcell,NJcell,NKcell,*ptrRange;
	double *slabPointer;
	
	iface = facePtr->Direction;
	ptrRange = facePtr->Point_Range;
	
	NIcell = *(ptrRange+3) - *(ptrRange+0);
	NJcell = *(ptrRange+4) - *(ptrRange+1);
	NKcell = *(ptrRange+5) - *(ptrRange+2);
	
	slabPointer=FindSlabPtrIface(facePtr,variableTag);
	
	switch(iface)
	{
		case 1://Imax face
		case 2://Imin face
			totalCells = NJcell * NKcell;
			break;
		case 3://Jmax face
		case 4://Jmin face
			totalCells = NIcell * NKcell;
			break;
		case 5://Kmax face
		case 6://Kmin face
			totalCells = NIcell * NJcell;
			break;
	}
	
	for(targetIndices=0; targetIndices<totalCells; targetIndices++)
	{
		*(slabPointer+targetIndices) = 0.0;
	}
	
	facePtr->ConjugateRms = 0.0;//<<<< Narsimha
	return 0;
}


