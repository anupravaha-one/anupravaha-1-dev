#include "APKHomogeneousNeumannIC.h"
/**This function is used to apply Homogeneous Neumann BC at the SOLID-FLUID zone interface having ORTHOGONAL GRID ONLY.  
In general this function can be called for any variable which needs homogeneous Neumann BC to be applied at the
SOLID-FLUID interface having ORTHOGONAL GRID ONLY.
 */
int HomogeneousNeumannIC(Zone *zone,Interface *facePtr,int variableTag)
{
	int i,j,k,indexP,targetIndices=0;
	int ni,nj,nk;
	double *slabPointer,*domainPointer;
	
	ni=*(zone->Size+0);
	nj=*(zone->Size+1);
	nk=*(zone->Size+2);
	
	slabPointer=FindSlabPtrIface(facePtr,variableTag);
	domainPointer=MakeDomainArray(zone,variableTag);
	
	switch (facePtr->Direction)
	{
		case 1://Imax face
			i=ni-2;
			for (k=0;k<nk-1;k++) 
			for (j=0;j<nj-1;j++)
			{
				indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
				*(slabPointer+targetIndices)=*(domainPointer+indexP);
				targetIndices++;
			}
		break;

		case 2://Imin face
			i=0;
			for (k=0;k<nk-1;k++) 
			for (j=0;j<nj-1;j++)
			{
				indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
				*(slabPointer+targetIndices)=*(domainPointer+indexP);
				targetIndices++;
			}
		break;

		case 3://Jmax face
			j=nj-2;
			for (k=0;k<nk-1;k++) 
			for (i=0;i<ni-1;i++)
			{
				indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
				*(slabPointer+targetIndices)=*(domainPointer+indexP);
				targetIndices++;
			}
		break;

		case 4://Jmin face
			j=0;
			for (k=0;k<nk-1;k++) 
			for (i=0;i<ni-1;i++)
			{
				indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
				*(slabPointer+targetIndices)=*(domainPointer+indexP);
				targetIndices++;
			}
		break;

		case 5://Kmax face
			k=nk-2;
			for (j=0;j<nj-1;j++) 
			for (i=0;i<ni-1;i++)
			{
				indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
				*(slabPointer+targetIndices)=*(domainPointer+indexP);
				targetIndices++;
			}
		break;

		case 6://Kmin face
			k=0;
			for (j=0;j<nj-1;j++) 
			for (i=0;i<ni-1;i++)
			{
				indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
				*(slabPointer+targetIndices)=*(domainPointer+indexP);
				targetIndices++;
			}
		break;
	}
	
	facePtr->ConjugateRms = 0.0;
	return (0);
}

