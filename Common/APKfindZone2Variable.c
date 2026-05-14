#include "APKfindZone2Variable.h"
/**
This function, findZone2Variable() is for transfering the values of a variable from the neighbouring zone sharing face cells into the
interface zone slabs. These values for a variable are required when the solver is solving at the interface.
*/

int findZone2Variable(int variableTag,Interface *facePtr,Zone *zone,int index2[3],int indexofSlabs)
{
	int k,j,i,NI,NJ,NK;
	int *sizePtrNeighbor;
	double *slabPointer,*slabPointerMass,*domainPointer,*domainPointerMass;
	
	sizePtrNeighbor=zone->Size;
	i=index2[0]-1;
	j=index2[1]-1;
	k=index2[2]-1;
	
	NI=*(sizePtrNeighbor+0);
	NJ=*(sizePtrNeighbor+1);
	NK=*(sizePtrNeighbor+2);
	
// 	Groundwater Sanz..
	/*if((zone->Difftag==1 && (zone->Flowsol->flagTestSaturation!=4))||(zone->Difftag==0 && (zone->Flowsol->flagTestSaturation!=4)))
	{  
		if(variableTag >= 700 && variableTag<=709)//For Scalars
		{ 
			slabPointer=FindSlabPtrIface(facePtr,variableTag);
			domainPointer=MakeDomainArray(zone,variableTag);
	   		for(c=0;c<9;c++)
	  			*(slabPointer+9*indexofSlabs+c)= *(domainPointer+9*((NI-1)*(NJ-1)*k+(NI-1)*j+i)+c);
	  	}
	   	else
	   	{
			slabPointer=FindSlabPtrIface(facePtr,variableTag);
			domainPointer=MakeDomainArray(zone,variableTag);
			*(slabPointer+indexofSlabs)= *(domainPointer+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
	  	}
	}  
	else */
	{
		if((variableTag==1)||(variableTag==2)||(variableTag==3)||(variableTag==31)||(variableTag==32)||(variableTag==33)||(variableTag==62)||(variableTag==63)||(variableTag==64))
		{
			slabPointer=FindSlabPtrIface(facePtr,variableTag);
			slabPointerMass=FindSlabPtrIface(facePtr,variableTag*100);
			domainPointer=MakeDomainArray(zone,variableTag);
			domainPointerMass=MakeDomainArray(zone,variableTag*100);
		
			*(slabPointer+indexofSlabs)= *(domainPointer+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			*(slabPointerMass+indexofSlabs)= *(domainPointerMass+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		}
		else
		{
			slabPointer=FindSlabPtrIface(facePtr,variableTag);
			domainPointer=MakeDomainArray(zone,variableTag);
			*(slabPointer+indexofSlabs)= *(domainPointer+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		}
	}
	return 0;
}
