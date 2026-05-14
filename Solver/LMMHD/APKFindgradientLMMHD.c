#include "APKFindgradientLMMHD.h"


double FindgradientLMMHD(Base *base,Zone *zone,int index[3],double gradx[3],double grady[3],double gradz[3])
{
	//double gradphi;
	double *gradBxSlab[6],*gradBySlab[6],*gradBzSlab[6];
	double *B_x,* B_y,* B_z;
	//int indexP,NI,NJ,NK,*ptrSize,k,j,i;

// 	ptrSize= zone->Size;
// 	NI=*(ptrSize+0);
// 	NJ=*(ptrSize+1);
// 	NK=*(ptrSize+2);
// 	k=index[0],j=index[1],i=index[2];
// 	indexP=(NJ-1)*(NI-1)*k+(NI-1)*j+i;

	B_x = (zone->PrevFlowsol->TotalMfieldX);
	B_y = (zone->PrevFlowsol->TotalMfieldY);
	B_z = (zone->PrevFlowsol->TotalMfieldZ);

	MakeSlabPtrArray(zone,gradBxSlab,56);
	MakeSlabPtrArray(zone,gradBySlab,57);
	MakeSlabPtrArray(zone,gradBzSlab,58);

	CalcGradVar(base,zone,B_x,gradBxSlab,index,gradx);
	CalcGradVar(base,zone,B_y,gradBySlab,index,grady);
	CalcGradVar(base,zone,B_z,gradBzSlab,index,gradz);
	
	return 0 ;
}
