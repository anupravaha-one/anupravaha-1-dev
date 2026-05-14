#include "APKFindgradientPhiMHD.h"
/**
The function FindgradientPhiMHD() calls the function CalcGradVar() to compute the gradients of PotentialMHD, which arise in source terms of the momentum equations in case of Liquidmetal MHD flows.
Please refer Eqn() - Eqn() in Thesis Report of Vetcha Naveen.
*/

double FindgradientPhiMHD(Base *base,Zone *zone,int index[3],double Dphi[3])
{
	//double gradphi;
  double *gradPhiSlab[6], *phi,*volume,*SurfArea,vol,*slabVolume[6];
  int *ptrSize;

                            	phi=MakeDomainArray(zone,17);
                         	MakeSlabPtrArray(zone,gradPhiSlab,17);
	                        volume=MakeDomainArray(zone,5);
                          	SurfArea=zone->ZGeom->SurfaceArea;
                               	ptrSize= zone->Size;
	
				MakeSlabPtrArray(zone,slabVolume,5);
                                Dphi[0]=CalcPressureGrad(phi,gradPhiSlab,volume,slabVolume,SurfArea,index,ptrSize,1);
                                Dphi[1]=CalcPressureGrad(phi,gradPhiSlab,volume,slabVolume,SurfArea,index,ptrSize,2);
                                Dphi[2]=CalcPressureGrad(phi,gradPhiSlab,volume,slabVolume,SurfArea,index,ptrSize,3);

                                Dphi[0]= -Dphi[0];
                                Dphi[1]= -Dphi[1];
                                Dphi[2]= -Dphi[2];



				//	CalcGradVar(base,zone,phi,gradPhiSlab,index,Dphi);
	
	return 0 ;
}
