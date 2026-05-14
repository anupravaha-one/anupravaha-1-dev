#include "APKsolverLMMHD.h"
int solverLMMHD(Base *base)
{
	Zone *zone;
	int g,NI,NJ,NK,*ptrSize;//rajendra
	double BstarResidual[3];
	int iZone,numZones;
	numZones=base->Number_of_zones;
	//totalBaseCells = base->TotalCellFluid;
	for(iZone=0;iZone<numZones;iZone++)
	{
		zone=base->CgnsZone+iZone;
		g = 0;
		ptrSize = (zone->grid+g)->GSize;
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		(zone->grid+g)->GRHS=(double *)malloc((NI-1)*(NJ-1)*(NK-1)*sizeof(double));
	}
	
	if(base->BaseParameters->flag2d==1)
	{
		if(base->flaginduced==1)
		{
			
 			savemagneticflux(base);
			CallCalcCoefficientsForVariables(base,0,56);
			Conv_Coeff_Calculation(base,0,56);
			CDcoeffCalculation(base,56);
			
			Induced_Source_Calculation(base);
			
				CalcSource_opti(base,56);

				BstarResidual[0] = CalcVariable(base,56);
			
				CalcSource_opti(base,57);

				BstarResidual[1] = CalcVariable(base,57);
			
				BstarResidual[2] = 0.0;
				
			Magneticflux2d(base);
			
			CalcPsiSource(base);
			CalcConjugate(base,59);
			 			
 			UpdateinducedB2d(base);
 			
 			CalcNewMagneticflux2d(base);

	            }
		else
		{
			CalcUcrossB2d(base);
			PotentialMHDpoissonRHS2d(base);
			CalcSource_opti(base,17);
			PotentialMHDIter(base);
		     
		}
		
	}

	else
	{
		if(base->flaginduced==1)
		{
			savemagneticflux(base);
			CallCalcCoefficientsForVariables(base,0,56);
			Conv_Coeff_Calculation(base,0,56);
			CDcoeffCalculation(base,56);

			Induced_Source_Calculation(base);

			CalcSource_opti(base,56);
			BstarResidual[0] = CalcVariable(base,56);


			CalcSource_opti(base,57);	
			BstarResidual[1] = CalcVariable(base,57);

			CalcSource_opti(base,58);
			BstarResidual[2] = CalcVariable(base,58);

			Magneticflux(base);
 			
			CalcPsiSource(base);

		       	CalcConjugate(base,59);
                        			
			UpdateinducedB(base);
			
			CalcNewMagneticflux(base);

                     }
		else
		{
			CalcUcrossB(base);

			PotentialMHDpoissonRHS(base);
			CalcSource_opti(base,17);
			PotentialMHDIter(base);
		}
		
	}
// 	
	//RefreshSurfaceCells(59,base);
	
	for(iZone=0;iZone<numZones;iZone++)
	{
		zone=base->CgnsZone+iZone;
		free((zone->grid)->GRHS);
	}
	return 0;
}

