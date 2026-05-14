#include"APKEvaluateVarCoefficients.h"
/**
The function  EvaluateVarCoefficients() loops over all the zones and computes the values of variable coefficients
i.e; Density, Viscosity, Thermal Conductivity etc, by calling their individual user defined functions.
This function is called only when the variable coefficients are choosen not to be constant.
*/
int EvaluateVarCoefficients(Base *base)
{
	int i,j,k,NI,NJ,NK,numZones,cntZ,*ptrSize;
	Zone *zone;
	double *variable;
	int totalcells,index;
	double *coordX,*coordY,*coordZ;
	double *ptrPoint1X,*ptrPoint1Y,*ptrPoint1Z,centroid[3];
	numZones= base->Number_of_zones;
	
	for(cntZ=0;cntZ<numZones;cntZ++)
	{
		zone = (base->CgnsZone+cntZ);
		ptrSize = zone->Size;
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		totalcells = (NI-1)*(NJ-1)*(NK-1);
		
		if(zone->Flowsol->flagRho == 1)
		{
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				index  = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
				variable = zone->Flowsol->Rho;

				coordX=(zone->CgnsGrid)->X;
				coordY=(zone->CgnsGrid)->Y; 
				coordZ=(zone->CgnsGrid)->Z;

				ptrPoint1X=coordX+NI*NJ*k+NI*j+i;
				ptrPoint1Y=coordY+NI*NJ*k+NI*j+i;
				ptrPoint1Z=coordZ+NI*NJ*k+NI*j+i;

				FindCellCentre(ptrSize,ptrPoint1X,ptrPoint1Y,ptrPoint1Z,centroid);
				*(variable+index) = userfuncDensity(zone,index,centroid[0],centroid[1],centroid[2],base);
			}
		}
		
		/*if(zone->Flowsol->flagelecond == 1)//rajendra
		{
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				index  = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
				variable = zone->Flowsol->Eleconductivity;

				coordX=(zone->CgnsGrid)->X;
				coordY=(zone->CgnsGrid)->Y; 
				coordZ=(zone->CgnsGrid)->Z;

				ptrPoint1X=coordX+NI*NJ*k+NI*j+i;
				ptrPoint1Y=coordY+NI*NJ*k+NI*j+i;
				ptrPoint1Z=coordZ+NI*NJ*k+NI*j+i;

				FindCellCentre(ptrSize,ptrPoint1X,ptrPoint1Y,ptrPoint1Z,centroid);
				*(variable+index) = userfuncEleconductivity(zone,index,centroid[0],centroid[1],centroid[2],base);
			}
		}
		
		if(zone->Flowsol->flagpermeability == 1)//rajendra
		{
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				index  = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
				variable = zone->Flowsol->Permeability;

				coordX=(zone->CgnsGrid)->X;
				coordY=(zone->CgnsGrid)->Y; 
				coordZ=(zone->CgnsGrid)->Z;

				ptrPoint1X=coordX+NI*NJ*k+NI*j+i;
				ptrPoint1Y=coordY+NI*NJ*k+NI*j+i;
				ptrPoint1Z=coordZ+NI*NJ*k+NI*j+i;

				FindCellCentre(ptrSize,ptrPoint1X,ptrPoint1Y,ptrPoint1Z,centroid);
				*(variable+index) = userfuncPermeability(zone,index,centroid[0],centroid[1],centroid[2],base);
			}
		}*/
		
		if(zone->Flowsol->flagViscoeff == 1)
		{
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				index  = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
				variable = zone->Flowsol->Viscoeff;

				coordX=(zone->CgnsGrid)->X;
				coordY=(zone->CgnsGrid)->Y; 
				coordZ=(zone->CgnsGrid)->Z;

				ptrPoint1X=coordX+NI*NJ*k+NI*j+i;
				ptrPoint1Y=coordY+NI*NJ*k+NI*j+i;
				ptrPoint1Z=coordZ+NI*NJ*k+NI*j+i;

				FindCellCentre(ptrSize,ptrPoint1X,ptrPoint1Y,ptrPoint1Z,centroid);

				*(variable+index) = userfuncViscosity(zone,index,centroid[0],centroid[1],centroid[2],base);
			}
		}
// 		printf("\nzone->Flowsol->flagRho == %d",zone->Flowsol->flagRho);
// 		printf("\nzone->Flowsol->flagViscoeff == %d",zone->Flowsol->flagViscoeff);
// 		printf("\nzone->Flowsol->flagSpecificHeat == %d",zone->Flowsol->flagSpecificHeat);
// 		printf("\nzone->Flowsol->flagScalcoeff[0] == %d",zone->Flowsol->flagScalcoeff[0]);
// 		printf("\nzone->Flowsol->flagScalcoeff[1] == %d",zone->Flowsol->flagScalcoeff[1]);
// 		printf("\nzone->Flowsol->flagScalcoeff[2] == %d",zone->Flowsol->flagScalcoeff[2]);
// 		printf("\nzone->Flowsol->flagScalcoeff[3] == %d",zone->Flowsol->flagScalcoeff[3]);
// 		printf("\nzone->Flowsol->flagScalcoeff[4] == %d",zone->Flowsol->flagScalcoeff[4]);
// 		printf("\nzone->Flowsol->flagScalcoeff[5] == %d",zone->Flowsol->flagScalcoeff[5]);
// 		printf("\nzone->Flowsol->flagScalcoeff[6] == %d",zone->Flowsol->flagScalcoeff[6]);
// 		printf("\nzone->Flowsol->flagScalcoeff[7] == %d",zone->Flowsol->flagScalcoeff[7]);
// 		printf("\nzone->Flowsol->flagScalcoeff[8] == %d",zone->Flowsol->flagScalcoeff[8]);
// 		printf("\nzone->Flowsol->flagScalcoeff[9] == %d",zone->Flowsol->flagScalcoeff[9]);
// 		printf("\nzone->Flowsol->flagChargeDiffCoeff == %d",zone->Flowsol->flagChargeDiffCoeff);
// 		printf("\nzone->Flowsol->flagCPdrop == %d",zone->Flowsol->flagCPdrop);
// 		printf("\nNum of turbulent scalars is %d",base->BaseParameters->Turb->numTurbScalars);
		
		if(zone->Flowsol->flagDiffcoeff == 1)
		{
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				index  = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
				variable = zone->Flowsol->Diffcoeff;

				coordX=(zone->CgnsGrid)->X;
				coordY=(zone->CgnsGrid)->Y; 
				coordZ=(zone->CgnsGrid)->Z;

				ptrPoint1X=coordX+NI*NJ*k+NI*j+i;
				ptrPoint1Y=coordY+NI*NJ*k+NI*j+i;
				ptrPoint1Z=coordZ+NI*NJ*k+NI*j+i;

				FindCellCentre(ptrSize,ptrPoint1X,ptrPoint1Y,ptrPoint1Z,centroid);

				*(variable+index) = userfuncThermalConductivity(zone,index,centroid[0],centroid[1],centroid[2],base);
			}
		}
		
		if(zone->Flowsol->flagSpecificHeat == 1)
		{
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				index  = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
				variable = zone->Flowsol->SpecificHeat;

				coordX=(zone->CgnsGrid)->X;
				coordY=(zone->CgnsGrid)->Y; 
				coordZ=(zone->CgnsGrid)->Z;

				ptrPoint1X=coordX+NI*NJ*k+NI*j+i;
				ptrPoint1Y=coordY+NI*NJ*k+NI*j+i;
				ptrPoint1Z=coordZ+NI*NJ*k+NI*j+i;

				FindCellCentre(ptrSize,ptrPoint1X,ptrPoint1Y,ptrPoint1Z,centroid);

				*(variable+index) = userfuncSpecificHeat(zone,index,centroid[0],centroid[1],centroid[2],base);
			}
		}
		
		if(zone->Flowsol->flagScalcoeff[0] == 1 && base->BaseParameters->Turb->numTurbScalars < 1)
		{
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				index  = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
				variable = zone->Flowsol->Scalcoeff[0];

				coordX=(zone->CgnsGrid)->X;
				coordY=(zone->CgnsGrid)->Y; 
				coordZ=(zone->CgnsGrid)->Z;

				ptrPoint1X=coordX+NI*NJ*k+NI*j+i;
				ptrPoint1Y=coordY+NI*NJ*k+NI*j+i;
				ptrPoint1Z=coordZ+NI*NJ*k+NI*j+i;

				FindCellCentre(ptrSize,ptrPoint1X,ptrPoint1Y,ptrPoint1Z,centroid);

				*(variable+index) = userfuncScalCoeff1(zone,index,centroid[0],centroid[1],centroid[2],base);
			}
		}
		
		if(zone->Flowsol->flagScalcoeff[1] == 1 && base->BaseParameters->Turb->numTurbScalars < 2)
		{
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				index  = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
				variable = zone->Flowsol->Scalcoeff[1];

				coordX=(zone->CgnsGrid)->X;
				coordY=(zone->CgnsGrid)->Y; 
				coordZ=(zone->CgnsGrid)->Z;

				ptrPoint1X=coordX+NI*NJ*k+NI*j+i;
				ptrPoint1Y=coordY+NI*NJ*k+NI*j+i;
				ptrPoint1Z=coordZ+NI*NJ*k+NI*j+i;

				FindCellCentre(ptrSize,ptrPoint1X,ptrPoint1Y,ptrPoint1Z,centroid);

				*(variable+index) = userfuncScalCoeff2(zone,index,centroid[0],centroid[1],centroid[2],base);
			}
		}
		if(zone->Flowsol->flagScalcoeff[2] == 1 && base->BaseParameters->Turb->numTurbScalars < 3)
		{
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				index  = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
				variable = zone->Flowsol->Scalcoeff[2];

				coordX=(zone->CgnsGrid)->X;
				coordY=(zone->CgnsGrid)->Y; 
				coordZ=(zone->CgnsGrid)->Z;

				ptrPoint1X=coordX+NI*NJ*k+NI*j+i;
				ptrPoint1Y=coordY+NI*NJ*k+NI*j+i;
				ptrPoint1Z=coordZ+NI*NJ*k+NI*j+i;

				FindCellCentre(ptrSize,ptrPoint1X,ptrPoint1Y,ptrPoint1Z,centroid);

				*(variable+index) = userfuncScalCoeff3(zone,index,centroid[0],centroid[1],centroid[2],base);
			}
		}
		
		if(zone->Flowsol->flagScalcoeff[3] == 1 && base->BaseParameters->Turb->numTurbScalars < 4)
		{
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				index  = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
				variable = zone->Flowsol->Scalcoeff[3];

				coordX=(zone->CgnsGrid)->X;
				coordY=(zone->CgnsGrid)->Y; 
				coordZ=(zone->CgnsGrid)->Z;

				ptrPoint1X=coordX+NI*NJ*k+NI*j+i;
				ptrPoint1Y=coordY+NI*NJ*k+NI*j+i;
				ptrPoint1Z=coordZ+NI*NJ*k+NI*j+i;

				FindCellCentre(ptrSize,ptrPoint1X,ptrPoint1Y,ptrPoint1Z,centroid);

				*(variable+index) = userfuncScalCoeff4(zone,index,centroid[0],centroid[1],centroid[2],base);
			}
		}
		
		if(zone->Flowsol->flagScalcoeff[4] == 1)
		{
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				index  = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
				variable = zone->Flowsol->Scalcoeff[4];

				coordX=(zone->CgnsGrid)->X;
				coordY=(zone->CgnsGrid)->Y; 
				coordZ=(zone->CgnsGrid)->Z;

				ptrPoint1X=coordX+NI*NJ*k+NI*j+i;
				ptrPoint1Y=coordY+NI*NJ*k+NI*j+i;
				ptrPoint1Z=coordZ+NI*NJ*k+NI*j+i;

				FindCellCentre(ptrSize,ptrPoint1X,ptrPoint1Y,ptrPoint1Z,centroid);

				*(variable+index) = userfuncScalCoeff5(zone,index,centroid[0],centroid[1],centroid[2],base);
			}
		}
		
		if(zone->Flowsol->flagScalcoeff[5] == 1)
		{
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				index  = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
				variable = zone->Flowsol->Scalcoeff[5];

				coordX=(zone->CgnsGrid)->X;
				coordY=(zone->CgnsGrid)->Y; 
				coordZ=(zone->CgnsGrid)->Z;

				ptrPoint1X=coordX+NI*NJ*k+NI*j+i;
				ptrPoint1Y=coordY+NI*NJ*k+NI*j+i;
				ptrPoint1Z=coordZ+NI*NJ*k+NI*j+i;

				FindCellCentre(ptrSize,ptrPoint1X,ptrPoint1Y,ptrPoint1Z,centroid);

				*(variable+index) = userfuncScalCoeff6(zone,index,centroid[0],centroid[1],centroid[2],base);
			}
		}
		
		if(zone->Flowsol->flagScalcoeff[6] == 1)
		{
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				index  = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
				variable = zone->Flowsol->Scalcoeff[6];

				coordX=(zone->CgnsGrid)->X;
				coordY=(zone->CgnsGrid)->Y; 
				coordZ=(zone->CgnsGrid)->Z;

				ptrPoint1X=coordX+NI*NJ*k+NI*j+i;
				ptrPoint1Y=coordY+NI*NJ*k+NI*j+i;
				ptrPoint1Z=coordZ+NI*NJ*k+NI*j+i;

				FindCellCentre(ptrSize,ptrPoint1X,ptrPoint1Y,ptrPoint1Z,centroid);

				*(variable+index) = userfuncScalCoeff7(zone,index,centroid[0],centroid[1],centroid[2],base);
			}
		}
		
		if(zone->Flowsol->flagScalcoeff[7] == 1)
		{
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				index  = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
				variable = zone->Flowsol->Scalcoeff[7];

				coordX=(zone->CgnsGrid)->X;
				coordY=(zone->CgnsGrid)->Y; 
				coordZ=(zone->CgnsGrid)->Z;

				ptrPoint1X=coordX+NI*NJ*k+NI*j+i;
				ptrPoint1Y=coordY+NI*NJ*k+NI*j+i;
				ptrPoint1Z=coordZ+NI*NJ*k+NI*j+i;

				FindCellCentre(ptrSize,ptrPoint1X,ptrPoint1Y,ptrPoint1Z,centroid);

				*(variable+index) = userfuncScalCoeff8(zone,index,centroid[0],centroid[1],centroid[2],base);
			}
		}
		
		if(zone->Flowsol->flagScalcoeff[8] == 1)
		{
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				index  = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
				variable = zone->Flowsol->Scalcoeff[8];

				coordX=(zone->CgnsGrid)->X;
				coordY=(zone->CgnsGrid)->Y; 
				coordZ=(zone->CgnsGrid)->Z;

				ptrPoint1X=coordX+NI*NJ*k+NI*j+i;
				ptrPoint1Y=coordY+NI*NJ*k+NI*j+i;
				ptrPoint1Z=coordZ+NI*NJ*k+NI*j+i;

				FindCellCentre(ptrSize,ptrPoint1X,ptrPoint1Y,ptrPoint1Z,centroid);

				*(variable+index) = userfuncScalCoeff9(zone,index,centroid[0],centroid[1],centroid[2],base);
			}
		}
		
		if(zone->Flowsol->flagScalcoeff[9] == 1)
		{
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				index  = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
				variable = zone->Flowsol->Scalcoeff[9];

				coordX=(zone->CgnsGrid)->X;
				coordY=(zone->CgnsGrid)->Y; 
				coordZ=(zone->CgnsGrid)->Z;

				ptrPoint1X=coordX+NI*NJ*k+NI*j+i;
				ptrPoint1Y=coordY+NI*NJ*k+NI*j+i;
				ptrPoint1Z=coordZ+NI*NJ*k+NI*j+i;

				FindCellCentre(ptrSize,ptrPoint1X,ptrPoint1Y,ptrPoint1Z,centroid);

				*(variable+index) = userfuncScalCoeff10(zone,index,centroid[0],centroid[1],centroid[2],base);
			}
		}
		
		if(zone->Flowsol->flagChargeDiffCoeff == 1)
		{
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				index  = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
				variable = zone->Flowsol->ChargeDiffcoeff;

				coordX=(zone->CgnsGrid)->X;
				coordY=(zone->CgnsGrid)->Y; 
				coordZ=(zone->CgnsGrid)->Z;

				ptrPoint1X=coordX+NI*NJ*k+NI*j+i;
				ptrPoint1Y=coordY+NI*NJ*k+NI*j+i;
				ptrPoint1Z=coordZ+NI*NJ*k+NI*j+i;

				FindCellCentre(ptrSize,ptrPoint1X,ptrPoint1Y,ptrPoint1Z,centroid);

				*(variable+index) = userfuncChargeDiffusion(zone,index,centroid[0],centroid[1],centroid[2],base);
			}
		}
		if(zone->Flowsol->flagCPdrop == 1)//shash
		{
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				index  = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
				variable = zone->Flowsol->CPdrop;

				coordX=(zone->CgnsGrid)->X;
				coordY=(zone->CgnsGrid)->Y; 
				coordZ=(zone->CgnsGrid)->Z;

				ptrPoint1X=coordX+NI*NJ*k+NI*j+i;
				ptrPoint1Y=coordY+NI*NJ*k+NI*j+i;
				ptrPoint1Z=coordZ+NI*NJ*k+NI*j+i;

				FindCellCentre(ptrSize,ptrPoint1X,ptrPoint1Y,ptrPoint1Z,centroid);

				*(variable+index) = userfuncCPdrop(zone,index,centroid[0],centroid[1],centroid[2],base);
				
			}
		}
	}
	
return 0;
}

		
