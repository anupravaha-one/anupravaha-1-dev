#include "APKTURBCalcTurbconductivity.h"
/**This fuction calculates the turbulent thermal conductivity that is required to solve the energy equation for a turbulent flow situation.*/

int CalcTurbconductivity(Base *base)
{
	int NI,NJ,NK,nZones,numZones,cellindex,flagDiffcoeff;
	int i,j,k,*ptrSize,flagrho,flagspheat,numIface,numZbc,dir;
	int iface,zBc,*ptrRange,index;
	Zone *zone;
	double *rho,*slabrho,*spheat,*slabspheat;
	double *domainPtr,*slabIptr,*slabBcptr;//,*domainTvisc,*slabIvisc,*slabBcvisc;
	//Interface *IfacePtr;
	nZones=((base)->Number_of_zones);
	for(numZones=0;numZones<nZones;numZones++)
	{
		zone = base->CgnsZone+numZones;
		if(zone->flagSolidFluid==0)
		{
			ptrSize = zone->Size;
			NI=*(ptrSize+0);
			NJ=*(ptrSize+1);
			NK=*(ptrSize+2);
			(zone->Flowsol->Tcoeff->Prt)=1.0;
			flagDiffcoeff=zone->Flowsol->flagDiffcoeff;
			
			flagrho = zone->Flowsol->flagRho;
			flagspheat=zone->Flowsol->flagSpecificHeat;
			domainPtr= zone->Flowsol->Tcoeff->TurbConductivity;
			spheat=zone->Flowsol->SpecificHeat;
			rho = zone->Flowsol->Rho;
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				cellindex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
				*(domainPtr+cellindex)= *(spheat+flagspheat*cellindex) * *(zone->Flowsol->Tcoeff->TViscoeff+cellindex)/(zone->Flowsol->Tcoeff->Prt);
				*(zone->Flowsol->Tcoeff->EffTurbConductivity+cellindex)=*(zone->Flowsol->Diffcoeff+flagDiffcoeff*cellindex)+*(domainPtr+cellindex);
			}
						
			numIface=(zone)->Number_of_Interfaces;
			numZbc=(zone)->Number_of_BC_Patches;
			
			for(dir=1;dir<=6;dir++)
			{
				for(iface=0;iface<numIface;iface++)//interfaces
				{
					if(((zone->Iface+iface)->Direction)==dir)
					{
						
						ptrRange=(zone->Iface+iface)->Point_Range;
						slabIptr =(zone->Iface+iface)->ZSlab->TurbConductivity;
						slabspheat=(zone->Iface+iface)->ZSlab->SpecificHeat;
						slabrho = (zone->Iface+iface)->ZSlab->Rho;
						if(dir==1||dir==2)
						{
							for(k=*(ptrRange+2);k < *(ptrRange+5);k++)
								for(j=*(ptrRange+1);j < *(ptrRange+4);j++)
							{
								index=(k-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1))+(j-*(ptrRange+1));
								*(slabIptr+index) =(*(slabspheat+flagspheat*index))*(*((zone->Iface+iface)->ZSlab->TViscoeff+index))/(zone->Flowsol->Tcoeff->Prt); 
								*((zone->Iface+iface)->ZSlab->EffTurbConductivity+index)=*((zone->Iface+iface)->ZSlab->Diffcoeff+flagDiffcoeff*index)+*(slabIptr+index);
								
							}
						}
						if(dir==3||dir==4)
						{
							for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
								for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
							{
								index=(k-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
								*(slabIptr+index) =(*(slabspheat+flagspheat*index))*(*((zone->Iface+iface)->ZSlab->TViscoeff+index))/(zone->Flowsol->Tcoeff->Prt); 
								*((zone->Iface+iface)->ZSlab->EffTurbConductivity+index)=*((zone->Iface+iface)->ZSlab->Diffcoeff+flagDiffcoeff*index)+*(slabIptr+index);
								
							}
	
						}
						if(dir==5||dir==6)
						{
							for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
								for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
							{
								index=(j-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
								*(slabIptr+index) =(*(slabspheat+flagspheat*index))*(*((zone->Iface+iface)->ZSlab->TViscoeff+index))/(zone->Flowsol->Tcoeff->Prt);  
								*((zone->Iface+iface)->ZSlab->EffTurbConductivity+index)=*((zone->Iface+iface)->ZSlab->Diffcoeff+flagDiffcoeff*index)+*(slabIptr+index);
								
							}
						}
						
					}
				}
				for(zBc=0;zBc<numZbc;zBc++)//boundaries
				{
					if((zone->ZoneBc+zBc)->Direction==dir)
					{
						ptrRange=(zone->ZoneBc+zBc)->Point_Range;
						
						slabBcptr =(zone->ZoneBc+zBc)->ZSlab->TurbConductivity;
						slabspheat=(zone->ZoneBc+zBc)->ZSlab->SpecificHeat;
						slabrho = (zone->ZoneBc+zBc)->ZSlab->Rho;
	
						if(dir==1||dir==2)
						{
							for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
								for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
							{
								index=(k-*(ptrRange+2))*(*(ptrRange+4)-*(ptrRange+1))+(j-*(ptrRange+1));
								*(slabBcptr+index)=(*(slabspheat+flagspheat*index))*(*((zone->ZoneBc+zBc)->ZSlab->TViscoeff+index))/(zone->Flowsol->Tcoeff->Prt);
								*((zone->ZoneBc+zBc)->ZSlab->EffTurbConductivity+index)=*((zone->ZoneBc+zBc)->ZSlab->Diffcoeff+flagDiffcoeff*index)+*(slabBcptr+index);
								
							}	
						}
						if(dir==3||dir==4)
						{
							for(k=*(ptrRange+2);k<*(ptrRange+5);k++)
								for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
							{
								index=(k-*(ptrRange+2))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
								*(slabBcptr+index)=(*(slabspheat+flagspheat*index))*(*((zone->ZoneBc+zBc)->ZSlab->TViscoeff+index))/(zone->Flowsol->Tcoeff->Prt);
								*((zone->ZoneBc+zBc)->ZSlab->EffTurbConductivity+index)=*((zone->ZoneBc+zBc)->ZSlab->Diffcoeff+flagDiffcoeff*index)+*(slabBcptr+index);
								
							}
						}
						if(dir==5||dir==6)
						{
							for(j=*(ptrRange+1);j<*(ptrRange+4);j++)
								for(i=*(ptrRange+0);i<*(ptrRange+3);i++)
							{
								index=(j-*(ptrRange+1))*(*(ptrRange+3)-*(ptrRange+0))+(i-*(ptrRange+0));
								*(slabBcptr+index)=(*(slabspheat+flagspheat*index))*(*((zone->ZoneBc+zBc)->ZSlab->TViscoeff+index))/(zone->Flowsol->Tcoeff->Prt);
								*((zone->ZoneBc+zBc)->ZSlab->EffTurbConductivity+index)=*((zone->ZoneBc+zBc)->ZSlab->Diffcoeff+flagDiffcoeff*index)+*(slabBcptr+index);
								
							}
						}
						
					}
				}
			}
		}
	}
	
	//Updating Turbulent Conductivity in the Interfaces of all FLUID Zones
// 	for(numZones=0;numZones<nZones;numZones++)
// 	{
// 		zone = base->CgnsZone+numZones;
	// 
// 		if(zone->flagSolidFluid==0)
// 		{
// 			numIface = zone->Number_of_Interfaces;
	// 
// 			for(iface=0;iface<numIface;iface++)
// 			{
// 				IfacePtr = zone->Iface+iface;
// 				ptrZoneNeighbor = (Zone *)FindZonePtr(IfacePtr->Adjoining_Zone_Name,base);
	// 
// 				if(ptrZoneNeighbor->flagSolidFluid == 0)
// 				{
// 					FillerInterfaceVariables(718,IfacePtr,ptrZoneNeighbor,(IfacePtr->Direction)-1);
// 				}
// 				else
// 				{
// 					//If the adjacent zone is SOLID, boundary condition has to be applied for turbulent viscosity
// 					//i.e., Turbulent Conductivity can be made equal to 'Zero' for all turbulence models except standard and RNG k-epsilon models.
// 					if(base->flagTurbModel==1 || base->flagTurbModel==3)	//Check the Interface Condition for "log version of standard k-epsilon model"
// 						HomogeneousNeumannIC(zone,IfacePtr,718);//NonOrthogonalHomogeneousNeumannIC(base,zone,IfacePtr,718);
// 					else
// 						HomogeneousDirichletIC(718,IfacePtr);
// 				}
// 			}
// 		}
// 	}

	return 0;
}
			





