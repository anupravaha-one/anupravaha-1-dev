#include "APKExtractVariables.h" 

/**
The function ExtractVariables updates all the interface slab values in the current zone for the input 
variable, the variables here are Velocities, temperature, Scalars and so on. 
It takes variable tag, Base pointer and Zone pointer as input. It loops over all the faces of a Zone and finds 
Interface pointers, and then with the help of function FillerInterfaceVariables() it updates the interface 
slabs. It also checks whether the adjoining zone is also of the same type as the current zone, i.e; whether the adjoining
zone is  a solid or fluid if the current zone is a fluid etc, and calls FillerConjugateInterface() to update and implement
interface conditions if it is a conjugate heat transfer case. When it encounters boundary in the Zone then it calls
the function  update_boundary() or update_boundaryTurb(),according to whether it is a laminar or turbulent flow problem,
which implement the boundary conditions.
 */
int ExtractVariables(int variableTag,Base *base,Zone *zone)
{
	int iface,Tag=1,tag,masstag=0;
	void *ptr;
	double  ConjugateRms=0.0;
	Interface *facePtr;
	ZoneBoundary *bcPtr;
	Zone *ptrZoneNeighbor;
	
	if(variableTag==3100 || variableTag==3200 || variableTag==3300)
		masstag=1;
	if(variableTag==100 || variableTag==200 || variableTag==300)
		masstag=1;
	if(variableTag==6200 || variableTag==6300 || variableTag==6400)
		masstag=1;

	for(iface=0;iface < 6;iface++)
	{
		ptr = zone->ZONE_FACES[iface];
		tag = zone->BC_IFACE_TAG[iface];
		
		if(tag==2)
		{
			facePtr=(Interface *)ptr;
			
			ptrZoneNeighbor= base->CgnsZone+(zone->ZONE_NEIGHBOUR_NUM[iface]);
			
			Tag=strcmp(((zone)->SolidFluidName),(ptrZoneNeighbor->SolidFluidName));
			if(Tag!=0) //Conjugate Heat Transfer Case
			{
				switch(variableTag)
				{
					case 100:
					case 200:
					case 300:
					case 1:
					case 2:
					case 3:
					case 31:	case 32:	case 33:
					case 41:	case 42:	case 43:
					case 44:	case 45:	case 46:
					case 47:	case 48:	case 49:
					case 50:	case 51:	case 52:
					case 53:	case 54:	case 55:
						if(((zone->flagSolidFluid)-(ptrZoneNeighbor->flagSolidFluid))!=0) //May be '<' symbol should be used in place of !=0 Check!!! <<<< Narsimha.
						{
							if(zone->flagSolidFluid == 0 && (base->flagTurbModel ==1 || base->flagTurbModel ==3))
							{
								
								FillerConjugateInterfaceTurb(base,variableTag,facePtr,iface+1,zone);
							}
							else
								FillerConjugateInterface(variableTag,facePtr,iface);
						}
						else 
						{
							//for Slip BC Implimentation
							FillerIfaceConjugateVariable(base,variableTag,zone,ptrZoneNeighbor,facePtr,iface);
						}
						break;
					case 4://updating Pressure at interfaces
					case 34:	case 35:	case 36:	case 37:	case 38:
						if(zone->flagSolidFluid == 0)
						{
							if(ptrZoneNeighbor->flagSolidFluid == 1)
								ConjugateRms = NonOrthogonalHomogeneousNeumannIC(base,zone,facePtr,variableTag);//NonOrthogonalNeumannIC(base,zone,facePtr,variableTag,0.0);
							else
								FillerInterfaceVariables(variableTag,facePtr,ptrZoneNeighbor,iface);
						}
						else 
						{
							HomogeneousDirichletIC(variableTag,facePtr);
						}
						break;
					case 990:
						ConjugateRms = NonOrthogonalHomogeneousNeumannIC(base,zone,facePtr,variableTag);
						break;
					case 6:
						if((zone->flagSolidFluid - ptrZoneNeighbor->flagSolidFluid != 0) && (base->flagTurbModel ==1 || base->flagTurbModel ==3))
						{
							FillerIfaceConjugateTemperatureTurb(base,variableTag,zone,ptrZoneNeighbor,facePtr);
// 							printf("ConjugateRms for Temperature = %e\n", facePtr->ConjugateRms);
						}
						else
						{
							FillerIfaceConjugateVariable(base,variableTag,zone,ptrZoneNeighbor,facePtr,iface);
						}
						break;
					case 600:
// 						if(base->flagTurb == 1 && zone->flagSolidFluid == 0 && ptrZoneNeighbor->flagSolidFluid ==1)
						if(base->BaseParameters->Turb->numTurbScalars >= 1 && facePtr->ConjugateIfaceTag == 1)//for turbulence case with interface corresponding to FLUID zone
						{
							if(base->flagTurbModel ==1 || base->flagTurbModel ==3)//for standard or RNG k-epsilon models
							{
								//HomogeneousNeumannIC(zone,facePtr,variableTag);
								ConjugateRms = NonOrthogonalHomogeneousNeumannIC(base,zone,facePtr,variableTag);
							}
							else//for all other turbulence models
							{
								HomogeneousDirichletIC(variableTag,facePtr);
							}
						}
						else
						{
							FillerInterfaceVariables(variableTag,facePtr,ptrZoneNeighbor,iface);
						}
						break;
					case 601:
// 						if(base->flagTurb == 1 && facePtr->ConjugateIfaceTag == 1)
						if(base->BaseParameters->Turb->numTurbScalars >= 2 && facePtr->ConjugateIfaceTag == 1)
						{
							FillerConjugateInterfaceTurb(base,variableTag,facePtr,iface+1,zone);
						}
						else
						{
							FillerInterfaceVariables(variableTag,facePtr,ptrZoneNeighbor,iface);
						}
						break;
					case 17://MHD Potential, naveen
					case 56://Magnetic Field X, rajendra
					case 57://Magnetic Field X, rajendra
					case 58://Magnetic Field X, rajendra
				                                  
						//Check whether variableTag 17 is included in find_convcri() function
						FillerIfaceConjugateVariable(base,variableTag,zone,ptrZoneNeighbor,facePtr,iface);
						break;
			               	case 59: //Psi
                                                FillerInterfaceVariables(variableTag,facePtr,ptrZoneNeighbor,iface);
						break;

					default:
						FillerInterfaceVariables(variableTag,facePtr,ptrZoneNeighbor,iface);
						break;
				}
			}
			else
			{
				if(variableTag!=16)
					FillerInterfaceVariables(variableTag,facePtr,ptrZoneNeighbor,iface);
				else
					FillerInterfaceRadiationVariables(base,variableTag,facePtr,ptrZoneNeighbor,iface+1);
			}
		}
		else if(tag==1)
		{
			bcPtr=(ZoneBoundary *)ptr;
			if(masstag==0)
				update_boundary(base,variableTag,zone,bcPtr);
			else
				velocitymassBc(base,variableTag/100,zone,bcPtr);
			//In Conjugate Heat Transfer case, if the zone is SOLID update_boundary is called for the boundaries??
			
		}
		
	}
	
	
	return 0;
}



//--------------------------- Buffer Function used for temporary purpose <<<<< by Narsimha -------------
// #include "APKExtractVariables.h"
// 
// /**
// The function ExtractVariables updates all the interface slab values in the current zone for the input 
// variable, the variables here are Velocities, temperature, Scalars and so on. 
// It takes variable tag, Base pointer and Zone pointer as input. It loops over all the faces of a Zone and finds 
// Interface pointers, and then with the help of function FillerInterfaceVariables() it updates the interface 
// slabs. It also checks whether the adjoining zone is also of the same type as the current zone, i.e; whether the adjoining
// zone is  a solid or fluid if the current zone is a fluid etc, and calls FillerConjugateInterface() to update and implement
// interface conditions if it is a conjugate heat transfer case. When it encounters boundary in the Zone then it calls
// the function  update_boundary() or update_boundaryTurb(),according to whether it is a laminar or turbulent flow problem,
// which implement the boundary conditions.
//  */
// int ExtractVariables(int variableTag,Base *base,Zone *zone)
// {
// 	int iface,Tag=1;
// 	char tag,name[33];
// 	void *ptr;
// 	double  ConjugateRms=0.0;
// 	Interface *facePtr;
// 	ZoneBoundary *bcPtr;
// 	Zone *ptrZoneNeighbor;
// 
// 	for(iface=0;iface < 6;iface++)
// 	{
// 		ptr = (void *)FindInterfaceBCPtrMod(zone,iface+1,name,&tag);
// 		if(tag=='i')
// 		{
// 			facePtr=(Interface *)ptr;
// 			
// 			ptrZoneNeighbor = (Zone *)FindZonePtr(name,base);/* working properly*/
// 			
// 			Tag=strcmp(((zone)->SolidFluidName),(ptrZoneNeighbor->SolidFluidName));
// 			if(Tag!=0) //Two zones are of different material (may not be Conjugate Heat Transfer Case)
// 			{
// 				switch(variableTag)
// 				{
// //					case 100:
// //					case 200:
// //					case 300:
// 					case 1:
// 					case 2:
// 					case 3:
// 						//Conjugate Heat Transfer Case
// 						if(((zone->flagSolidFluid)-(ptrZoneNeighbor->flagSolidFluid))!=0) //May be '<' symbol should be used in place of !=0 Check!!! <<<< Narsimha.
// 						{
// 							if(zone->flagSolidFluid == 0 && (base->flagTurbModel ==1 || base->flagTurbModel ==3))
// 								FillerConjugateInterfaceTurb(base,variableTag,facePtr,iface+1,zone);
// 							else
// 								FillerConjugateInterface(variableTag,facePtr,iface);
// 						}
// 						else //Two zones are of different fluids or solids
// 						{
// 							//for Slip BC Implementation
// 							FillerIfaceConjugateVariable(base,variableTag,zone,ptrZoneNeighbor,facePtr,iface);
// 						}
// 						break;
// 					case 4://updating Pressure at interfaces
// 						if(zone->flagSolidFluid == 0)
// 						{
// 							if(ptrZoneNeighbor->flagSolidFluid == 1)
// 								ConjugateRms = NonOrthogonalHomogeneousNeumannIC(base,zone,facePtr,variableTag);//NonOrthogonalNeumannIC(base,zone,facePtr,variableTag,0.0);
// 							else
// 								FillerInterfaceVariables(variableTag,facePtr,ptrZoneNeighbor,iface);
// 						}
// 						else 
// 						{
// 							HomogeneousDirichletIC(variableTag,facePtr);
// 						}
// 						break;
// 					case 990:
// 						ConjugateRms = NonOrthogonalHomogeneousNeumannIC(base,zone,facePtr,variableTag);
// 						break;
// 					case 6:
// 						if((zone->flagSolidFluid - ptrZoneNeighbor->flagSolidFluid != 0) && (base->flagTurbModel ==1 || base->flagTurbModel ==3))
// 						{
// 							FillerIfaceConjugateTemperatureTurb(base,variableTag,zone,ptrZoneNeighbor,facePtr);
// // 							printf("ConjugateRms for Temperature = %e\n", facePtr->ConjugateRms);
// 						}
// 						else
// 						{
// 							FillerIfaceConjugateVariable(base,variableTag,zone,ptrZoneNeighbor,facePtr,iface);
// 						}
// 						break;
// 					case 600:
// // 						if(base->flagTurb == 1 && zone->flagSolidFluid == 0 && ptrZoneNeighbor->flagSolidFluid ==1)
// 						if(base->flagTurb == 1 && facePtr->ConjugateIfaceTag == 1)//for turbulence case with interface corresponding to FLUID zone
// 						{
// 							if(base->flagTurbModel ==1 || base->flagTurbModel ==3)//for standard or RNG k-epsilon models
// 							{
// 								//HomogeneousNeumannIC(zone,facePtr,variableTag);
// 								ConjugateRms = NonOrthogonalHomogeneousNeumannIC(base,zone,facePtr,variableTag);
// 							}
// 							else//for all other turbulence models
// 							{
// 								HomogeneousDirichletIC(variableTag,facePtr);
// 							}
// 						}
// 						else
// 						{
// 							FillerInterfaceVariables(variableTag,facePtr,ptrZoneNeighbor,iface);
// 						}
// 						break;
// 					case 601:
// // 						if(base->flagTurb == 1 && zone->flagSolidFluid == 0 && ptrZoneNeighbor->flagSolidFluid ==1)
// 						if(base->flagTurb == 1 && facePtr->ConjugateIfaceTag == 1)
// 						{
// 							FillerConjugateInterfaceTurb(base,variableTag,facePtr,iface+1,zone);
// 						}
// 						else
// 						{
// 							FillerInterfaceVariables(variableTag,facePtr,ptrZoneNeighbor,iface);
// 						}
// 						break;
// 					case 17://MHD Potential, naveen
// 					case 56://Magnetic Field X, rajendra
// 					case 57://Magnetic Field X, rajendra
// 					case 58://Magnetic Field X, rajendra
// 					case 59://Psi, rajendra
// 						//Check whether variableTag 17 is included in find_convcri() function
// 						FillerIfaceConjugateVariable(base,variableTag,zone,ptrZoneNeighbor,facePtr,iface);
// 						break;
// 					default:
// 						FillerInterfaceVariables(variableTag,facePtr,ptrZoneNeighbor,iface);
// 						break;
// 				}
// 			}
// 			else
// 			{
// 				if(variableTag!=16)
// 					FillerInterfaceVariables(variableTag,facePtr,ptrZoneNeighbor,iface);
// 				else
// 					FillerInterfaceRadiationVariables(base,variableTag,facePtr,ptrZoneNeighbor,iface+1);
// 			}
// 		}
// 		else
// 		if(tag=='b')
// 		{
// 			bcPtr=(ZoneBoundary *)ptr;
// 			update_boundary(base,variableTag,zone,bcPtr);
// 			//In Conjugate Heat Transfer case, if the zone is SOLID update_boundary is called for the boundaries??
// 		}
// 	}
// 	return 0;
// }
