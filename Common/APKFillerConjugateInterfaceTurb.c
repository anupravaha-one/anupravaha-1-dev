#include "APKFillerConjugateInterfaceTurb.h"

int FillerConjugateInterfaceTurb(Base *base,int variableTag,Interface *facePtr,int iface,Zone *zone)
{
	int i,j,k,targetIndices;
	int *ptrRange,Face,flagRho,flagViscoeff;
	int totalCells,domainCellIndex,NIcell,NJcell,NKcell,NI,NJ,NK;
	double *slabPointer,ConjugateRms=0.0;//,*slabPointerMass,*domainPointer;
	double beta=3.0/40.0,rho,mu,y1;

	iface = facePtr->Direction;
	ptrRange = facePtr->Point_Range;
	
	NIcell = *(ptrRange+3) - *(ptrRange+0);
	NJcell = *(ptrRange+4) - *(ptrRange+1);
	NKcell = *(ptrRange+5) - *(ptrRange+2);
	
	NI = *(zone->Size+0);
	NJ = *(zone->Size+1);
	NK = *(zone->Size+2);
	
	slabPointer = FindSlabPtrIface(facePtr,variableTag);
	//domainPointer=MakeDomainArray(zone,variableTag);
	
	switch(variableTag)
	{
		case 100:
		case 200:
		case 300:
		case 1:
		case 2:
		case 3:
			//slabPointerMass=FindSlabPtrIface(facePtr,variableTag*100);
			//for standard and RNG k-epsilon turbulence models only
			switch(iface)
			{
				case 1://Imax face
					i=NI-2;
					Face=1;//East face
					targetIndices = 0;
					for(k=0; k<NKcell; k++)
					{
						for(j=0; j<NJcell; j++)
						{
							domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
							*(slabPointer+targetIndices) = WallFunctionDirichletConjugate(zone,facePtr,Face,variableTag,targetIndices,domainCellIndex,0.0);
							//*(slabPointerMass+targetIndices) = WallFunctionDirichletConjugate(zone,facePtr,Face,variableTag*100,targetIndices,domainCellIndex,0.0);
							targetIndices++;
						}
					}
					break;
				case 2://Imin face
					i=0;
					Face = 0;//West face
					targetIndices = 0;
					for(k=0; k<NKcell; k++)
					{
						for(j=0; j<NJcell; j++)
						{
							domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
							*(slabPointer+targetIndices) = WallFunctionDirichletConjugate(zone,facePtr,Face,variableTag,targetIndices,domainCellIndex,0.0);
							//*(slabPointerMass+targetIndices) = WallFunctionDirichletConjugate(zone,facePtr,Face,variableTag*100,targetIndices,domainCellIndex,0.0);
							targetIndices++;
						}
					}
					break;
				case 3:
					j = NJ-2;
					Face = 3;
					targetIndices = 0;
					for(k=0; k<NKcell; k++)
					{
						for(i=0; i<NIcell; i++)
						{
							domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
							*(slabPointer+targetIndices) = WallFunctionDirichletConjugate(zone,facePtr,Face,variableTag,targetIndices,domainCellIndex,0.0);
							//*(slabPointerMass+targetIndices) = WallFunctionDirichletConjugate(zone,facePtr,Face,variableTag*100,targetIndices,domainCellIndex,0.0);
							targetIndices++;
						}
					}
					break;
				case 4:
					j = 0;
					Face = 2;
					targetIndices = 0;
					for(k=0; k<NKcell; k++)
					{
						for(i=0; i<NIcell; i++)
						{
							domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
							*(slabPointer+targetIndices) = WallFunctionDirichletConjugate(zone,facePtr,Face,variableTag,targetIndices,domainCellIndex,0.0);
							//*(slabPointerMass+targetIndices) = WallFunctionDirichletConjugate(zone,facePtr,Face,variableTag*100,targetIndices,domainCellIndex,0.0);
							targetIndices++;
						}
					}
					break;
				case 5:
					k = NK-2;
					Face = 5;
					targetIndices = 0;
					for(j=0; j<NJcell; j++)
					{
						for(i=0; i<NIcell; i++)
						{
							domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
							*(slabPointer+targetIndices) = WallFunctionDirichletConjugate(zone,facePtr,Face,variableTag,targetIndices,domainCellIndex,0.0);
							//*(slabPointerMass+targetIndices) = WallFunctionDirichletConjugate(zone,facePtr,Face,variableTag*100,targetIndices,domainCellIndex,0.0);
							targetIndices++;
						}
					}
					break;
				case 6:
					k = 0;
					Face = 4;
					targetIndices = 0;
					for(j=0; j<NJcell; j++)
					{
						for(i=0; i<NIcell; i++)
						{
							domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
							*(slabPointer+targetIndices) = WallFunctionDirichletConjugate(zone,facePtr,Face,variableTag,targetIndices,domainCellIndex,0.0);
							//*(slabPointerMass+targetIndices) = WallFunctionDirichletConjugate(zone,facePtr,Face,variableTag*100,targetIndices,domainCellIndex,0.0);
							targetIndices++;
						}
					}
					break;
			}
			facePtr->ConjugateRms = 0;//<<<< Narsimha
			break;
// 		case 600:
// 			if((base->flagTurbModel ==1)||(base->flagTurbModel ==3))//for standard k-epsilon and RNG models
// 			{
// 				ConjugateRms = NonOrthogonalHomogeneousNeumannIC(base,zone,facePtr,variableTag);
// 			}
// 			else
// 			{
// 				HomogeneousDirichletIC(variableTag,facePtr);
// 			}
// 			break;
		case 601:
			if(base->flagTurbModel == 4)//for Launder-Sharma model
			{
				HomogeneousDirichletIC(variableTag,facePtr);
			}
			else if(base->flagTurbModel == 5)//for Lam-Bremhorst turbulence model
			{
				ConjugateRms = NonOrthogonalHomogeneousNeumannIC(base,zone,facePtr,variableTag);
			}
			else
			{
				switch(iface)
				{
					case 1://Imax face
						totalCells = NJcell * NKcell;
						i=NI-2;
						Face=1;//East face
						targetIndices = 0;
						if(base->flagTurbModel ==1 || base->flagTurbModel ==3)//for standard k-epsilon and RNG models
						{
							for(k=0; k<NKcell; k++)
							{
								for(j=0; j<NJcell; j++)
								{
									domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
									*(slabPointer+targetIndices) = WallFunctionDirichletConjugate(zone,facePtr,Face,variableTag,targetIndices,domainCellIndex,0.0);
									targetIndices++;
								}
							}
						}
						else if(base->flagTurbModel==2 || base->flagTurbModel==6)//for standard k-omega and SST models
						{
							flagViscoeff = zone->Flowsol->flagViscoeff;
							flagRho = zone->Flowsol->flagRho;
							for(targetIndices=0; targetIndices<totalCells; targetIndices++)
							{
								rho = *(facePtr->ZSlab->Rho+targetIndices*flagRho);
								mu  = *(facePtr->ZSlab->Viscoeff+targetIndices*flagViscoeff);
								y1 = *(facePtr->deltaP+targetIndices);
								*(slabPointer+targetIndices) = 60.0*mu/(rho*beta*y1*y1);
							}
						}
						break;
					case 2://Imin face
						totalCells = NJcell * NKcell;
						i=0;
						Face = 0;//West face
						targetIndices = 0;
						if(base->flagTurbModel==1 || base->flagTurbModel==3)//for standard k-epsilon and RNG models
						{
							for(k=0; k<NKcell; k++)
							{
								for(j=0; j<NJcell; j++)
								{
									domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
									*(slabPointer+targetIndices) = WallFunctionDirichletConjugate(zone,facePtr,Face,variableTag,targetIndices,domainCellIndex,0.0);
									targetIndices++;
								}
							}
						}
						else if(base->flagTurbModel==2 || base->flagTurbModel==6)
						{
							flagViscoeff = zone->Flowsol->flagViscoeff;
							flagRho = zone->Flowsol->flagRho;
							for(targetIndices=0; targetIndices<totalCells; targetIndices++)
							{
								rho = *(facePtr->ZSlab->Rho+targetIndices*flagRho);
								mu  = *(facePtr->ZSlab->Viscoeff+targetIndices*flagViscoeff);
								y1 = *(facePtr->deltaP+targetIndices);
								*(slabPointer+targetIndices) = 60.0*mu/(rho*beta*y1*y1);
							}
						}
						break;
					case 3:
						totalCells = NIcell * NKcell;
						j = NJ-2;
						Face = 3;
						targetIndices = 0;
						if((base->flagTurbModel ==1)||(base->flagTurbModel ==3))//for standard k-epsilon and RNG models
						{
							for(k=0; k<NKcell; k++)
							{
								for(i=0; i<NIcell; i++)
								{
									domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
									*(slabPointer+targetIndices) = WallFunctionDirichletConjugate(zone,facePtr,Face,variableTag,targetIndices,domainCellIndex,0.0);
									targetIndices++;
								}
							}
						}
						else if(base->flagTurbModel==2 || base->flagTurbModel==6)
						{
							flagViscoeff = zone->Flowsol->flagViscoeff;
							flagRho = zone->Flowsol->flagRho;
							for(targetIndices=0; targetIndices<totalCells; targetIndices++)
							{
								rho = *(facePtr->ZSlab->Rho+targetIndices*flagRho);
								mu  = *(facePtr->ZSlab->Viscoeff+targetIndices*flagViscoeff);
								y1 = *(facePtr->deltaP+targetIndices);
								*(slabPointer+targetIndices) = 60.0*mu/(rho*beta*y1*y1);
							}
						}
						break;
					case 4:
						totalCells = NIcell * NKcell;
						j = 0;
						Face = 2;
						targetIndices = 0;
						if((base->flagTurbModel ==1)||(base->flagTurbModel ==3))//for standard k-epsilon and RNG models
						{
							for(k=0; k<NKcell; k++)
							{
								for(i=0; i<NIcell; i++)
								{
									domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
									*(slabPointer+targetIndices) = WallFunctionDirichletConjugate(zone,facePtr,Face,variableTag,targetIndices,domainCellIndex,0.0);
									targetIndices++;
								}
							}
						}
						else if(base->flagTurbModel==2 || base->flagTurbModel==6)
						{
							flagViscoeff = zone->Flowsol->flagViscoeff;
							flagRho = zone->Flowsol->flagRho;
							for(targetIndices=0; targetIndices<totalCells; targetIndices++)
							{
								rho = *(facePtr->ZSlab->Rho+targetIndices*flagRho);
								mu  = *(facePtr->ZSlab->Viscoeff+targetIndices*flagViscoeff);
								y1 = *(facePtr->deltaP+targetIndices);
								*(slabPointer+targetIndices) = 60.0*mu/(rho*beta*y1*y1);
							}
						}
						break;
					case 5:
						totalCells = NIcell * NJcell;
						k = NK-2;
						Face = 5;
						targetIndices = 0;
						if((base->flagTurbModel ==1)||(base->flagTurbModel ==3))//for standard k-epsilon and RNG models
						{
							for(j=0; j<NJcell; j++)
							{
								for(i=0; i<NIcell; i++)
								{
									domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
									*(slabPointer+targetIndices) = WallFunctionDirichletConjugate(zone,facePtr,Face,variableTag,targetIndices,domainCellIndex,0.0);
									targetIndices++;
								}
							}
						}
						else if(base->flagTurbModel==2 || base->flagTurbModel==6)
						{
							flagViscoeff = zone->Flowsol->flagViscoeff;
							flagRho = zone->Flowsol->flagRho;
							for(targetIndices=0; targetIndices<totalCells; targetIndices++)
							{
								rho = *(facePtr->ZSlab->Rho+targetIndices*flagRho);
								mu  = *(facePtr->ZSlab->Viscoeff+targetIndices*flagViscoeff);
								y1 = *(facePtr->deltaP+targetIndices);
								*(slabPointer+targetIndices) = 60.0*mu/(rho*beta*y1*y1);
							}
						}
						break;
					case 6:
						totalCells = NIcell * NJcell;
						k = 0;
						Face = 4;
						targetIndices = 0;
						if((base->flagTurbModel ==1)||(base->flagTurbModel ==3))//for standard k-epsilon and RNG models
						{
							for(j=0; j<NJcell; j++)
							{
								for(i=0; i<NIcell; i++)
								{
									domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
									*(slabPointer+targetIndices) = WallFunctionDirichletConjugate(zone,facePtr,Face,variableTag,targetIndices,domainCellIndex,0.0);
									targetIndices++;
								}
							}
						}
						else if(base->flagTurbModel==2 || base->flagTurbModel==6)
						{
							flagViscoeff = zone->Flowsol->flagViscoeff;
							flagRho = zone->Flowsol->flagRho;
							for(targetIndices=0; targetIndices<totalCells; targetIndices++)
							{
								rho = *(facePtr->ZSlab->Rho+targetIndices*flagRho);
								mu  = *(facePtr->ZSlab->Viscoeff+targetIndices*flagViscoeff);
								y1 = *(facePtr->deltaP+targetIndices);
								*(slabPointer+targetIndices) = 60.0*mu/(rho*beta*y1*y1);
							}
						}
						break;
				}
			}
			break;	
	}
	return(0);
}

