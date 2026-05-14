#include "APKTURBCalcUtau.h"
/**
This function CalcUtau calculates the value of uTau(Friction Velocity). The function loops over all the boundaries and calculates uTau
only for those boundaries having a boundary Tag of 22 (WallViscous).
 */
int CalcUtau(Base *base)
{
	int i,j,k,domainCellIndex,NI,NJ,NK,slabCellIndex,Face;
	int numZones,countZone,numIface,countIface,numBc,countBc,flagrho,flagviscoeff;
	double EffViscosity,rho,u,y,Cmu,K;
	Zone *z;
	ZoneBoundary *zBc;
	Interface *IfacePtr;

	numZones=(base->Number_of_zones);

	if(base->flagTurbModel == 1 || base->flagTurbModel == 3)
	{
		if(base->flagLogKEps == 0)
		{
			for(countZone=0;countZone<numZones;countZone++)
			{
				z = (base->CgnsZone+countZone);//printf("\n zone name-%s",z->zonename);
				if(z->flagSolidFluid == 0)
				{
					NI = *(z->Size+0);
					NJ = *(z->Size+1);
					NK = *(z->Size+2);
	
					numBc = z->Number_of_BC_Patches;
					numIface = z->Number_of_Interfaces;
	
					Cmu = z->Flowsol->Tcoeff->Cmu;
					for(countBc=0;countBc<numBc;countBc++)
					{
						zBc = z->ZoneBc+countBc;
						if(zBc->Boundary_Tag == 22)
						{
							slabCellIndex = 0;
							switch(zBc->Direction)
							{
								case 1: //IMax face,East face
									i=NI-2;
									for (k=0;k<NK-1;k++) 
									for (j=0;j<NJ-1;j++)
									{
										domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
										K = *(z->Flowsol->Scal[0]+domainCellIndex);
										*(zBc->uTau+slabCellIndex) = sqrt(sqrt(Cmu) * K);
										//printf("\n BC Name-%s$Utau-%lf",zBc->BCName,*(zBc->uTau+slabCellIndex));
										slabCellIndex++;
									}
									break;
								case 2: //Imin face, West face
									i=0;
									for (k=0;k<NK-1;k++) 
									for (j=0;j<NJ-1;j++)
									{
										domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
										K = *(z->Flowsol->Scal[0]+domainCellIndex);
										*(zBc->uTau+slabCellIndex) = sqrt(sqrt(Cmu) * K);
										//printf("\n BC Name-%s$Utau-%lf",zBc->BCName,*(zBc->uTau+slabCellIndex));
										slabCellIndex++;
									}
									break;
								case 3: //JMax face, North face
									j=NJ-2;
									for (k=0;k<NK-1;k++)
									for (i=0;i<NI-1;i++)
									{
										domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
										K = *(z->Flowsol->Scal[0]+domainCellIndex);
										*(zBc->uTau+slabCellIndex) = sqrt(sqrt(Cmu) * K);
										//printf("\n BC Name-%s$Utau-%lf",zBc->BCName,*(zBc->uTau+slabCellIndex));
										slabCellIndex++;
									}
									break;
								case 4: //Jmin face,South face
									j=0;
									for (k=0;k<NK-1;k++)
									for (i=0;i<NI-1;i++)
									{
										domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
										K = *(z->Flowsol->Scal[0]+domainCellIndex);
										*(zBc->uTau+slabCellIndex) = sqrt(sqrt(Cmu) * K);
										//printf("\n BC Name-%s$Utau-%lf",zBc->BCName,*(zBc->uTau+slabCellIndex));
										slabCellIndex++;
									}
									break;
								case 5: //KMax face, Top face
									k=NK-2;
									for (j=0;j<NJ-1;j++)
									for (i=0;i<NI-1;i++)
									{
										domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
										K = *(z->Flowsol->Scal[0]+domainCellIndex);
										*(zBc->uTau+slabCellIndex) = sqrt(sqrt(Cmu) * K);
										//printf("\n BC Name-%s$Utau-%lf",zBc->BCName,*(zBc->uTau+slabCellIndex));
										slabCellIndex++;
									}
									break;
								case 6: //Kmin face, Bottom face
									k=0;
									for (j=0;j<NJ-1;j++)
									for (i=0;i<NI-1;i++)
									{
										domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
										K = *(z->Flowsol->Scal[0]+domainCellIndex);
										*(zBc->uTau+slabCellIndex) = sqrt(sqrt(Cmu) * K);
										//printf("\n BC Name-%s$Utau-%lf",zBc->BCName,*(zBc->uTau+slabCellIndex));
										slabCellIndex++;
									}
									break;
							}
						}//End of if loop
					}//End of looping over Boundary Patches
					
					for(countIface=0; countIface<numIface; countIface++)
					{
						IfacePtr = z->Iface+countIface;
						if(IfacePtr->ConjugateIfaceTag == 1)
						{
							slabCellIndex = 0;
							switch(IfacePtr->Direction)
							{
								case 1: //IMax face,East face
									i=NI-2;
									for (k=0;k<NK-1;k++) 
									for (j=0;j<NJ-1;j++)
									{
										domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
										K = *(z->Flowsol->Scal[0]+domainCellIndex);
										*(IfacePtr->uTau+slabCellIndex) = sqrt(sqrt(Cmu) * K);
										slabCellIndex++;
									}
									break;
								case 2: //Imin face, West face
									i=0;
									for (k=0;k<NK-1;k++) 
									for (j=0;j<NJ-1;j++)
									{
										domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
										K = *(z->Flowsol->Scal[0]+domainCellIndex);
										*(IfacePtr->uTau+slabCellIndex) = sqrt(sqrt(Cmu) * K);
										slabCellIndex++;
									}
									break;
								case 3: //JMax face, North face
									j=NJ-2;
									for (k=0;k<NK-1;k++)
									for (i=0;i<NI-1;i++)
									{
										domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
										K = *(z->Flowsol->Scal[0]+domainCellIndex);
										*(IfacePtr->uTau+slabCellIndex) = sqrt(sqrt(Cmu) * K);
										slabCellIndex++;
									}
									break;
								case 4: //Jmin face,South face
									j=0;
									for (k=0;k<NK-1;k++)
									for (i=0;i<NI-1;i++)
									{
										domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
										K = *(z->Flowsol->Scal[0]+domainCellIndex);
										*(IfacePtr->uTau+slabCellIndex) = sqrt(sqrt(Cmu) * K);
										slabCellIndex++;
									}
									break;
								case 5: //KMax face, Top face
									k=NK-2;
									for (j=0;j<NJ-1;j++)
									for (i=0;i<NI-1;i++)
									{
										domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
										K = *(z->Flowsol->Scal[0]+domainCellIndex);
										*(IfacePtr->uTau+slabCellIndex) = sqrt(sqrt(Cmu) * K);
										slabCellIndex++;
									}
									break;
								case 6: //Kmin face, Bottom face
									k=0;
									for (j=0;j<NJ-1;j++)
									for (i=0;i<NI-1;i++)
									{
										domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
										K = *(z->Flowsol->Scal[0]+domainCellIndex);
										*(IfacePtr->uTau+slabCellIndex) = sqrt(sqrt(Cmu) * K);
										slabCellIndex++;
									}
									break;
							}
						}
					}
				}//End of flagSolidFluid check
			}//End of Zone looping
		}
		else if(base->flagLogKEps == 1)
		{	//Calculating Utau using Logarithmic value of Turbulent Kinetic Energy
			for(countZone=0;countZone<numZones;countZone++)
			{
				z = (base->CgnsZone+countZone);
				if(z->flagSolidFluid == 0)
				{
					NI = *(z->Size+0);
					NJ = *(z->Size+1);
					NK = *(z->Size+2);
	
					numBc = z->Number_of_BC_Patches;
					numIface = z->Number_of_Interfaces;
	
					for(countBc=0;countBc<numBc;countBc++)
					{
						zBc = z->ZoneBc+countBc;
						if(zBc->Boundary_Tag == 22)
						{
							switch(zBc->Direction)
							{
								case 1: //IMax face,East face
									for (k=0;k<NK-1;k++) 
										for (j=0;j<NJ-1;j++)
									{
										i=NI-2;
										domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
										slabCellIndex = (NJ-1)*k+j;
										*(zBc->uTau+slabCellIndex) = sqrt(sqrt(z->Flowsol->Tcoeff->Cmu) * exp(*(z->Flowsol->Scal[0]+domainCellIndex)));
									}
									break;
									case 2: //Imin face, West face
										for (k=0;k<NK-1;k++) 
										for (j=0;j<NJ-1;j++)
										{
										i=0;
										domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
										slabCellIndex = (NJ-1)*k+j;
										*(zBc->uTau+slabCellIndex) = sqrt(sqrt(z->Flowsol->Tcoeff->Cmu) * exp(*(z->Flowsol->Scal[0]+domainCellIndex)));
										}
										break;
										case 3: //JMax face, North face
										for (k=0;k<NK-1;k++)
										for (i=0;i<NI-1;i++)
										{
										j=NJ-2;
										domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
										slabCellIndex = (NI-1)*k+i;
										*(zBc->uTau+slabCellIndex) = sqrt(sqrt(z->Flowsol->Tcoeff->Cmu) * exp(*(z->Flowsol->Scal[0]+domainCellIndex)));
										}
										break;
										case 4: //Jmin face,South face
										for (k=0;k<NK-1;k++)
										for (i=0;i<NI-1;i++)
										{
										j=0;
										domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
										slabCellIndex = (NI-1)*k+i;
										*(zBc->uTau+slabCellIndex) = sqrt(sqrt(z->Flowsol->Tcoeff->Cmu) * exp(*(z->Flowsol->Scal[0]+domainCellIndex)));
										}
										break;
										case 5: //KMax face, Top face
										for (j=0;j<NJ-1;j++)
										for (i=0;i<NI-1;i++)
										{
										k=NK-2;
										domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
										slabCellIndex = (NI-1)*j+i;
										*(zBc->uTau+slabCellIndex) = sqrt(sqrt(z->Flowsol->Tcoeff->Cmu) * exp(*(z->Flowsol->Scal[0]+domainCellIndex)));
										}
										break;
										case 6: //Kmin face, Bottom face
										for (j=0;j<NJ-1;j++)
										for (i=0;i<NI-1;i++)
										{
										k=0;
										domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
										slabCellIndex = (NI-1)*j+i;
										*(zBc->uTau+slabCellIndex) = sqrt(sqrt(z->Flowsol->Tcoeff->Cmu) * exp(*(z->Flowsol->Scal[0]+domainCellIndex)));
										}
										break;
							}
						}//End of if loop
					}//End of looping over Boundary Patches
					for(countIface=0; countIface<numIface; countIface++)
					{
						IfacePtr = z->Iface+countIface;
						if(IfacePtr->ConjugateIfaceTag == 1)
						{
							switch(IfacePtr->Direction)
							{
								case 1: //IMax face,East face
									i=NI-2;
									for (k=0;k<NK-1;k++) 
										for (j=0;j<NJ-1;j++)
									{
										domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
										slabCellIndex = (NJ-1)*k+j;
										*(IfacePtr->uTau+slabCellIndex) = sqrt(sqrt(z->Flowsol->Tcoeff->Cmu) * exp(*(z->Flowsol->Scal[0]+domainCellIndex)));
									}
									break;
									case 2: //Imin face, West face
										i=0;
										for (k=0;k<NK-1;k++) 
										for (j=0;j<NJ-1;j++)
										{
										domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
										slabCellIndex = (NJ-1)*k+j;
										*(IfacePtr->uTau+slabCellIndex) = sqrt(sqrt(z->Flowsol->Tcoeff->Cmu) * exp(*(z->Flowsol->Scal[0]+domainCellIndex)));
										}
										break;
										case 3: //JMax face, North face
										j=NJ-2;
										for (k=0;k<NK-1;k++)
										for (i=0;i<NI-1;i++)
										{
										domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
										slabCellIndex = (NI-1)*k+i;
										*(IfacePtr->uTau+slabCellIndex) = sqrt(sqrt(z->Flowsol->Tcoeff->Cmu) * exp(*(z->Flowsol->Scal[0]+domainCellIndex)));
										}
										break;
										case 4: //Jmin face,South face
										j=0;
										for (k=0;k<NK-1;k++)
										for (i=0;i<NI-1;i++)
										{
										domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
										slabCellIndex = (NI-1)*k+i;
										*(IfacePtr->uTau+slabCellIndex) = sqrt(sqrt(z->Flowsol->Tcoeff->Cmu) * exp(*(z->Flowsol->Scal[0]+domainCellIndex)));
										}
										break;
										case 5: //KMax face, Top face
										k=NK-2;
										for (j=0;j<NJ-1;j++)
										for (i=0;i<NI-1;i++)
										{
										domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
										slabCellIndex = (NI-1)*j+i;
										*(IfacePtr->uTau+slabCellIndex) = sqrt(sqrt(z->Flowsol->Tcoeff->Cmu) * exp(*(z->Flowsol->Scal[0]+domainCellIndex)));
										}
										break;
										case 6: //Kmin face, Bottom face
										k=0;
										for (j=0;j<NJ-1;j++)
										for (i=0;i<NI-1;i++)
										{
										domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
										slabCellIndex = (NI-1)*j+i;
										*(IfacePtr->uTau+slabCellIndex) = sqrt(sqrt(z->Flowsol->Tcoeff->Cmu) * exp(*(z->Flowsol->Scal[0]+domainCellIndex)));
										}
										break;
							}
						}
					}
				}//End of flagSolidFluid check
			}//End of Zone looping
		}
	}
	else //if(base->flagTurbModel == 2 || base->flagTurbModel == 6 || base->flagTurbModel == 9)
	{
		for(countZone=0;countZone<numZones;countZone++)
		{
			z = (base->CgnsZone+countZone);
			if(z->flagSolidFluid == 0)
			{
				NI = *(z->Size+0);
				NJ = *(z->Size+1);
				NK = *(z->Size+2);

				numBc = z->Number_of_BC_Patches;
				numIface = z->Number_of_Interfaces;
				
				flagrho = z->Flowsol->flagRho;
				flagviscoeff = z->Flowsol->flagViscoeff;
				for(countBc=0;countBc<numBc;countBc++)
				{
					zBc = z->ZoneBc+countBc;
					if(zBc->Boundary_Tag == 22)
					{
						switch(zBc->Direction)
						{
							case 1: //IMax face,East face
								i=NI-2;
								for (k=0;k<NK-1;k++) 
								for (j=0;j<NJ-1;j++)
								{
									domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
									slabCellIndex = (NJ-1)*k+j;
									
									//EffViscosity = *(z->Flowsol->Tcoeff->EffViscoeff+domainCellIndex);
									EffViscosity = *(z->Flowsol->Viscoeff+flagviscoeff*domainCellIndex);
									rho = *(z->PrevFlowsol->Rho+domainCellIndex*flagrho);
									y = *(zBc->deltaP+slabCellIndex);
									//u = *(z->Flowsol->v+domainCellIndex);		//This is NOT a generalised formulation
									
									Face = 1;
									u = CalcmoduPt(z,domainCellIndex,Face);	//This is a generalised formulation

									//This is true for a non moving solid wall boundary
									//if wall has some velocity the u/y must be replaced with (u-uwall)/y;
									*(zBc->uTau+slabCellIndex) = sqrt(EffViscosity/rho * u/y);
								}
								break;
							case 2: //Imin face, West face
								i=0;
								for (k=0;k<NK-1;k++) 
								for (j=0;j<NJ-1;j++)
								{
									domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
									slabCellIndex = (NJ-1)*k+j;
								
									//EffViscosity = *(z->Flowsol->Tcoeff->EffViscoeff+domainCellIndex);
									EffViscosity = *(z->Flowsol->Viscoeff+flagviscoeff*domainCellIndex);
									rho = *(z->PrevFlowsol->Rho+domainCellIndex*flagrho);
									y = *(zBc->deltaP+slabCellIndex);
									//u = *(z->Flowsol->v+domainCellIndex);
								
									Face = 0;
									u = CalcmoduPt(z,domainCellIndex,Face);
								
									//This is true for a non moving solid wall boundary
									//if wall has some velocity the u/y must be replace with (u-uwall)/y;
									*(zBc->uTau+slabCellIndex) = sqrt(EffViscosity/rho * u/y);
								}
								break;
							case 3: //JMax face, North face
								j=NJ-2;
								for (k=0;k<NK-1;k++)
								for (i=0;i<NI-1;i++)
								{
									domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
									slabCellIndex = (NI-1)*k+i;
								
									//EffViscosity = *(z->Flowsol->Tcoeff->EffViscoeff+domainCellIndex);
									EffViscosity = *(z->Flowsol->Viscoeff+flagviscoeff*domainCellIndex);
									rho = *(z->PrevFlowsol->Rho+domainCellIndex*flagrho);
									y = *(zBc->deltaP+slabCellIndex);
									u = *(z->Flowsol->u+domainCellIndex);
								
									Face = 3;
									u = CalcmoduPt(z,domainCellIndex,Face);
									
									//This is true for a non moving solid wall boundary
									//if wall has some velocity the u/y must be replace with (u-uwall)/y;
									*(zBc->uTau+slabCellIndex) = sqrt(EffViscosity/rho * u/y);
								}
								break;
							case 4: //Jmin face,South face
								j=0;
								for (k=0;k<NK-1;k++)
								for (i=0;i<NI-1;i++)
								{
									domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
									slabCellIndex = (NI-1)*k+i;
								
									//EffViscosity = *(z->Flowsol->Tcoeff->EffViscoeff+domainCellIndex);
									EffViscosity = *(z->Flowsol->Viscoeff+flagviscoeff*domainCellIndex);
									rho = *(z->PrevFlowsol->Rho+domainCellIndex*flagrho);
									y = *(zBc->deltaP+slabCellIndex);
									u = *(z->Flowsol->u+domainCellIndex);
								
									Face = 2;
									u = CalcmoduPt(z,domainCellIndex,Face);
								
									//This is true for a non moving solid wall boundary
									//if wall has some velocity the u/y must be replace with (u-uwall)/y;
									*(zBc->uTau+slabCellIndex) = sqrt(EffViscosity/rho * u/y);
								}
								break;
							case 5: //KMax face, Top face
								k=NK-2;
								for (j=0;j<NJ-1;j++)
								for (i=0;i<NI-1;i++)
								{
									domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
									slabCellIndex = (NI-1)*j+i;
								
									//EffViscosity = *(z->Flowsol->Tcoeff->EffViscoeff+domainCellIndex);
									EffViscosity = *(z->Flowsol->Viscoeff+flagviscoeff*domainCellIndex);
									rho = *(z->PrevFlowsol->Rho+domainCellIndex*flagrho);
									y = *(zBc->deltaP+slabCellIndex);
									//u = *(z->Flowsol->u+domainCellIndex);
								
									Face = 5;
									u = CalcmoduPt(z,domainCellIndex,Face);
								
									//This is true for a non moving solid wall boundary
									//if wall has some velocity the u/y must be replace with (u-uwall)/y;
									*(zBc->uTau+slabCellIndex) = sqrt(EffViscosity/rho * u/y);
								}
								break;
							case 6: //Kmin face, Bottom face
								k=0;
								for (j=0;j<NJ-1;j++)
								for (i=0;i<NI-1;i++)
								{
									domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
									slabCellIndex = (NI-1)*j+i;
									
									//EffViscosity = *(z->Flowsol->Tcoeff->EffViscoeff+domainCellIndex);
									EffViscosity = *(z->Flowsol->Viscoeff+flagviscoeff*domainCellIndex);
									rho = *(z->PrevFlowsol->Rho+domainCellIndex*flagrho);
									y = *(zBc->deltaP+slabCellIndex);
									//u = *(z->Flowsol->u+domainCellIndex);
								
									Face = 4;
									u = CalcmoduPt(z,domainCellIndex,Face);
									
									//This is true for a non moving solid wall boundary
									//if wall has some velocity the u/y must be replace with (u-uwall)/y;
									*(zBc->uTau+slabCellIndex) = sqrt(EffViscosity/rho * u/y);
								}
								break;
						}
					}//End of if loop
				}//End of looping over Boundary Patches
				
				for(countIface=0; countIface<numIface; countIface++)
				{
					IfacePtr = z->Iface+countIface;
					if(IfacePtr->ConjugateIfaceTag == 1)
					{
						switch(IfacePtr->Direction)
						{
							case 1: //IMax face,East face
								i=NI-2;
								for (k=0;k<NK-1;k++) 
								for (j=0;j<NJ-1;j++)
								{
									domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
									slabCellIndex = (NJ-1)*k+j;
									
									//EffViscosity = *(z->Flowsol->Tcoeff->EffViscoeff+domainCellIndex);
									EffViscosity = *(z->Flowsol->Viscoeff+flagviscoeff*domainCellIndex);
									rho = *(z->PrevFlowsol->Rho+domainCellIndex*flagrho);
									y = *(IfacePtr->deltaP+slabCellIndex);
									//u = *(z->Flowsol->v+domainCellIndex);		//This is NOT a generalised formulation
									
									Face = 1;
									u = CalcmoduPt(z,domainCellIndex,Face);	//This is a generalised formulation

									//This is true for a non moving solid wall boundary
									//if wall has some velocity the u/y must be replaced with (u-uwall)/y;
									*(IfacePtr->uTau+slabCellIndex) = sqrt(EffViscosity/rho * u/y);
								}
								break;
							case 2: //Imin face, West face
								i=0;
								for (k=0;k<NK-1;k++) 
								for (j=0;j<NJ-1;j++)
								{
									domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
									slabCellIndex = (NJ-1)*k+j;
									
									//EffViscosity = *(z->Flowsol->Tcoeff->EffViscoeff+domainCellIndex);
									EffViscosity = *(z->Flowsol->Viscoeff+flagviscoeff*domainCellIndex);
									rho = *(z->PrevFlowsol->Rho+domainCellIndex*flagrho);
									y = *(IfacePtr->deltaP+slabCellIndex);
									//u = *(z->Flowsol->v+domainCellIndex);
									
									Face = 0;
									u = CalcmoduPt(z,domainCellIndex,Face);
									
									//This is true for a non moving solid wall boundary
									//if wall has some velocity the u/y must be replace with (u-uwall)/y;
									*(IfacePtr->uTau+slabCellIndex) = sqrt(EffViscosity/rho * u/y);
								}
								break;
							case 3: //JMax face, North face
								j=NJ-2;
								for (k=0;k<NK-1;k++)
								for (i=0;i<NI-1;i++)
								{
									domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
									slabCellIndex = (NI-1)*k+i;
								
									//EffViscosity = *(z->Flowsol->Tcoeff->EffViscoeff+domainCellIndex);
									EffViscosity = *(z->Flowsol->Viscoeff+flagviscoeff*domainCellIndex);
									rho = *(z->PrevFlowsol->Rho+domainCellIndex*flagrho);
									y = *(IfacePtr->deltaP+slabCellIndex);
									u = *(z->Flowsol->u+domainCellIndex);
							
									Face = 3;
									u = CalcmoduPt(z,domainCellIndex,Face);
								
									//This is true for a non moving solid wall boundary
									//if wall has some velocity the u/y must be replace with (u-uwall)/y;
									*(IfacePtr->uTau+slabCellIndex) = sqrt(EffViscosity/rho * u/y);
								}
								break;
							case 4: //Jmin face,South face
								j=0;
								for (k=0;k<NK-1;k++)
								for (i=0;i<NI-1;i++)
								{
									domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
									slabCellIndex = (NI-1)*k+i;
									
									//EffViscosity = *(z->Flowsol->Tcoeff->EffViscoeff+domainCellIndex);
									EffViscosity = *(z->Flowsol->Viscoeff+flagviscoeff*domainCellIndex);
									rho = *(z->PrevFlowsol->Rho+domainCellIndex*flagrho);
									y = *(IfacePtr->deltaP+slabCellIndex);
									u = *(z->Flowsol->u+domainCellIndex);
									
									Face = 2;
									u = CalcmoduPt(z,domainCellIndex,Face);
									
									//This is true for a non moving solid wall boundary
									//if wall has some velocity the u/y must be replace with (u-uwall)/y;
									*(IfacePtr->uTau+slabCellIndex) = sqrt(EffViscosity/rho * u/y);
								}
								break;
							case 5: //KMax face, Top face
								k=NK-2;
								for (j=0;j<NJ-1;j++)
								for (i=0;i<NI-1;i++)
								{
									domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
									slabCellIndex = (NI-1)*j+i;
									
									//EffViscosity = *(z->Flowsol->Tcoeff->EffViscoeff+domainCellIndex);
									EffViscosity = *(z->Flowsol->Viscoeff+flagviscoeff*domainCellIndex);
									rho = *(z->PrevFlowsol->Rho+domainCellIndex*flagrho);
									y = *(IfacePtr->deltaP+slabCellIndex);
									//u = *(z->Flowsol->u+domainCellIndex);
									
									Face = 5;
									u = CalcmoduPt(z,domainCellIndex,Face);
									
									//This is true for a non moving solid wall boundary
									//if wall has some velocity the u/y must be replace with (u-uwall)/y;
									*(IfacePtr->uTau+slabCellIndex) = sqrt(EffViscosity/rho * u/y);
								}
								break;
							case 6: //Kmin face, Bottom face
								k=0;
								for (j=0;j<NJ-1;j++)
								for (i=0;i<NI-1;i++)
								{
									domainCellIndex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
									slabCellIndex = (NI-1)*j+i;
									
									//EffViscosity = *(z->Flowsol->Tcoeff->EffViscoeff+domainCellIndex);
									EffViscosity = *(z->Flowsol->Viscoeff+flagviscoeff*domainCellIndex);
									rho = *(z->PrevFlowsol->Rho+domainCellIndex*flagrho);
									y = *(IfacePtr->deltaP+slabCellIndex);
									//u = *(z->Flowsol->u+domainCellIndex);
									
									Face = 4;
									u = CalcmoduPt(z,domainCellIndex,Face);
									
									//This is true for a non moving solid wall boundary
									//if wall has some velocity the u/y must be replace with (u-uwall)/y;
									*(IfacePtr->uTau+slabCellIndex) = sqrt(EffViscosity/rho * u/y);
								}
								break;
						}
					}
				}
			}//End of flagSolidFluid check
		}//End of Zone looping
	}
	return 0;
}



