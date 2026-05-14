#include "APKTURBCalcMuEff.h"

/**
This function calculates Mu_eff = Mu + constant * Mu_t; constant value is 1 for velocities (momentums in case of Two-Phase flow )and (1/sigmaK) and (1/sigmaEpsilon)
for k and Epsilon respectively.
 */
int CalcMuEff(Base *base, int vartag)
{
	int NI,NJ,NK,nZones,numZones,cellindex;
	int i,j,k,*ptrSize,flagviscoeff,numIface,numZbc,dir;
	int iface,zBc,*ptrRange,index,DiffCoeffVarTag;
	int totalCells,NIcell,NJcell,NKcell,flagr1,flagr2;
	double constant,*domainPtr,*slabIptr,*slabBcptr,dummyViscosity=0.0;
	double F1,mu,mu_t,phi1,phi2,RhoC,RhoD,sigma_f,sigma_r,Ratio,Tp,Te,dia,Kg,Eg,St_t,Kd=0.01,Kp;
	Interface *IfacePtr;
	ZoneBoundary *zBcPtr;
	Zone *zone,*ptrZoneNeighbor;
	FILE *fp;

#ifdef FILE_OPEN
	{
		fp=fopen("Settings/DummyViscocity.dat","r");
		fscanf(fp,"%lf", &dummyViscosity);
		fclose(fp);
	}
#endif
	
	nZones=(base->Number_of_zones);
	if(base->flagTurbModel != 6)	//Other than SST turbulence model
	{
		for(numZones=0;numZones<nZones;numZones++)
		{
			zone = base->CgnsZone+numZones;
			
			if(zone->flagSolidFluid==0)
			{
				ptrSize = zone->Size;
				NI=*(ptrSize+0);
				NJ=*(ptrSize+1);
				NK=*(ptrSize+2);

				if(vartag == 1 || vartag == 2 || vartag == 3||vartag==31||vartag==32||vartag==33)
				{
					constant = 1.0;
					domainPtr = zone->Flowsol->Tcoeff->EffViscoeff;
				}
				else if(vartag == 600)
				{
					switch(base->flagTurbModel)
					{
						case 1:
							constant = 1.0/(zone->Flowsol->Tcoeff->sigmaK);
							break;
						case 2:
							constant = zone->Flowsol->Tcoeff->sigmaK;
							break;
						case 3:
							constant = 1.0/(zone->Flowsol->Tcoeff->sigmaK);
							break;
						case 4:
							constant = 1.0/(zone->Flowsol->Tcoeff->sigmaK);
							break;
						case 5:
							constant = 1.0/(zone->Flowsol->Tcoeff->sigmaK);
							break;
						default:
							printf("Improper turbulence model tag in CalcMuEff.c, vartag=%d\n", vartag);
							exit(78);
							break;
					}
					
					domainPtr = zone->Flowsol->Scalcoeff[0];
				}
				else if(vartag == 601)
				{
					switch(base->flagTurbModel)
					{
						case 1:
							constant = 1.0/(zone->Flowsol->Tcoeff->sigmaEpsilon);
							break;
						case 2:
							constant = zone->Flowsol->Tcoeff->sigmaOmega;
							break;
						case 3:
							constant = 1.0/(zone->Flowsol->Tcoeff->sigmaEpsilon);
							break;
						case 4:
							constant = 1.0/(zone->Flowsol->Tcoeff->sigmaEpsilon);
							break;
						case 5:
							constant = 1.0/(zone->Flowsol->Tcoeff->sigmaEpsilon);
							break;
						default:
							printf("Improper turbulence model tag in CalcMuEff.c, vartag=%d\n", vartag);
							exit(78);
							break;
					}
					
					domainPtr = zone->Flowsol->Scalcoeff[1];
				}
				else if(vartag == 7 || vartag == 14 || vartag == 15)
				{
					constant = 1.0/(zone->Flowsol->Tcoeff->sigmacombvar);
					domainPtr = zone->Flowsol->Tcoeff->EffScalcoeff;
				}
				else if(vartag==41||vartag==42||vartag==43)//chandrasekhar
				{
					constant = 1.0;
					domainPtr=zone->Flowsol->Tcoeff->DummyeffViscoeffofphase[0];
				}
				else if(vartag==44||vartag==45||vartag==46)//chandrasekhar
				{
					constant = 1.0;
					domainPtr=zone->Flowsol->Tcoeff->DummyeffViscoeffofphase[1];
				}
				else if(vartag==35)//chandrasekhar
				{
					constant=1.0/(zone->Flowsol->Tcoeff->sigmaVolfrac);
					domainPtr=zone->Flowsol->Tcoeff->DiffVolfrac[1];
				}
				else if(vartag==62||vartag==63||vartag==64)
				{
					constant=1.0;//1.0/(zone->Flowsol->Tcoeff->sigma_f);//sigma_f (Moukalled et al.)	
					domainPtr=zone->Flowsol->Tcoeff->DummyViscoeff_disperse;
				}
				else if(vartag==65)
				{
					sigma_f=(zone->Flowsol->Tcoeff->sigma_f);
					sigma_r=(zone->Flowsol->Tcoeff->sigma_r);
					constant=1.0/0.7;//1.0/(sigma_f*sigma_r);//(sigma_f*sigma_r)(Moukalled et al.)	
					domainPtr=zone->Flowsol->Tcoeff->DiffVolfrac_disperse;
				}
				flagviscoeff = zone->Flowsol->flagViscoeff;
				mu = *(zone->Flowsol->Viscoeff);
				
				cellindex = 0;
				for(k=0;k<NK-1;k++)
				for(j=0;j<NJ-1;j++)
				for(i=0;i<NI-1;i++)
				{
					if(vartag==41||vartag==42||vartag==43)//chandrasekhar
					{	
						mu_t = *(zone->Flowsol->Tcoeff->TViscoeffofphase[0]+cellindex);
						*(domainPtr+cellindex)=*(zone->Flowsol->ViscoeffDummy[0]+cellindex)+(mu_t* *(zone->Flowsol->VarthetaPhase[0]+cellindex));
						//printf("The vartag is %d and zone is %s\n",vartag,zone->zonename);
						//printf("The value of eff dummy viscosity of the 0 th phase is %lf\n",*(domainPtr+cellindex));
					}
					else if(vartag==44||vartag==45||vartag==46)//chandrasekhar
					{	
						mu_t = *(zone->Flowsol->Tcoeff->TViscoeffofphase[1]+cellindex);
						*(domainPtr+cellindex)=*(zone->Flowsol->ViscoeffDummy[1]+cellindex)+(mu_t* *(zone->Flowsol->VarthetaPhase[1]+cellindex));
						//printf("The vartag is %d and zone is %s\n",vartag,zone->zonename);
						//printf("The value of eff dummy viscosity of the 1st phase is %e\n",*(domainPtr+cellindex));
					}
					else if(vartag==35)
					{
						flagr1=base->BaseParameters->flagRhoPhase[1];
						*(domainPtr+cellindex)= (*(zone->Flowsol->Tcoeff->TViscoeffofphase[1]+cellindex)/(*(zone->Flowsol->RhoofPhase[1]+(cellindex*flagr1))))*constant;
						
						//printf("The vartag is %d and zone is %s\n",vartag,zone->zonename);
						//printf("The value of eff dummy olfrac of the 1st phase is %e\n",*(domainPtr+cellindex));
					}
					else if(vartag==62||vartag==63||vartag==64)//chandrasekhar
					{	
						flagr1=zone->Flowsol->flagRho;//Gas phase 
						flagr2=base->BaseParameters->flagRhoDisperse;//Disperse phase 
						mu_t = *(zone->Flowsol->Tcoeff->TViscoeff+cellindex);
						dia= *(zone->PrevFlowsol->ParticleDia);
						RhoC=(*(zone->Flowsol->Rho+(cellindex*flagr1)));
						RhoD=(*(zone->Flowsol->RhoDisperse+(cellindex*flagr2)));
						Tp=(dia*dia*RhoD)/(18.0*mu);
						Kg=*(zone->Flowsol->Scal[0]+cellindex);
						Eg=*(zone->Flowsol->Scal[1]+cellindex);
						Te=0.125*(Kg/Eg);
						
						St_t=1.0/(1.0+(Tp/Te));
						
						if(Kd>St_t)
							Kp=Kd;
						else
							Kp=St_t;
						
						*(zone->Flowsol->Tcoeff->TViscoeff_disperse+cellindex)=Kp*(mu_t/RhoC);
						//printf("\n *(zone->Flowsol->Tcoeff->TViscoeff_disperse+cellindex)=%e ",*(zone->Flowsol->Tcoeff->TViscoeff_disperse+cellindex));
						*(domainPtr+cellindex)=RhoD* *(zone->Flowsol->VarthetaDisperse+cellindex)* *(zone->Flowsol->Tcoeff->TViscoeff_disperse+cellindex);
						//*(domainPtr+cellindex)=mu_t*constant*(RhoD/RhoC)* *(zone->Flowsol->VarthetaDisperse+cellindex);
						//printf("The vartag is %d and zone is %s\n",vartag,zone->zonename);
						//printf("Zone-%s$The value of eff dummy viscosity of the 1st phase is %e\n",zone->zonename,*(domainPtr+cellindex));
					}
					else if(vartag==65)
					{
						//flagr1=zone->Flowsol->flagRho;
						//mu_t = *(zone->Flowsol->Tcoeff->TViscoeff+cellindex);
						//*(domainPtr+cellindex)=(mu_t /(*(zone->Flowsol->Rho+(cellindex*flagr1))))*constant;
						*(domainPtr+cellindex)=*(zone->Flowsol->Tcoeff->TViscoeff_disperse+cellindex)*constant;
						//printf("The vartag is %d and zone is %s\n",vartag,zone->zonename);
						//printf("The value of eff dummy olfrac of the 1st phase is %e\n",*(domainPtr+cellindex));
					}
					else
					{
						if(flagviscoeff == 1)
						{
							mu   = *(zone->Flowsol->Viscoeff+cellindex);
						}
						mu_t = *(zone->Flowsol->Tcoeff->TViscoeff+cellindex);
						
						*(domainPtr+cellindex) = mu + constant * mu_t + dummyViscosity;
						//printf("\nZone-%s$Vartag-%d$$*(domainPtr+cellindex)=%e",zone->zonename,vartag,*(domainPtr+cellindex));
					}
					cellindex++;
				}
				//printf("The vartag is %d and zone is %s\n",vartag,zone->zonename);
				//printf("%d %d %d Mu_effective Domain = %lf mu=%lf mu_t=%lf dummy=%lf con=%lf\n",k,j,i,*(domainPtr),mu,mu_t,dummyViscosity,constant);
				
				numIface=(zone)->Number_of_Interfaces;
				numZbc=(zone)->Number_of_BC_Patches;

				for(iface=0;iface<numIface;iface++)
				{
					IfacePtr = zone->Iface+iface;
					dir = IfacePtr->Direction;
					
					if(vartag == 1 || vartag == 2 || vartag == 3||vartag==31||vartag==32||vartag==33)
					{
						slabIptr = IfacePtr->ZSlab->EffViscoeff;
					}
					else if(vartag == 600)
					{
						slabIptr=IfacePtr->ZSlab->Scalcoeff[0];
					}
					else if(vartag == 601)
					{
						slabIptr=IfacePtr->ZSlab->Scalcoeff[1];
					}
					else if(vartag == 7 || vartag == 14 || vartag == 15) //Shyam
					{
						slabIptr = IfacePtr->ZSlab->EffScalcoeff;
					}
					else if(vartag==41||vartag==42||vartag==43)
					{
						slabIptr=IfacePtr->ZSlab->DummyeffViscoeffofphase[0];
					}
					else if(vartag==44||vartag==45||vartag==46)
					{
						slabIptr=IfacePtr->ZSlab->DummyeffViscoeffofphase[1];
					}
					else if(vartag==35)
					{
						slabIptr=IfacePtr->ZSlab->DiffVolfrac[1];
						flagr1=base->BaseParameters->flagRhoPhase[1];
					}
					else if(vartag==62||vartag==63||vartag==64)
					{
						slabIptr=IfacePtr->ZSlab->DummyViscoeff_disperse;
					}
					else if(vartag==65)
					{
						slabIptr=IfacePtr->ZSlab->DiffVolfrac_disperse;
					}
					ptrRange=IfacePtr->Point_Range;
					
					switch(dir)
					{
						case 1://Imax face
						case 2://Imin face
							totalCells = (NJ-1) * (NK-1);
							break;
						case 3://Jmax face
						case 4://Jmin face
							totalCells = (NI-1) * (NK-1);
							break;
						case 5://Kmax face
						case 6://Kmin face
							totalCells = (NI-1) * (NJ-1);
							break;
					}
	
					for(index=0; index < totalCells; index++)
					{
						if(vartag==41||vartag==42||vartag==43)
						{
							mu_t = *(IfacePtr->ZSlab->TViscoeffofphase[0]+index);	
							*(slabIptr+index)=*(IfacePtr->ZSlab->ViscoeffDummy[0]+index)+(mu_t* *(IfacePtr->ZSlab->VarthetaPhase[0]+index));	
							//printf("\nDummy Eff(0)=%lf--zone=%s--Iface=%s",*(slabIptr+index),zone->zonename,IfacePtr->Interface_Name);
						}
						else if(vartag==44||vartag==45||vartag==46)
						{
							mu_t = *(IfacePtr->ZSlab->TViscoeffofphase[1]+index);
							*(slabIptr+index)=*(IfacePtr->ZSlab->ViscoeffDummy[1]+index)+(mu_t* *(IfacePtr->ZSlab->VarthetaPhase[1]+index));	
							//printf("\nDummy Eff(1)=%lf--zone=%s--Iface=%s",*(slabIptr+index),zone->zonename,IfacePtr->Interface_Name);
						}
						else if(vartag==35)
						{
							flagr1=base->BaseParameters->flagRhoPhase[1];
							*(slabIptr+index)=(*(IfacePtr->ZSlab->TViscoeffofphase[1]+index)/(*(IfacePtr->ZSlab->RhoofPhase[1]+index*flagr1)))*constant;
						}
						else if(vartag==62||vartag==63||vartag==64)
						{
							flagr1=zone->Flowsol->flagRho;//Gas phase 
							flagr2=base->BaseParameters->flagRhoDisperse;//Disperse phase 
							mu_t = *(IfacePtr->ZSlab->TViscoeff+index);
							dia= *(zone->PrevFlowsol->ParticleDia);
							RhoC=(*(IfacePtr->ZSlab->Rho+index*flagr1));
							RhoD=(*(IfacePtr->ZSlab->RhoDisperse+index*flagr2));
							Tp=(dia*dia*RhoD)/(18.0*mu);
							Kg=*(IfacePtr->ZSlab->Scal[0]+index);
							Eg=*(IfacePtr->ZSlab->Scal[1]+index);
							Te=0.125*(Kg/Eg);
							St_t=1.0/(1.0+(Tp/Te));
						
							if(Kd>St_t)
								Kp=Kd;
							else
								Kp=St_t;
							
							*(IfacePtr->ZSlab->TViscoeff_disperse+index)=Kp*(mu_t/RhoC);
							*(slabIptr+index)=RhoD * *(IfacePtr->ZSlab->VarthetaDisperse+index)* *(IfacePtr->ZSlab->TViscoeff_disperse+index);
							//*(slabIptr+index)=mu_t*	constant *(RhoD/RhoC)* *(IfacePtr->ZSlab->VarthetaDisperse+index);
							//printf("\nDummy Eff_P=%e--zone=%s--Iface=%s",*(slabIptr+index),zone->zonename,IfacePtr->Interface_Name);
						}
						else if(vartag==65)
						{
							//flagr1=zone->Flowsol->flagRho;
							//mu_t = *(IfacePtr->ZSlab->TViscoeff+index);
							//*(slabIptr+index)=(mu_t/(*(IfacePtr->ZSlab->Rho+index*flagr1)))*constant;
							*(slabIptr+index)=*(IfacePtr->ZSlab->TViscoeff_disperse+index)*constant;
							//printf("\nDVolfrac_P=%e--zone=%s--Iface=%s",*(slabIptr+index),zone->zonename,IfacePtr->Interface_Name);
						}
						else 
						{	
							if(flagviscoeff == 1)
								{
									mu   = *(IfacePtr->ZSlab->Viscoeff+index);
								}
								mu_t = *(IfacePtr->ZSlab->TViscoeff+index);
									
								*(slabIptr+index) = mu + constant * mu_t + dummyViscosity;
								//printf("\nzone-%s$$*(slabIptr+index)=%e",zone->zonename,*(slabIptr+index));
								
						}					
					}
					//printf("The vartag is %d and zone is %s\n",vartag,zone->zonename);	
					//printf("dir=%d k=%d j=%d i=%d vartag=%d Mu_effective slab = %e mu=%e constant=%e mu_t=%e\n",dir,k,j,i,vartag,*(slabIptr),mu,constant,mu_t);
				}
				for(zBc=0;zBc<numZbc;zBc++)
				{
					zBcPtr = zone->ZoneBc+zBc;
					dir = zBcPtr->Direction;
					
					if(vartag == 1 || vartag == 2 || vartag == 3||vartag == 31 || vartag == 32 || vartag == 33)
						slabBcptr=zBcPtr->ZSlab->EffViscoeff;
					else if(vartag == 600)
						slabBcptr=zBcPtr->ZSlab->Scalcoeff[0];
					else if(vartag == 601)
						slabBcptr=zBcPtr->ZSlab->Scalcoeff[1];
					else if(vartag == 7 || vartag == 14 || vartag == 15)
						slabBcptr=zBcPtr->ZSlab->EffScalcoeff;
					else if(vartag==41||vartag==42||vartag==43)
						slabBcptr=zBcPtr->ZSlab->DummyeffViscoeffofphase[0];
					else if(vartag==44||vartag==45||vartag==46)
						slabBcptr=zBcPtr->ZSlab->DummyeffViscoeffofphase[1];
					else if(vartag==35)
					{
						slabBcptr=zBcPtr->ZSlab->DiffVolfrac[1];
						flagr1=base->BaseParameters->flagRhoPhase[1];
					}
					else if(vartag==62||vartag==63||vartag==64)
						slabBcptr=zBcPtr->ZSlab->DummyViscoeff_disperse;
					else if(vartag==65)
						slabBcptr=zBcPtr->ZSlab->DiffVolfrac_disperse;
						
					
					
					ptrRange=zBcPtr->Point_Range;
					
					switch(dir)
					{
						case 1://Imax face
						case 2://Imin face
							totalCells = (NJ-1) * (NK-1);
							break;
						case 3://Jmax face
						case 4://Jmin face
							totalCells = (NI-1) * (NK-1);
							break;
						case 5://Kmax face
						case 6://Kmin face
							totalCells = (NI-1) * (NJ-1);
							break;
					}
	
					for(index=0; index < totalCells; index++)
					{
						
						if(vartag==41||vartag==42||vartag==43)
						{
							mu_t = *(zBcPtr->ZSlab->TViscoeffofphase[0]+index);
							*(slabBcptr+index)=*(zBcPtr->ZSlab->ViscoeffDummy[0]+index)+(mu_t* *(zBcPtr->ZSlab->VarthetaPhase[0]+index));
							//printf("\nDummy Eff(0)=%lf--zone=%s--BC=%s",*(slabBcptr+index),zone->zonename,zBcPtr->BCName);
						}
						else if(vartag==44||vartag==45||vartag==46)
						{
							mu_t = *(zBcPtr->ZSlab->TViscoeffofphase[1]+index);
							*(slabBcptr+index)=*(zBcPtr->ZSlab->ViscoeffDummy[1]+index)+(mu_t* *(zBcPtr->ZSlab->VarthetaPhase[1]+index));
							//printf("\nDummy Eff(1)=%e--zone=%s--BC=%s",*(slabBcptr+index),zone->zonename,zBcPtr->BCName);
						}
						else if(vartag==35)
						{
							flagr1=base->BaseParameters->flagRhoPhase[1];
							*(slabBcptr+index)=(*(zBcPtr->ZSlab->TViscoeffofphase[1]+index)/(*(zBcPtr->ZSlab->RhoofPhase[1]+index*flagr1)))*constant;
							//printf("\n*(slabBcptr+index)=%e$Bc-%s",*(slabBcptr+index),zBcPtr->BCName);
						}
						else if(vartag==62||vartag==63||vartag==64)
						{
							flagr1=zone->Flowsol->flagRho;//Gas phase 
							flagr2=base->BaseParameters->flagRhoDisperse;//Disperse phase 
							mu_t = *(zBcPtr->ZSlab->TViscoeff+index);
							dia= *(zone->PrevFlowsol->ParticleDia);
							RhoC=(*(zBcPtr->ZSlab->Rho+index*flagr1));
							RhoD=(*(zBcPtr->ZSlab->RhoDisperse+index*flagr2));
							Tp=(dia*dia*RhoD)/(18.0*mu);
							Kg=*(zBcPtr->ZSlab->Scal[0]+index);
							Eg=*(zBcPtr->ZSlab->Scal[1]+index);
							if((zone->ZoneBc+zBc)->Boundary_Tag == 22)
								Te=*((zone->ZoneBc+zBc)->ZSlab->T_e+index);
							else	
								Te=0.125*(Kg/Eg);
							
							St_t=1.0/(1.0+(Tp/Te));
						
							if(Kd>St_t)
								Kp=Kd;
							else
								Kp=St_t;
							
							*(zBcPtr->ZSlab->TViscoeff_disperse+index)=Kp*(mu_t/RhoC);
							*(slabBcptr+index)=RhoD* *(zBcPtr->ZSlab->VarthetaDisperse+index)**(zBcPtr->ZSlab->TViscoeff_disperse+index);
							//*(slabBcptr+index)=constant * mu_t *(RhoD/RhoC)* *(zBcPtr->ZSlab->VarthetaDisperse+index);
							//printf("\nDummy Eff_P=%e--zone=%s--BC=%s",*(slabBcptr+index),zone->zonename,zBcPtr->BCName);
						}
						else if(vartag==65)
						{
							//flagr1=zone->Flowsol->flagRho;
							//mu_t = *(zBcPtr->ZSlab->TViscoeff+index);//printf("\nmu_t=%e$constant=%e$Density_p=%e",mu_t,constant,(*(zBcPtr->ZSlab->RhoDisperse+index*flagr1)));
							//*(slabBcptr+index)=(mu_t/(*(zBcPtr->ZSlab->Rho+index*flagr1)))*constant;
							*(slabBcptr+index)=*(zBcPtr->ZSlab->TViscoeff_disperse+index)*constant;
							//printf("\nDVolfrac_P=%e--zone=%s--BC=%s",*(slabBcptr+index),zone->zonename,zBcPtr->BCName);
						}
						else
						{
							if(flagviscoeff == 1)
								{
									mu   = *(zBcPtr->ZSlab->Viscoeff+index);
								}
								mu_t = *(zBcPtr->ZSlab->TViscoeff+index);
									
								*(slabBcptr+index) = mu + constant * mu_t + dummyViscosity;
							
							//printf("\nZone-%s$$BC-%s*(slabBcptr+index)=%e",zone->zonename,zBcPtr->BCName,*(slabBcptr+index));
						}
					}
					
					//printf("The vartag is %d and zone is %s\n",vartag,zone->zonename);
					//printf("dir=%d k=%d j=%d i=%d vartag=%d Mu_effective slab = %e mu=%e constant=%e mu_t=%e\n",dir,k,j,i,vartag,*(slabBcptr),mu,constant,mu_t);
				}
			}
		}
	}
	else	//SST turbulence model
	{
		for(numZones=0;numZones<nZones;numZones++)
		{
			zone = base->CgnsZone+numZones;

			if(zone->flagSolidFluid==0)
			{
				ptrSize = zone->Size;
				NI=*(ptrSize+0);
				NJ=*(ptrSize+1);
				NK=*(ptrSize+2);

				flagviscoeff = zone->Flowsol->flagViscoeff;

				if(vartag == 1 || vartag == 2 || vartag == 3)
				{
					constant = 1.0;
					domainPtr = zone->Flowsol->Tcoeff->EffViscoeff;

					for(k=0;k<NK-1;k++)
						for(j=0;j<NJ-1;j++)
							for(i=0;i<NI-1;i++)
					{
						cellindex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
						mu   = *(zone->Flowsol->Viscoeff+flagviscoeff*cellindex);
						mu_t = *(zone->Flowsol->Tcoeff->TViscoeff+cellindex);
						
						*(domainPtr+cellindex) = mu + constant * mu_t + dummyViscosity;
						
						//printf("%d %d %d Mu_effective Domain = %e\n",k,j,i,*(domainPtr+cellindex));
					}
				}
				else
				{
					if(vartag == 600)
					{
						domainPtr = zone->Flowsol->Scalcoeff[0];
						phi1 = zone->Flowsol->Tcoeff->SigmaK1;
						phi2 = zone->Flowsol->Tcoeff->SigmaK2;
					}
					else if(vartag == 601)
					{
						domainPtr = zone->Flowsol->Scalcoeff[1];
						phi1 = zone->Flowsol->Tcoeff->SigmaOmega1;
						phi2 = zone->Flowsol->Tcoeff->SigmaOmega2;
					}

					for(k=0;k<NK-1;k++)
					for(j=0;j<NJ-1;j++)
					for(i=0;i<NI-1;i++)
					{
						cellindex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
						F1 = *(zone->Flowsol->Tcoeff->F1+cellindex);
						constant = F1 * phi1 + (1.0-F1) * phi2;
						mu   = *(zone->Flowsol->Viscoeff+flagviscoeff*cellindex);
						mu_t = *(zone->Flowsol->Tcoeff->TViscoeff+cellindex);
						
						*(domainPtr+cellindex) = mu + constant * mu_t + dummyViscosity;
						
						//printf("%d %d %d Mu_effective Domain = %e\n",k,j,i,*(domainPtr+cellindex));
					}
				}

				numIface=(zone)->Number_of_Interfaces;
				numZbc=(zone)->Number_of_BC_Patches;

				for(zBc=0;zBc<numZbc;zBc++)
				{
					zBcPtr = zone->ZoneBc+zBc;
					dir = zBcPtr->Direction;
					
					if(vartag == 1 || vartag == 2 || vartag == 3)
					{
						constant = 1.0;
						slabBcptr=zBcPtr->ZSlab->EffViscoeff;
					}

					else if(vartag == 600)
					{
						slabBcptr=zBcPtr->ZSlab->Scalcoeff[0];
						if(zBcPtr->Boundary_Tag == 22)
							constant = zone->Flowsol->Tcoeff->SigmaK1;
						else
							constant = zone->Flowsol->Tcoeff->SigmaK2;
					}

					else if(vartag == 601)
					{
						slabBcptr=zBcPtr->ZSlab->Scalcoeff[1];
						if(zBcPtr->Boundary_Tag == 22)
							constant = zone->Flowsol->Tcoeff->SigmaOmega1;
						else
							constant = zone->Flowsol->Tcoeff->SigmaOmega2;
					}
					
					ptrRange=zBcPtr->Point_Range;
					
					switch(dir)
					{
						case 1://Imax face
						case 2://Imin face
							totalCells = (NJ-1) * (NK-1);
							break;
						case 3://Jmax face
						case 4://Jmin face
							totalCells = (NI-1) * (NK-1);
							break;
						case 5://Kmax face
						case 6://Kmin face
							totalCells = (NI-1) * (NJ-1);
							break;
					}
					
					for(index=0; index < totalCells; index++)
					{
						mu   = *(zBcPtr->ZSlab->Viscoeff+flagviscoeff*index);
						mu_t = *(zBcPtr->ZSlab->TViscoeff+index);
							
						*(slabBcptr+index) = mu + constant * mu_t + dummyViscosity;
						//printf("dir=%d k=%d j=%d i=%d vartag=%d Mu_effective slab = %e mu=%e constant=%e mu_t=%e\n",dir,k,j,i,vartag,*(slabBcptr+index),mu,constant,mu_t);
					}
				}
			}
		}
		//Updating Turbulent Viscosities in the Interfaces of all Zones
		for(numZones=0;numZones<nZones;numZones++)
		{
			zone = base->CgnsZone+numZones;

			if(zone->flagSolidFluid==0)
			{
				numIface = zone->Number_of_Interfaces;

				for(iface=0;iface<numIface;iface++)
				{
					IfacePtr = zone->Iface+iface;
					ptrZoneNeighbor = (Zone *)FindZonePtr(IfacePtr->Adjoining_Zone_Name,base);

					if(vartag == 1 || vartag == 2 || vartag == 3)
						DiffCoeffVarTag = 714;
					else if(vartag == 600)
						DiffCoeffVarTag = 700;
					else if(vartag == 601)
						DiffCoeffVarTag = 701;

					if(ptrZoneNeighbor->flagSolidFluid == 0)
					{
						FillerInterfaceVariables(DiffCoeffVarTag,IfacePtr,ptrZoneNeighbor,(IfacePtr->Direction)-1);
					}
					else
					{
						//If the adjacent zone is SOLID, boundary condition has to be applied for turbulent viscosity
						//i.e., Turbulent Viscosity can be made equal to 'Zero'.
						if(base->flagTurbModel==1 || base->flagTurbModel==3)	//Check the Interface Condition for "log version of standard k-epsilon model"
						{
							//printf("Updating turbulent viscosity in Zone: %s \n", zone->zonename);
							HomogeneousNeumannIC(zone,IfacePtr,DiffCoeffVarTag);//NonOrthogonalHomogeneousNeumannIC(base,zone,IfacePtr,715);//HomogeneousNeumannIC(zone,IfacePtr,715);
							//printf("Updation is over ...\n");
						}
						else
						{
							ptrRange = IfacePtr->Point_Range;
							
							NIcell = *(ptrRange+3) - *(ptrRange+0);
							NJcell = *(ptrRange+4) - *(ptrRange+1);
							NKcell = *(ptrRange+5) - *(ptrRange+2);
							
							slabIptr = FindSlabPtrIface(IfacePtr,DiffCoeffVarTag);
							flagviscoeff = zone->Flowsol->flagViscoeff;
							
							switch(IfacePtr->Direction)
							{
								case 1://Imax face
								case 2://Imin face
									totalCells = NJcell * NKcell;
									break;
								case 3://Jmax face
								case 4://Jmin face
									totalCells = NIcell * NKcell;
									break;
								case 5://Kmax face
								case 6://Kmin face
									totalCells = NIcell * NJcell;
									break;
							}
							
							mu = *(IfacePtr->ZSlab->Viscoeff);
							for(index=0; index < totalCells; index++)
							{
								if(flagviscoeff == 1)
								{
									mu = *(IfacePtr->ZSlab->Viscoeff+index);
								}
								*(slabIptr+index) = mu + dummyViscosity;//Since on viscous wall: mu_t = 0.0;
							}
						}
					}
				}
			}
		}
	}
	return 0;
}
