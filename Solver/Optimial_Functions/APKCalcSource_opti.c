#include "APKCalcSource_opti.h"
/**
The function CalcSource_Opti() is used to calculate the various source terms occuring in the governing equations of various variables. This function, hence, takes the variable tag as input and depending on the variable tag, either computes or calls necessary functions to calculate the source terms corresponding to that particular variable.  This is function is written more efficiently when compared to earlier versions of this function.
 */


double CalcSource_opti(Base *base,int variableTag)
{
	int NI,NJ,NK,k,j,i,*ptrSize,indexP;
	double diagonalCoef=0.0;
	double *domain_prev,*variable,*variable_DIAG,sumMassFlux,velP,*massFlux;
	Zone *zone;

	int nZones,numZones,temp,flag;
	
	///This is for storing the gradient location according to the variable direction in case of u,v,w variables
	int locate_grad,visc_tag=0,flagrho,index[3];
	
	///For Multiphase variables.
	int temp1,flagLimProduction;
	double *BETA,r0,*u,c,con,*B,*drag,*Lift;
	double *DiffFlux,*Volfrac,*presgrad,rho,*Rho,*r,*r1;
	
	///adhi+shash for disperse phase flow
	double *Drag,*MomentumD,*MomentumG,Ud,Vd,Wd;
	double t_p,t_e,t_L;
	double K_Source;
	int flagspecificheat,flagCPDrop,flagDia;//shash
	double dropmass;//shash
	double Tsurf,Cpv;//shash
	//shash
	double EPSILON;
	EPSILON = base->BaseParameters->EPSILONDrop;
	//mrunalini
	double VolumeFrac;
	int flagrhod;
	///All pointers are related to Buoyancy source term.
	double Tref,beta,gravity,*ptrVolume,*ptrTemp; 
	
 	int flag_energysource; ///Parameter related to energy source

	///Turbulence variables
	double Kold,Epsold,sumGradSquare,MuEffK,Dissipation,Gamma,betastar,gamma,beta1;
	double mue_t,TProduction,OmegaOld,delt,vol;
	double F1,Gamma1,Gamma2,Beta1,Beta2,CrossDiffusion;
	double C0epsi,C1epsi,C2epsi,MuEffEps,K,S,ProductionMaxLimit;
	
	///Prarameters realted to Solidification and melting problems and HEM
	double P_o,dele,dele_prev,energySource,*prevfluxes;
	double C,liq_frac,A,*CP,C_p;
	
	///For ESP Problems
	double DiElec;

	///Parameters realted to Turbulent combustion module
	double MinFrac,Stoich,FMassFrac,OxMassFrac,PMassFrac,MinProd,KE,epsilon;

	///For LMMHD Problems
	double B_x, B_y, B_z,sigma,magperm,magdiffu,*poten_source;
	int flagelectcond,flagpermeability,flagMagdiff;
	

	
	if (variableTag==1 || variableTag==31 ||variableTag==100|| variableTag==3100|| variableTag==62 || variableTag==6200)
		locate_grad=0;
	
	if (variableTag==2 || variableTag==32 ||variableTag==200|| variableTag==3200|| variableTag==63 || variableTag==6300)
		locate_grad=1;
	
	if (variableTag==3 || variableTag==33 ||variableTag==300|| variableTag==3300|| variableTag==64|| variableTag==6400)
		locate_grad=2;
	
	if(base->BaseParameters->flagBuoyancy==1)
	{
		Tref=base->BaseParameters->TReference;
		beta=base->BaseParameters->VolumetricExpCoeff;
	}
	if(base->BaseParameters->flagGravity==1)
	{
		if (locate_grad==0)
			gravity = base->BaseParameters->gravityX;
		if (locate_grad==1)
			gravity = base->BaseParameters->gravityY;
		if (locate_grad==2)
			gravity = base->BaseParameters->gravityZ;
	}
	
	if(variableTag>=41&&variableTag<=55)
	{
		temp = (variableTag-41)%3;
		if(temp==0)
			gravity = base->BaseParameters->gravityX;
		if(temp==1)
			gravity = base->BaseParameters->gravityY;
		if(temp==2)	
			gravity = base->BaseParameters->gravityZ;
	}
	
	C = base->BaseParameters->C;
	
	nZones = base->Number_of_zones;
	
	if(base->flagTwoPhase==1)
		P_o = base->BaseParameters->P_o;//reference pressure to handle {(P/Rho).(del.U)} term
	
	flagLimProduction = base->BaseParameters->Turb->flagLimProduction;
	base->BaseParameters->Turb->TProductionMaxLimitCount = 0;
	
	delt = base->BaseParameters->TimeStep;
	
	if((variableTag==600 || variableTag==601))
		delt=base->BaseParameters->Turb->DualTimeStep;//delt=delt/10.0;

	for(numZones=0;numZones<nZones;numZones++)
	{
		zone = (base->CgnsZone+numZones);
		ptrSize= zone->Size;
		
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		
		variable = zone->Flowsol->Source_term;
		variable_DIAG = zone->Flowsol->Diag_Source;
		if(base->flagMPFV==0)
			domain_prev = MakeDomainArray(zone,variableTag*10);
		
		if(base->flagSolidification==1)
		{
			CP =zone->Flowsol->SpecificHeat;
			ptrTemp=zone->Flowsol->enthalpy;
		}
		else
		{
			ptrTemp=zone->Flowsol->T;
			C_p = 1.0;
		}
		
		ptrVolume= zone->ZGeom->Volume;
		flagelectcond=zone->Flowsol->flageleccond;
		flagrho = zone->Flowsol->flagRho;
		
		flagspecificheat = zone->Flowsol->flagSpecificHeat;//shash
		flagCPDrop = zone->Flowsol->flagCPdrop;//shash
		flagDia=base->BaseParameters->flagVariableDia; //Mrunalini
		flagrhod = base->BaseParameters->flagRhoDisperse;
		
		if(base->BaseParameters->TimeStepingSchemeTag==1)
		{
			if (locate_grad==0)
				prevfluxes = zone->Flowsol->PrevTStepFluxes_u;
			if (locate_grad==1)
				prevfluxes = zone->Flowsol->PrevTStepFluxes_v;
			if (locate_grad==2)
				prevfluxes = zone->Flowsol->PrevTStepFluxes_w;
			if(variableTag==7 || variableTag==6)
				prevfluxes = zone->Flowsol->PrevTStepFluxes_T;
		}
		
		//adhi+shash
		if((base->flagParticle==1)|| (base->flagDroplet==1))
		{
			if (locate_grad==0)
			{
				Drag = zone->DragX ;
				MomentumD = zone->PrevFlowsol->U_Disperse; //Changed by mrunalini
				MomentumG = zone->PrevFlowsol->u;
			}
			if (locate_grad==1)
			{
				Drag = zone->DragY ;
				MomentumD = zone->PrevFlowsol->V_Disperse; //Changed by mrunalini
				MomentumG = zone->PrevFlowsol->v;
			} 
			if (locate_grad==2)
			{
				Drag = zone->DragZ ;
				MomentumD = zone->PrevFlowsol->W_Disperse; //Changed by mrunalini
				MomentumG = zone->PrevFlowsol->w;
			}
			
			
		}
		
		if(zone->flagHeatGeneration==1 && zone->flagHeat_con_var==1)
			flag_energysource=1;
		else
			flag_energysource=0;
		
		if(!(base->flagTurb == 0 && zone->Flowsol->flagViscoeff == 0))
			visc_tag=1;
		
		if(variableTag==1 || variableTag==100 || variableTag==2 || variableTag==200 || variableTag==3 || variableTag==300)
		{
			indexP=0;
			if(zone->flagSolidFluid == 0)
			{
				for(k=0;k<NK-1;k++)
				for(j=0;j<NJ-1;j++)
				for(i=0;i<NI-1;i++)
				{
					*(variable+indexP)=0.0;//printf("\n It is %d",variableTag);
					if(variableTag<4)
						*(variable+indexP)=*(zone->Flowsol->Pressure_Grad[variableTag-1]+indexP);
					
					massFlux = zone->Flowsol->Massflux+indexP*6;
					velP = *(domain_prev+indexP);
					sumMassFlux = *(massFlux)+*(massFlux+1)+*(massFlux+2)+*(massFlux+3)+*(massFlux+4)+*(massFlux+5);
					
 					*(variable+indexP) = *(variable+indexP) + velP*sumMassFlux;
					*(variable_DIAG+indexP) = 0.0;
					if(variableTag==1 || variableTag==100)
					  // *(variable_DIAG+indexP) = *(variable_DIAG+indexP) + (base->BaseParameters->u_source)*(*(zone->ZGeom->Volume+indexP));
					    *(variable+indexP) = *(variable+indexP) - velP*(base->BaseParameters->u_source)*(*(zone->ZGeom->Volume+indexP));
				
					if(variableTag==2 || variableTag==200)
					  // *(variable_DIAG+indexP) = *(variable_DIAG+indexP) + (base->BaseParameters->v_source)*(*(zone->ZGeom->Volume+indexP));
					   *(variable+indexP) = *(variable+indexP) - velP*(base->BaseParameters->v_source)*(*(zone->ZGeom->Volume+indexP));
         			
					rho = *(zone->PrevFlowsol->Rho+indexP*flagrho);
					
					if(base->flagSolidification==1)
					{
						liq_frac = *(zone->Flowsol->f+indexP);
						A  = -C * (1.0 - liq_frac) * (1.0-liq_frac);
						A  = A / ((liq_frac * liq_frac * liq_frac) + 0.001);
						A  = A/rho;
						*(variable_DIAG+indexP) = *(variable_DIAG+indexP) - A;
						C_p = *(CP+indexP);
					}
					
					if(base->flagEfield==1)
						*(variable+indexP) = *(variable+indexP) +*(zone->Flowsol->ESP_Source[locate_grad]+indexP);
					
					if(base->flagLMMHD==1)
					{
						
												
						if(base->flaginduced==1)
						{
							*(variable+indexP) = *(variable+indexP)+*(zone->Flowsol->LMMHD_Source[locate_grad]+indexP);	
						}
						else
						{
                                               	sigma = *(zone->PrevFlowsol->Eleconductivity+indexP*flagelectcond);
						rho = *(zone->Flowsol->Rho+indexP*flagrho);

							B_x = *(zone->Flowsol->MfieldX+indexP);
							B_y = *(zone->Flowsol->MfieldY+indexP);
							B_z = *(zone->Flowsol->MfieldZ+indexP);

		
							//	if(variableTag==100 ||variableTag==200 || variableTag==300)

							*(variable+indexP)=*(variable+indexP)+*(zone->Flowsol->LMMHDmass_Source[locate_grad]+indexP);

						//	else if (variableTag==1 ||variableTag==2 || variableTag==3)						
						//	*(variable+indexP) = *(variable+indexP)+*(zone->Flowsol->LMMHD_Source[locate_grad]+indexP);	

						//	if 	(variableTag==1)
						//	*(variable_DIAG+indexP) =  *(variable_DIAG+indexP)+*(zone->ZGeom->Volume+indexP)*(sigma/rho) *(B_y*B_y+B_z*B_z) ;
						//	else if (variableTag==2)
						//	*(variable_DIAG+indexP) = *(variable_DIAG+indexP)+*(zone->ZGeom->Volume+indexP)*(sigma/rho) *(B_x*B_x+B_z*B_z) ;
						//	else if (variableTag==3)
						//	*(variable_DIAG+indexP) = *(variable_DIAG+indexP)+*(zone->ZGeom->Volume+indexP)*(sigma/rho) *(B_x*B_x+B_y*B_y) ;

						}
					

					
						
					}
					if((base->flagParticle==1)|| (base->flagDroplet==1))
					{
						if((*(zone->PrevFlowsol->VarthetaDisperse+indexP))<=EPSILON)
						{
							*(variable+indexP)=*(variable+indexP) + 0.0;
						}
						else
						{
							
							VolumeFrac = *(zone->PrevFlowsol->VarthetaDisperse+indexP)*base->BaseParameters->InitialVolumeFraction;
							*(variable+indexP)=*(variable+indexP) - (1.0/(*(zone->PrevFlowsol->Rho+indexP*flagrho)*(1-VolumeFrac)))*(*(ptrVolume+indexP))*VolumeFrac*(*(Drag+indexP)); //
						}
					}
					if((base->BaseParameters->flagVariableDia == 1)&&(base->flagDroplet==1))
					{  if((*(zone->PrevFlowsol->VarthetaDisperse+indexP))<=EPSILON)
						{
							*(variable+indexP)=*(variable+indexP) + 0.0;
							
						}
						else
						{
							VolumeFrac = *(zone->PrevFlowsol->VarthetaDisperse+indexP)*base->BaseParameters->InitialVolumeFraction;
							*(variable_DIAG+indexP) = *(variable_DIAG+indexP) + (1.0/(*(zone->PrevFlowsol->Rho+indexP*flagrho)*(1-VolumeFrac)))*(*(ptrVolume+indexP))*(*(zone->PrevFlowsol->Evaporate+indexP)); //
							*(variable+indexP) = *(variable+indexP)+ (1.0/(*(zone->PrevFlowsol->Rho+indexP*flagrho)*(1-VolumeFrac)))*(*(ptrVolume+indexP))*(*(zone->PrevFlowsol->Evaporate+indexP))*(*(MomentumD+indexP)); //
						}
					}	
					if(base->BaseParameters->flagGravity==1)
					{ 
						if(base->BaseParameters->flagBuoyancy==1)
						{
							*(variable+indexP)=*(variable+indexP)-(*(ptrVolume+indexP))*beta*gravity*((*(ptrTemp+indexP)/C_p) - Tref);
						}
						else
							*(variable+indexP) = *(variable+indexP) + gravity * (*(ptrVolume+indexP));
					}
					
					if(base->BaseParameters->TimeStepingSchemeTag==1)
					{
						*(variable+indexP) = *(variable+indexP) + *(prevfluxes+indexP);
					}
					
					if(!(base->flagTurb == 0 && zone->Flowsol->flagViscoeff == 0))
					{  
						index[0]= k ; index[1]=j; index[2]=i;
						*(variable+indexP) = *(variable+indexP) + (*(ptrVolume+indexP)/rho)*Findgradient(base,zone,index,variableTag);
					}
					indexP++;
				}
			}
			else
			{
				for(k=0;k<NK-1;k++)
				for(j=0;j<NJ-1;j++)
				for(i=0;i<NI-1;i++)
				{
					*(variable+indexP)=0.0;
					*(variable_DIAG+indexP) = 0.0;
					indexP++;
				}
			}

		}
		else if (variableTag==6)
		{
			indexP=0;
			
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				*(variable+indexP)=0.0;
				massFlux = zone->Flowsol->Massflux+indexP*6;
				velP = *(domain_prev+indexP);
				sumMassFlux = *(massFlux)+*(massFlux+1)+*(massFlux+2)+*(massFlux+3)+*(massFlux+4)+*(massFlux+5);
			
				//*(variable+indexP) = velP*sumMassFlux*(*(zone->PrevFlowsol->SpecificHeat+indexP*flagspecificheat));

				*(variable_DIAG+indexP) = 0.0;
				rho = *(zone->PrevFlowsol->Rho+indexP*flagrho);
				if(base->BaseParameters->flagEnergySource==1 && zone->flagHeatGeneration == 1)
				{
					energySource = *((zone->Flowsol)->Energysource+indexP*flag_energysource);
					*(variable+indexP)=*(variable+indexP) + (1.0/rho)*energySource*(*(ptrVolume+indexP));
				}

				if(base->flagLMMHD==1)
				  {
				    if(base->BaseParameters->flagJouledissipation==1){                      
				                {
                                   sigma = *(zone->PrevFlowsol->Eleconductivity+indexP*flagelectcond);
                                       if(sigma != 0)
                                          {
								                        *(variable+indexP)=*(variable+indexP) + *(zone->ZGeom->Volume+indexP)*(1./(sigma*rho)) *(*(zone->Flowsol-> ElecurrentdensityX+indexP)*(*(zone->Flowsol-> ElecurrentdensityX+indexP))
                                                                                              +	*(zone->Flowsol-> ElecurrentdensityY+indexP)*(*(zone->Flowsol-> ElecurrentdensityY+indexP))
											                                                  + *(zone->Flowsol-> ElecurrentdensityZ+indexP)*(*(zone->Flowsol-> ElecurrentdensityZ+indexP)))    ;
                                           }      
                                 }
				    }
                  }


				if(base->BaseParameters->TimeStepingSchemeTag==1)
					*(variable+indexP)=*(variable+indexP)+(*(prevfluxes+indexP));
				
				
				
				if((base->BaseParameters->flagVariableDia==1)&&(base->flagDroplet==1))
		  		{
         	  			if((*(zone->PrevFlowsol->VarthetaDisperse+indexP))<=EPSILON)
		   			{
						*(variable+indexP)=*(variable+indexP) + 0.0;
					}
		  		  	else
				   	{
						VolumeFrac = *(zone->PrevFlowsol->VarthetaDisperse+indexP)*base->BaseParameters->InitialVolumeFraction;

				    		Tsurf = *(zone->Flowsol->Tdrop+indexP);			
				    		Cpv = base->BaseParameters->CPvapour;
						
						*(variable+indexP) = *(variable+indexP) + (1.0/((*(zone->PrevFlowsol->Rho+indexP*flagrho))*(1-VolumeFrac)))*(*(ptrVolume+indexP))*(*(zone->PrevFlowsol->Evaporate+indexP))*(Cpv*Tsurf);	
				    //
				   		*(variable+indexP) = *(variable+indexP) -(1.0/((*(zone->PrevFlowsol->Rho+indexP*flagrho))*(1-VolumeFrac)))*(*(ptrVolume+indexP)*(*(zone->PrevFlowsol->QDroptoGas+indexP))); //
						*(variable_DIAG+indexP) = *(variable_DIAG+indexP)+ (1.0/(*(zone->PrevFlowsol->Rho+indexP*flagrho)*(1-VolumeFrac)))*(*(ptrVolume+indexP))*(*(zone->PrevFlowsol->Evaporate+indexP))*Cpv; //(*(zone->PrevFlowsol->SpecificHeat+indexP*flagspecificheat)); //
					}
				   	
				  }
				indexP++;
			}
		}
		else if (variableTag==7)
		{
			indexP=0;
	
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				massFlux = zone->Flowsol->Massflux+indexP*6;
				velP = *(domain_prev+indexP);
				sumMassFlux = *(massFlux)+*(massFlux+1)+*(massFlux+2)+*(massFlux+3)+*(massFlux+4)+*(massFlux+5);
			
				*(variable+indexP) = velP*sumMassFlux;

				*(variable_DIAG+indexP) = 0.0;
				
				rho = *(zone->PrevFlowsol->Rho+indexP*flagrho);
				if(base->flagTwoPhase==1)
					*(variable+indexP) = *(variable+indexP) -(P_o/rho)*sumMassFlux;


				if(base->BaseParameters->flagEnergySource==1 && zone->flagHeatGeneration == 1)
				{
					energySource = *((zone->Flowsol)->Energysource+indexP*flag_energysource);
					*(variable+indexP)=*(variable+indexP) + (1.0/rho)*energySource*(*(ptrVolume+indexP));
				}
				
					
				if(base->flagSolidification==1)
				{
					dele 	 = *(zone->Flowsol->Latent_Heat_Content+indexP);
					dele_prev= *(zone->PrevFlowsol->Latent_Heat_Content+indexP);
					*(variable+indexP)  = *(variable+indexP) - ((dele - dele_prev)*(*(ptrVolume+indexP))/delt);
					index[0]= k ; index[1]=j; index[2]=i;
					if(base->BaseParameters->flag_mushy==1)
						*(variable+indexP)=*(variable+indexP)-(Find_Mushy_Source(base,zone,indexP,index)/delt);
				}
				
				if(base->BaseParameters->TimeStepingSchemeTag==1)
				*(variable+indexP)=*(variable+indexP)+(*(prevfluxes+indexP));
				
				indexP++;
			}
		}
		else if (variableTag==8)
		{
			indexP=0;
			DiElec = *(zone->Flowsol->DiElecConst);
			for(k=0;k<NK-1;k++)
				for(j=0;j<NJ-1;j++)
					for(i=0;i<NI-1;i++)
			{
				*(variable+indexP)=*(zone->Flowsol->RhoIon+indexP) * (*(zone->ZGeom->Volume+indexP))/DiElec;
				*(variable_DIAG+indexP) = 0.0;
				indexP++;
			}
		}
		else if (variableTag==9)
		{
			indexP=0;
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				*(variable+indexP)=0.0;
				*(variable_DIAG+indexP) = 0.0;
				indexP++;
			}
		}
		else if (variableTag==17)
		{
			indexP=0;
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			  {	*(variable+indexP)=0.0;
				poten_source      =zone->Flowsol->PotentialMHDpoissonRHS+indexP*6;
				*(variable+indexP)  = *(variable+indexP) + *(poten_source)+*(poten_source+1)+*(poten_source+2)+*(poten_source+3)+*(poten_source+4)+*(poten_source+5);
				//*(variable+indexP) = *(zone->Flowsol->PotentialMHDpoissonRHS+indexP);
				
				*(variable_DIAG+indexP) = 0.0;
				indexP++;
			}
		}
		else if (variableTag==56 )
		{
			indexP=0;
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				*(variable+indexP)=0.0;

				*(variable+indexP) = *(variable+indexP)+*(zone->Flowsol->LMMHD_Source[0]+indexP);
				
				*(variable_DIAG+indexP) = 0.0;
				indexP++;
			}
		}
		else if (variableTag==57)
		{
			indexP=0;
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				*(variable+indexP)=0.0;
				*(variable+indexP) = *(variable+indexP)+*(zone->Flowsol->LMMHD_Source[1]+indexP);
				
				*(variable_DIAG+indexP) = 0.0;
				indexP++;
			}
		}
		else if (variableTag==58)
		{
			indexP=0;
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				
				*(variable+indexP)=0.0;
				*(variable+indexP) = *(variable+indexP)+*(zone->Flowsol->LMMHD_Source[2]+indexP);
				// 				
				*(variable_DIAG+indexP) = 0.0;
				indexP++;
			}
		}
		else if (variableTag==14)
		{
			indexP=0;
	
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				massFlux = zone->Flowsol->Massflux+indexP*6;
				velP = *(domain_prev+indexP);
				sumMassFlux = *(massFlux)+*(massFlux+1)+*(massFlux+2)+*(massFlux+3)+*(massFlux+4)+*(massFlux+5);
			
				*(variable+indexP) = velP*sumMassFlux;

				*(variable_DIAG+indexP) = 0.0;
				indexP++;
			}
		}
		else if (variableTag==15)
		{
			indexP=0;
	
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				massFlux = zone->Flowsol->Massflux+indexP*6;
				velP = *(domain_prev+indexP);
				sumMassFlux = *(massFlux)+*(massFlux+1)+*(massFlux+2)+*(massFlux+3)+*(massFlux+4)+*(massFlux+5);
			
				*(variable+indexP) = velP*sumMassFlux;

				*(variable_DIAG+indexP) = 0.0;
				
				
				MinFrac = 0.0;
				Stoich =     *(base->BaseParameters->Stoich);
				FMassFrac =  *(zone->PrevFlowsol->FuelMassFrac+indexP);
				OxMassFrac=  *(zone->PrevFlowsol->OxiMassFrac +indexP);
				PMassFrac =  *(zone->PrevFlowsol->ProdMassFrac+indexP);
				MinProd = 0.5*PMassFrac/(1.0+Stoich);
				
				MinFrac = FMassFrac;
				if(MinFrac > (OxMassFrac/Stoich))
					MinFrac = OxMassFrac/Stoich;

				if(MinFrac > MinProd)
					MinFrac = MinProd;
				
				KE = 	  *(zone->PrevFlowsol->Scal[0]+indexP);
				epsilon = *(zone->PrevFlowsol->Scal[1]+indexP);

				*(variable+indexP) = *(variable+indexP) - 4.0*(epsilon/(KE))*MinFrac*(*(ptrVolume+indexP));
				indexP++;
			}
		}
		else if(variableTag==31 || variableTag==32 || variableTag==33||variableTag==3100 || variableTag==3200 || variableTag==3300)
		{	
			Rho = zone->PrevFlowsol->Rho;
			r = (zone->Flowsol)->RhoofPhase[0];
			flag = base->BaseParameters->flagRhoPhase[0];
			indexP=0;
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				*(variable+indexP)=0.0;
				if(variableTag==31)
				{
					*(variable+indexP)=*(zone->Flowsol->Pressure_Grad[0]+indexP);
					//printf("\n*(variable+indexP)=%e ",*(zone->Flowsol->Pressure_Grad[0]+indexP));
					if(base->BaseParameters->flagGravity==1 && base->flagVariableDensityProblem==0)
					{	
						if(base->flagTurb==0)
						{
							*(variable+indexP)=*(variable+indexP)+((*(ptrVolume+indexP))*(*(r+indexP*flag)-*(Rho+indexP*flagrho)))*gravity;
						}
						else if(base->flagTurb==1)
						{
							*(variable+indexP)=*(variable+indexP)+((*(ptrVolume+indexP))*(*(r+indexP*flag)-*(Rho+indexP*flagrho)))*gravity;
							//*(variable+indexP)=*(variable+indexP)- *(Rho+indexP)* gravity * (*(ptrVolume+indexP));
//((*(ptrVolume+indexP))*((*(r+indexP*flag)* *(zone->Flowsol->VarthetaPhase[1]+indexP))-*(Rho+indexP*flagrho)))*gravity;//body force only on bubble
							//printf("\nVartag=%d",variableTag);
							//printf("\n The body force term is %lf",((*(ptrVolume+indexP))*((*(r+indexP*flag)* *(zone->Flowsol->VarthetaPhase[1]+indexP))-*(Rho+indexP*flagrho)))*gravity);
						}
					}
					if(base->flagTurb==1)//Added by chandrasekhar --- Turbulent MPFV
					{
						index[0]= k ; index[1]=j; index[2]=i;
						*(variable+indexP) = *(variable+indexP) + (*(ptrVolume+indexP))*Findgradient(base,zone,index,variableTag);
						//*(variable+indexP) = *(variable+indexP)+KSource(base,zone,index,variableTag);
						//printf("[Vartag=%d]$[zone %s]$[k=%d$j=%d$i=%d]$Source=%lf\n",variableTag,zone->zonename,k,j,i,*(variable+indexP));
					}
				}
				else if(variableTag==32)
				{	
					*(variable+indexP)=*(zone->Flowsol->Pressure_Grad[1]+indexP);
					if(base->BaseParameters->flagGravity==1 && base->flagVariableDensityProblem==0)
					{	
						
						if(base->flagTurb==0)
						{
							*(variable+indexP)=*(variable+indexP)+((*(ptrVolume+indexP))*(*(r+indexP*flag)-*(Rho+indexP*flagrho)))*gravity;
						}
						else if(base->flagTurb==1)
						{
							*(variable+indexP)=*(variable+indexP)+((*(ptrVolume+indexP))*(*(r+indexP*flag)-*(Rho+indexP*flagrho)))*gravity;
							//*(variable+indexP)=*(variable+indexP)- *(Rho+indexP)* gravity * (*(ptrVolume+indexP));
							//printf("\nVartag=%d",variableTag);
							//printf("\n The body force term is %lf",((*(ptrVolume+indexP))*((*(r+indexP*flag)* *(zone->Flowsol->VarthetaPhase[1]+indexP))-*(Rho+indexP*flagrho)))*gravity);
						}
					}
					if(base->flagTurb==1)//Added by chandrasekhar --- Turbulent MPFV
					{
						index[0]= k ; index[1]=j; index[2]=i;
						*(variable+indexP) = *(variable+indexP) + (*(ptrVolume+indexP))*Findgradient(base,zone,index,variableTag);
						//*(variable+indexP) = *(variable+indexP) + KSource(base,zone,index,variableTag);
						//printf("[Vartag=%d]$[zone %s]$[k=%d$j=%d$i=%d]$Source=%lf\n",variableTag,zone->zonename,k,j,i,*(variable+indexP));
					}
				}
				else if(variableTag==33)                                                                {	
					*(variable+indexP)=*(zone->Flowsol->Pressure_Grad[2]+indexP);
					if(base->BaseParameters->flagGravity==1 && base->flagVariableDensityProblem==0)
					{	
						if(base->flagTurb==0)
						{
							*(variable+indexP)=*(variable+indexP)+((*(ptrVolume+indexP))*(*(r+indexP*flag)-*(Rho+indexP*flagrho)))*gravity;
						}
						else if(base->flagTurb==1)
						{
							*(variable+indexP)=*(variable+indexP)+((*(ptrVolume+indexP))*(*(r+indexP*flag)-*(Rho+indexP*flagrho)))*gravity;
							//*(variable+indexP)=*(variable+indexP)- *(Rho+indexP)* gravity * (*(ptrVolume+indexP));
							//printf("\n The body force term is %lf",((*(ptrVolume+indexP))*((*(r+indexP*flag)* *(zone->Flowsol->VarthetaPhase[1]+indexP))-*(Rho+indexP*flagrho)))*gravity);
						}
					}
					if(base->flagTurb==1)//Added by chandrasekhar --- Turbulent MPFV
					{
						index[0]= k ; index[1]=j; index[2]=i;
						*(variable+indexP) = *(variable+indexP) + (*(ptrVolume+indexP))*Findgradient(base,zone,index,variableTag);
						//*(variable+indexP) = *(variable+indexP) + KSource(base,zone,index,variableTag);
						//printf("[Vartag=%d]$[zone %s]$[k=%d$j=%d$i=%d]$Source=%lf\n",variableTag,zone->zonename,k,j,i,*(variable+indexP));
					}
				}
				else if(variableTag==3100||variableTag==3200||variableTag==3300)
				{
					*(variable+indexP) =0.0;
					if(base->BaseParameters->flagGravity==1 && base->flagVariableDensityProblem==0)
					{	
						if(base->flagTurb==0)
						{
							*(variable+indexP)=*(variable+indexP)+((*(ptrVolume+indexP))*(*(r+indexP*flag)-*(Rho+indexP*flagrho)))*gravity;
						}
						else if(base->flagTurb==1)
						{
							*(variable+indexP)=*(variable+indexP)+((*(ptrVolume+indexP))*(*(r+indexP*flag)-*(Rho+indexP*flagrho)))*gravity;
							//*(variable+indexP)=*(variable+indexP)- *(Rho+indexP)* gravity * (*(ptrVolume+indexP));
							
							//printf("\nVartag=%d",variableTag);
							//printf("\n The body force term is %lf",((*(ptrVolume+indexP))*((*(r+indexP*flag)* *(zone->Flowsol->VarthetaPhase[1]+indexP))-*(Rho+indexP*flagrho)))*gravity);
						}
					}
					if(base->flagTurb==1)//Added by chandrasekhar --- Turbulent MPFV
					{
						index[0]= k ; index[1]=j; index[2]=i;
						*(variable+indexP) = *(variable+indexP) + (*(ptrVolume+indexP))*Findgradient(base,zone,index,variableTag);
						//*(variable+indexP) = *(variable+indexP) + KSource(base,zone,index,variableTag);
						//printf("[Vartag=%d]$[zone %s]$[k=%d$j=%d$i=%d]$Source=%lf\n",variableTag,zone->zonename,k,j,i,*(variable+indexP));
					}	
				}
				*(variable_DIAG+indexP) = diagonalCoef;
					
				indexP++;
			}	
		}
		else if(variableTag>=41 && variableTag<=55 )
		{	
			
			indexP=0;
			con=1.0;
			temp =(variableTag-41)/3;
			BETA = zone->Dragcoeff;
			r1 = (zone->Flowsol)->RhoofPhase[0];
			if(temp==1)
				c = 1.0;
			else
				c =-1.0;
			B=zone->Bx;	
			Volfrac = zone->PrevFlowsol->VarthetaPhase[temp];
			r = (zone->Flowsol)->RhoofPhase[temp];
			flag = base->BaseParameters->flagRhoPhase[temp]; 
			temp1 = (variableTag-41)%3;
			if(temp1==0)
			{	
				DiffFlux = zone->Flowsol->diffu;
				drag = zone->DragX;
				Lift = zone->LiftX;
				if(temp==0)
				{	
					
					u = ((zone->PrevFlowsol)->Momentumofphasex[temp+1]);
				}
				else
					u = ((zone->Flowsol)->Momentumofphasex[temp-1]);	
			}
			else if(temp1==1)
			{	
				DiffFlux = zone->Flowsol->diffv;
				drag = zone->DragY;
				Lift = zone->LiftY;
				if(temp==0)
				{	
					
					u = ((zone->PrevFlowsol)->Momentumofphasey[temp+1]);
				}
				else
					u = ((zone->Flowsol)->Momentumofphasey[temp-1]);
			}	
			else if(temp1==2)
			{	
				DiffFlux = zone->Flowsol->diffw;
				drag = zone->DragZ;
				Lift = zone->LiftZ;
				if(temp==0)
				{	
					
					u = ((zone->PrevFlowsol)->Momentumofphasez[temp+1]);
				}
				else
					u = ((zone->Flowsol)->Momentumofphasez[temp-1]);
			}
			

			presgrad=(zone->Flowsol->Pressure_Grad[temp1]);
			
			temp =(variableTag-41)/3;
			
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				*(variable+indexP)=0.0;
	
				*(variable+indexP)=(*(Volfrac+indexP)*(*(presgrad+indexP))-*(DiffFlux+indexP));
				//*(variable+indexP)=(*(Volfrac+indexP)*(*(presgrad+indexP)))-*(DiffFlux+indexP);//Both are same
				if(base->flagTurb==1)//Added by chandrasekhar --- Turbulent MPFV
				{
					index[0]= k ; index[1]=j; index[2]=i;
					*(variable+indexP) = *(variable+indexP) + (*(ptrVolume+indexP))*Findgradient(base,zone,index,variableTag);
					//*(variable+indexP) = *(variable+indexP)+KSource(base,zone,index,variableTag);
					//printf("[Vartag=%d]$[zone %s]$[k=%d$j=%d$i=%d]$Source=%lf\n",variableTag,zone->zonename,k,j,i,*(variable+indexP));
				}
				
				if(base->BaseParameters->flagGravity==1 && base->flagVariableDensityProblem==0)
				{	
					if(temp==0)
					{
 					//not used	
					}
					else if(temp==1)
					{	
						
						r0 = *((zone->Flowsol)->RhoofPhase[1]+indexP*flag);
						if(base->BaseParameters->Dragflag==0)
						{	
							*(variable+indexP)=*(variable+indexP)+(*(ptrVolume+indexP))*(*(r1+indexP*flag)-r0)*gravity*c;
							*(variable_DIAG+indexP)=0.0;
						}	
						else
						{
							if(base->flagTurb==0)
							{
							*(variable+indexP)=*(variable+indexP)+(*(drag+indexP))*(*(ptrVolume+indexP))+(*(ptrVolume+indexP))*(*(B+indexP)*gravity)*c;	
							}
							else if(base->flagTurb==1)
							{
							*(variable+indexP)=*(variable+indexP)+(*(drag+indexP))*(*(ptrVolume+indexP))+(*(Lift+indexP))*(*(ptrVolume+indexP))+(*(ptrVolume+indexP))*(*(B+indexP)*gravity)*c;
								//*(variable+indexP)=*(variable+indexP)+(*(drag+indexP))*(*(ptrVolume+indexP))+(*(Lift+indexP))*(*(ptrVolume+indexP))-(*(ptrVolume+indexP))*gravity*r0 * *(zone->PrevFlowsol->VarthetaPhase[1]+indexP);
							}
							*(variable_DIAG+indexP) = (*(ptrVolume+indexP)*(*(BETA+indexP)));
						}	
 						
					}
					
				}
				
 				
				indexP++;
			}
		}
		else if(variableTag>=34 && variableTag<=38)
		{	
			indexP=0;
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				*(variable+indexP) = 0.0;
				*(variable_DIAG+indexP) = 0.0;
				indexP++;
			}
		}
		else if (variableTag==600)
		{
			indexP=0;domain_prev = MakeDomainArray(zone,variableTag*10);
			if(zone->flagSolidFluid == 0)
			{
				if(base->flagTurbModel==8)
				{
					for(k=0;k<NK-1;k++)
					for(j=0;j<NJ-1;j++)
					for(i=0;i<NI-1;i++)
					{
						index[0]= k ; index[1]=j; index[2]=i;
						*(variable_DIAG+indexP) = SpalartAllmarasSource(base,zone,(variable+indexP),index);
						indexP++;
					}
				}
				else
				{
					for(k=0;k<NK-1;k++)
					for(j=0;j<NJ-1;j++)
					for(i=0;i<NI-1;i++)
					{
						massFlux = zone->Flowsol->Massflux+indexP*6;
						velP = *(domain_prev+indexP);
						sumMassFlux = *(massFlux)+*(massFlux+1)+*(massFlux+2)+*(massFlux+3)+*(massFlux+4)+*(massFlux+5);
						
						*(variable+indexP) = velP*sumMassFlux;
		
						*(variable_DIAG+indexP) = 0.0;
						indexP++;
					}
					
					if(base->flagTurb==1)
					{
						betastar = (zone->Flowsol->Tcoeff->betastar);
						gamma = zone->Flowsol->Tcoeff->gamma;
						beta1  = zone->Flowsol->Tcoeff->beta;
						
						switch(base->flagTurbModel)
						{
							case 2:
							case 6:
								indexP=0;
								for(k=0;k<NK-1;k++)
								for(j=0;j<NJ-1;j++)
								for(i=0;i<NI-1;i++)
								{
									Kold = *(zone->PrevFlowsol->Scal[0]+indexP);
									OmegaOld = *(zone->PrevFlowsol->Scal[1]+indexP);
									vol = *(zone->ZGeom->Volume+indexP);
									
									mue_t=*(zone->Flowsol->Tcoeff->TViscoeff+indexP);
									TProduction = *(zone->Flowsol->Tcoeff->TProduction+indexP);
									
									*(variable+indexP) = *(variable+indexP) + mue_t * TProduction;
									*(variable_DIAG+indexP)    = *(variable_DIAG+indexP) + vol * betastar * OmegaOld;
									indexP++;
								}
							break;
							
							default:
								if(base->flagLogKEps != 1)
								{
									
									Rho = zone->PrevFlowsol->Rho;
									indexP=0;
									for(k=0;k<NK-1;k++)
									for(j=0;j<NJ-1;j++)
									for(i=0;i<NI-1;i++)
									{
										Kold   = *(zone->PrevFlowsol->Scal[0]+indexP);
										Epsold = *(zone->PrevFlowsol->Scal[1]+indexP);
										//Gamma = (Epsold/Kold);
										Gamma = *(zone->Flowsol->Gamma+indexP);//With this also Gamma = (Epsold/Kold); but bounded by Durbin's Realizability constraints.
										//if(Gamma<0)
										//	Gamma = 0.0;
										vol = *(zone->ZGeom->Volume+indexP);
										rho = *(Rho+indexP*flagrho);
										
										mue_t=*(zone->Flowsol->Tcoeff->TViscoeff+indexP);
										S = *(zone->Flowsol->Tcoeff->TProduction+indexP);
										
										TProduction = mue_t*S;
										if(flagLimProduction == 1)
										{
											ProductionMaxLimit = vol*Kold*sqrt((S*rho/vol)/3.0);
											if(TProduction > ProductionMaxLimit)
											{
												TProduction = ProductionMaxLimit;
												base->BaseParameters->Turb->TProductionMaxLimitCount++;
											}
										}
										//TProduction = *(zone->Flowsol->Tcoeff->TProduction+indexP);
										
										//NOTE:WITHOUT Newton's linearization technique for dissipation term
										//Use the following statement if TProduction = mue_t*S; (presently active)
// 										*(variable+indexP) = *(variable+indexP) + TProduction - vol * Epsold;
										//Use the following statement if TProduction = S;
// 										*(variable+indexP) = *(variable+indexP) + mue_t * TProduction - vol * Epsold;
										
										//NOTE:WITH Newton's linearization technique for dissipation term
										//Use the following statement if TProduction = S;
// 										*(variable+indexP) = *(variable+indexP) + mue_t * TProduction;
										//Use the following statement if TProduction = mue_t*S; (presently active)
										*(variable+indexP) = *(variable+indexP) + TProduction;
										*(variable_DIAG+indexP)    = *(variable_DIAG+indexP) + vol * (Epsold/Kold);
// 										*(variable_DIAG+indexP)    = *(variable_DIAG+indexP) + vol * Gamma;
										
										if(zone->Flowsol->flagTurb == 4)
											*(variable+indexP)= *(variable+indexP) - *(zone->Flowsol->Tcoeff->PhiK+indexP);
										//if(indexP==0)
										//	printf("Source = %e\n",*(variable+indexP));
										
										if(base->flagModulation==1)
										{	
											if((base->flagParticle==1)|| (base->flagDroplet==1))
											{
											t_p=*(zone->t_p+indexP);
											t_e=0.125*(Kold/Epsold);
											t_L=Kold/Epsold;
											if((t_p==0.0)||(*(zone->PrevFlowsol->VarthetaDisperse+indexP)<EPSILON))
												K_Source=0.0;
											else
											{
											K_Source=-2.0*(*(zone->PrevFlowsol->VarthetaDisperse+indexP))*(*(zone->Flowsol->RhoDisperse))/t_p;
											
											if(base->flagTurbModel==1)
												K_Source=K_Source*(1.0-exp(-0.0825*t_p/t_e));
											else if(base->flagTurbModel==3)
												K_Source=K_Source*(1.0-exp(-0.09*t_p/t_L));
										
											K_Source=K_Source/(*(zone->PrevFlowsol->Rho+indexP*flagrho));
											*(variable+indexP)= *(variable+indexP)+(*(ptrVolume+indexP))*K_Source*Kold;
											//printf("\n K_Source*Kold=%e",K_Source*Kold);
											//*(variable_DIAG+indexP)=*(variable_DIAG+indexP)-(*(ptrVolume+indexP))*K_Source;
											//printf("\n K Source=%e",K_Source*(1-(*(zone->PrevFlowsol->VarthetaDisperse+indexP)))*Kold);
											}
											}
										}
										indexP++;
									}
								}
								else
								{	
									indexP=0;
									for(k=0;k<NK-1;k++)
									for(j=0;j<NJ-1;j++)
									for(i=0;i<NI-1;i++)
									{
										Kold   = *(zone->PrevFlowsol->Scal[0]+indexP);
										Epsold = *(zone->PrevFlowsol->Scal[1]+indexP);
										vol = *(zone->ZGeom->Volume+indexP);
										
										sumGradSquare = CalcGradSquare(base,zone,index,variableTag);
										MuEffK = *(zone->Flowsol->Scalcoeff[0]+indexP);
										mue_t=*(zone->Flowsol->Tcoeff->TViscoeff+indexP);
										TProduction = *(zone->Flowsol->Tcoeff->TProduction+indexP);
							
										*(variable+indexP) 	= *(variable+indexP) + vol * (MuEffK/rho) * sumGradSquare;
										*(variable+indexP) 	= *(variable+indexP) + mue_t * TProduction * exp(-Kold) * (1.0+Kold);
										Dissipation 		= vol * exp(Epsold - Kold) * (1.0+Kold);
										*(variable_DIAG+indexP) = - (vol * exp(Epsold - Kold)) + mue_t * TProduction * exp(-Kold);
										if(base->BaseParameters->TimeStepingSchemeTag==1)
											*(variable_DIAG+indexP) = *(variable_DIAG+indexP)/2.0;
							
										*(variable+indexP) = *(variable+indexP) - Dissipation;
										indexP++;
									}
								}
							break;
						}
					}
				}
			}
			else
			{
				for(k=0;k<NK-1;k++)
				for(j=0;j<NJ-1;j++)
				for(i=0;i<NI-1;i++)
				{
					*(variable+indexP)=0.0;
					*(variable_DIAG+indexP) = 0.0;
					indexP++;
				}
			}
		}
		else if (variableTag==601)
		{
			indexP=0;domain_prev = MakeDomainArray(zone,variableTag*10);
			if(zone->flagSolidFluid == 0)
			{
				for(k=0;k<NK-1;k++)
				for(j=0;j<NJ-1;j++)
				for(i=0;i<NI-1;i++)
				{
					massFlux = zone->Flowsol->Massflux+indexP*6;
					velP = *(domain_prev+indexP);
					sumMassFlux = *(massFlux)+*(massFlux+1)+*(massFlux+2)+*(massFlux+3)+*(massFlux+4)+*(massFlux+5);
				
					*(variable+indexP) = velP*sumMassFlux;
	
					*(variable_DIAG+indexP) = 0.0;
					indexP++;
				}
				if(base->flagTurb==1)
				{
					betastar = (zone->Flowsol->Tcoeff->betastar);
					gamma = zone->Flowsol->Tcoeff->gamma;
					beta1  = zone->Flowsol->Tcoeff->beta;
					
					switch(base->flagTurbModel)
					{
						case 2:
							indexP=0;
							for(k=0;k<NK-1;k++)
							for(j=0;j<NJ-1;j++)
							for(i=0;i<NI-1;i++)
							{
								TProduction = *(zone->Flowsol->Tcoeff->TProduction+indexP);
								rho = *(zone->PrevFlowsol->Rho+indexP*flagrho);
								*(variable+indexP) = *(variable+indexP) + gamma *rho *TProduction;
								
								OmegaOld = *(zone->PrevFlowsol->Scal[1]+indexP);
								vol = *(zone->ZGeom->Volume+indexP);
								
								*(variable_DIAG+indexP)    = *(variable_DIAG+indexP)+2.0*vol * beta1 * OmegaOld;
								*(variable+indexP) 	= *(variable+indexP) + vol * beta1 * OmegaOld *OmegaOld;
								indexP++;
							}
						break;
						case 6:
							indexP=0;
							Gamma1=zone->Flowsol->Tcoeff->Gamma1;
							Beta1=zone->Flowsol->Tcoeff->Beta1;
							Gamma2=zone->Flowsol->Tcoeff->Gamma2;
							Beta2=zone->Flowsol->Tcoeff->Beta2;
							for(k=0;k<NK-1;k++)
							for(j=0;j<NJ-1;j++)
							for(i=0;i<NI-1;i++)
							{
								OmegaOld = *(zone->PrevFlowsol->Scal[1]+indexP);
								vol = *(zone->ZGeom->Volume+indexP);
								F1= *(zone->Flowsol->Tcoeff->F1+indexP);
								gamma = F1 * Gamma1 + (1.0-F1) * Gamma2;
								beta1 = F1 * Beta1  + (1.0-F1) * Beta2;
								CrossDiffusion = *(zone->Flowsol->CrossDiff+indexP);
								TProduction = *(zone->Flowsol->Tcoeff->TProduction+indexP);
								rho = *(zone->PrevFlowsol->Rho+indexP*flagrho);
						
								*(variable+indexP) = *(variable+indexP) + gamma * rho * TProduction;
								*(variable+indexP) = *(variable+indexP) + CrossDiffusion;
						
								*(variable_DIAG+indexP)    = *(variable_DIAG+indexP)+2.0* vol * beta1 * OmegaOld;
								*(variable+indexP)  	 = *(variable+indexP) + vol * beta1 * OmegaOld *OmegaOld;
								
								indexP++;
							}
							break;
						
						default:
							if(base->flagLogKEps != 1)
							{
								Rho = zone->PrevFlowsol->Rho;
								indexP=0;
								for(k=0;k<NK-1;k++)
								for(j=0;j<NJ-1;j++)
								for(i=0;i<NI-1;i++)
								{
									C0epsi = *(zone->Flowsol->Tcoeff->C0epsi+indexP);//printf("\n C0epsi=%lf",C0epsi);
									C1epsi = *(zone->Flowsol->Tcoeff->C1epsi+indexP);//printf("\n C1epsi=%lf",C1epsi);
									C2epsi = *(zone->Flowsol->Tcoeff->C2epsi+indexP);//printf("\n C2epsi=%lf",C2epsi);
									
									Kold   = *(zone->PrevFlowsol->Scal[0]+indexP);
									Epsold = *(zone->PrevFlowsol->Scal[1]+indexP);
									if((base->flagParticle==1))
									{
									Gamma = (Epsold/Kold);
									}
									else
									{
									Gamma = *(zone->Flowsol->Gamma+indexP);
									}
									//if(Gamma<0)
									//	Gamma = 0.0;
									
									vol = *(zone->ZGeom->Volume+indexP);
									rho = *(Rho+indexP*flagrho);
									
									mue_t=*(zone->Flowsol->Tcoeff->TViscoeff+indexP);
									S = *(zone->Flowsol->Tcoeff->TProduction+indexP);
									
									TProduction = mue_t*S;
									
									if(flagLimProduction == 1)
									{
										ProductionMaxLimit = vol*Kold*sqrt((S*rho/vol)/3.0);
										if(TProduction > ProductionMaxLimit)
										{
											TProduction = ProductionMaxLimit;
											base->BaseParameters->Turb->TProductionMaxLimitCount++;
										}
									}
									//TProduction = *(zone->Flowsol->Tcoeff->TProduction+indexP);
									
									*(variable+indexP) = *(variable+indexP) + C1epsi * Gamma*TProduction;
									//*(variable+indexP) = *(variable+indexP) + mue_t *C1epsi * Gamma*TProduction;
									
									*(variable_DIAG+indexP)    = *(variable_DIAG+indexP)+vol * C2epsi*Gamma *2.0;
									
									
									Dissipation = -1.0* vol * C2epsi*Gamma * Epsold;
									
									if((zone->Flowsol->flagRNG_KEPS == 1) && (zone->Flowsol->flagTurb == 1))
									{
										*(variable_DIAG+indexP)    = *(variable_DIAG+indexP)    + vol * C0epsi * Gamma;///Ask Narsimha about C0epsi initiliazation.
									}
									*(variable+indexP) = *(variable+indexP) - Dissipation;
									if(zone->Flowsol->flagTurb == 4)
										*(variable+indexP)= *(variable+indexP) - *(zone->Flowsol->Tcoeff->PhiK+indexP);
									if(base->flagModulation==1)
									{
										if(base->flagTurbModel==1)
										{
											if((base->flagParticle==1)|| (base->flagDroplet==1))
											{
												t_p=*(zone->t_p+indexP);
												if((t_p==0.0)||(*(zone->PrevFlowsol->VarthetaDisperse+indexP)<EPSILON))
												K_Source=0.0;
												else
												{
												K_Source=-2.0*(*(zone->PrevFlowsol->VarthetaDisperse+indexP))*(*(zone->Flowsol->RhoDisperse))/t_p;
												K_Source=K_Source/(*(zone->PrevFlowsol->Rho+indexP*flagrho));
												*(variable+indexP)    = *(variable+indexP)+(*(ptrVolume+indexP))*K_Source*Epsold;
												//printf("\n K_Source*Kold=%e",K_Source*Epsold);
												//*(variable_DIAG+indexP)=*(variable_DIAG+indexP)-(*(ptrVolume+indexP))*K_Source;
												//printf("\n E source=%e",K_Source*(1-(*(zone->PrevFlowsol->VarthetaDisperse+indexP)))*Epsold);
												}
											}
										}
										if(base->flagTurbModel==3)
										{
											if((base->flagParticle==1)|| (base->flagDroplet==1))
											{
												t_p=*(zone->t_p+indexP);
												t_e=0.125*(Kold/Epsold);
												//t_L=Kold/Epsold;
												if((t_p==0.0)||(*(zone->PrevFlowsol->VarthetaDisperse+indexP)<EPSILON))
												K_Source=0.0;
												else
												{
												K_Source=-2.0*(*(zone->PrevFlowsol->VarthetaDisperse+indexP))*(*(zone->Flowsol->RhoDisperse))/t_p;
												//K_Source=K_Source*(1.0-exp(-0.0825*t_p/t_e));
												K_Source=K_Source*(1.0-exp(-0.4*t_p/t_e));
												K_Source=K_Source/(*(zone->PrevFlowsol->Rho+indexP*flagrho));
												*(variable+indexP)= *(variable+indexP)+(*(ptrVolume+indexP))*K_Source*Epsold;
												//printf("\n K_Source*Kold=%e",K_Source*Epsold);
												//*(variable_DIAG+indexP)=*(variable_DIAG+indexP)-(*(ptrVolume+indexP))*K_Source;
												//printf("\n K Source=%e",K_Source*(1-(*(zone->PrevFlowsol->VarthetaDisperse+indexP)))*Kold);
												}
											}
										}
									}
									

									indexP++;
								}
							}
							else
							{
								indexP=0;
								for(k=0;k<NK-1;k++)
								for(j=0;j<NJ-1;j++)
								for(i=0;i<NI-1;i++)
								{
									C0epsi = *(zone->Flowsol->Tcoeff->C0epsi+indexP);
									C1epsi = *(zone->Flowsol->Tcoeff->C1epsi+indexP);
									C2epsi = *(zone->Flowsol->Tcoeff->C2epsi+indexP);
									
									Kold   = *(zone->PrevFlowsol->Scal[0]+indexP);
									Epsold = *(zone->PrevFlowsol->Scal[1]+indexP);
									K = *(zone->Flowsol->Scal[0]+indexP);
									vol = *(zone->ZGeom->Volume+indexP);
									
									
									sumGradSquare = CalcGradSquare(base,zone,index,variableTag);
									MuEffEps = *(zone->Flowsol->Scalcoeff[1]+indexP);
									mue_t=*(zone->Flowsol->Tcoeff->TViscoeff+indexP);
									TProduction = *(zone->Flowsol->Tcoeff->TProduction+indexP);
									rho = *(zone->PrevFlowsol->Rho+indexP*flagrho);
									
									*(variable+indexP)= *(variable+indexP) + vol * (MuEffEps/rho)* sumGradSquare;
									*(variable+indexP) = *(variable+indexP) + C1epsi * mue_t * TProduction * exp(-K);
									Dissipation = vol * C2epsi * exp(Epsold-K) * (1.0-Epsold);
									*(variable_DIAG+indexP) = vol * C2epsi * exp(Epsold-K);
									if(base->BaseParameters->TimeStepingSchemeTag==1)
										*(variable_DIAG+indexP) = *(variable_DIAG+indexP)/2.0;
									*(variable+indexP) = *(variable+indexP) - Dissipation;
									//Fsum1 = Fsum1 + diacoeff * Eps;
									indexP++;
								}
							}
							break;
					}
				}
			}
			else
			{
				for(k=0;k<NK-1;k++)
				for(j=0;j<NJ-1;j++)
				for(i=0;i<NI-1;i++)
				{
					*(variable+indexP)=0.0;
					*(variable_DIAG+indexP) = 0.0;
					indexP++;
				}
			}
		}
		else if(variableTag==602)
		{
			indexP=0;
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				massFlux = zone->Flowsol->Massflux+indexP*6;
				velP = *(domain_prev+indexP);
				sumMassFlux = *(massFlux)+*(massFlux+1)+*(massFlux+2)+*(massFlux+3)+*(massFlux+4)+*(massFlux+5);
			
				*(variable+indexP) = velP*sumMassFlux;

				*(variable_DIAG+indexP) = 0.0;
				indexP++;
			}

		}
		else if(variableTag==603)
		{
			indexP=0;
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				massFlux = zone->Flowsol->Massflux+indexP*6;
				velP = *(domain_prev+indexP);
				sumMassFlux = *(massFlux)+*(massFlux+1)+*(massFlux+2)+*(massFlux+3)+*(massFlux+4)+*(massFlux+5);
			
				*(variable+indexP) = velP*sumMassFlux;

				*(variable_DIAG+indexP) = 0.0;
				indexP++;
			}

		}
		else if(variableTag==604)
		{
			indexP=0;
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				massFlux = zone->Flowsol->Massflux+indexP*6;
				velP = *(domain_prev+indexP);
				sumMassFlux = *(massFlux)+*(massFlux+1)+*(massFlux+2)+*(massFlux+3)+*(massFlux+4)+*(massFlux+5);
			
				*(variable+indexP) = velP*sumMassFlux;

				*(variable_DIAG+indexP) = 0.0;
				indexP++;
			}

		}
		else if(variableTag==605)
		{
			indexP=0;
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				massFlux = zone->Flowsol->Massflux+indexP*6;
				velP = *(domain_prev+indexP);
				sumMassFlux = *(massFlux)+*(massFlux+1)+*(massFlux+2)+*(massFlux+3)+*(massFlux+4)+*(massFlux+5);
			
				*(variable+indexP) = velP*sumMassFlux;

				*(variable_DIAG+indexP) = 0.0;
				indexP++;
			}

		}
		else if(variableTag==606)
		{
			indexP=0;
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				massFlux = zone->Flowsol->Massflux+indexP*6;
				velP = *(domain_prev+indexP);
				sumMassFlux = *(massFlux)+*(massFlux+1)+*(massFlux+2)+*(massFlux+3)+*(massFlux+4)+*(massFlux+5);
			
				*(variable+indexP) = velP*sumMassFlux;

				*(variable_DIAG+indexP) = 0.0;
				indexP++;
			}

		}
		else if(variableTag==607)
		{
			indexP=0;
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				massFlux = zone->Flowsol->Massflux+indexP*6;
				velP = *(domain_prev+indexP);
				sumMassFlux = *(massFlux)+*(massFlux+1)+*(massFlux+2)+*(massFlux+3)+*(massFlux+4)+*(massFlux+5);
			
				*(variable+indexP) = velP*sumMassFlux;

				*(variable_DIAG+indexP) = 0.0;
				indexP++;
			}

		}
		else if(variableTag==608)
		{
			indexP=0;
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				massFlux = zone->Flowsol->Massflux+indexP*6;
				velP = *(domain_prev+indexP);
				sumMassFlux = *(massFlux)+*(massFlux+1)+*(massFlux+2)+*(massFlux+3)+*(massFlux+4)+*(massFlux+5);
			
				*(variable+indexP) = velP*sumMassFlux;

				*(variable_DIAG+indexP) = 0.0;
				indexP++;
			}

		}
		else if(variableTag==609)
		{
			indexP=0;
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				massFlux = zone->Flowsol->Massflux+indexP*6;
				velP = *(domain_prev+indexP);
				sumMassFlux = *(massFlux)+*(massFlux+1)+*(massFlux+2)+*(massFlux+3)+*(massFlux+4)+*(massFlux+5);
			
				*(variable+indexP) = velP*sumMassFlux;

				*(variable_DIAG+indexP) = 0.0;
				indexP++;
			}

		}
		else if(variableTag== 62 || variableTag== 63 || variableTag== 64 || variableTag== 6200 || variableTag== 6300 || variableTag== 6400 )
		{
			
			indexP=0;

			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				*(variable+indexP)=0.0;//printf("\n It is %d",variableTag);
				*(variable_DIAG+indexP)=0.0;
				
				//Added by Mrunalini.. Gas Particle n Gas Droplet flows with evaporation
				massFlux = zone->Flowsol->MassfluxDisperse+indexP*6;
				velP = *(domain_prev+indexP);
				sumMassFlux = *(massFlux)+*(massFlux+1)+*(massFlux+2)+*(massFlux+3)+*(massFlux+4)+*(massFlux+5);
				*(variable+indexP) = *(variable+indexP) + velP*sumMassFlux;

				indexP++;
			}
			
			if(variableTag==62 || variableTag==63 || variableTag==64)
			{
				if(base->flagDroplet==1 && base->BaseParameters->flagVariableDia==1)
				{
					indexP=0;

					for(k=0;k<NK-1;k++)
					for(j=0;j<NJ-1;j++)
					for(i=0;i<NI-1;i++)
					{
						*(variable+indexP)=*(variable+indexP) + (1/(*(zone->Flowsol->RhoDisperse+indexP*flagrhod)))*(*(zone->Flowsol->Pressure_Grad[variableTag-62]+indexP));	
						indexP++;
					}
				}
			}
			
			/*if((base->flagParticle==1)|| (base->flagDroplet==1))//do we need this if at all? adhiraj
			{*/
				indexP=0;
				for(k=0;k<NK-1;k++)
				for(j=0;j<NJ-1;j++)
				for(i=0;i<NI-1;i++)
				{
					if((*(zone->PrevFlowsol->VarthetaDisperse+indexP))<=EPSILON)
		   			{
						*(variable+indexP)=*(variable+indexP) + 0.0;
					}
		  			else
					{	//changed by Mrunalini
						*(variable+indexP)=*(variable+indexP)
							+ (*(ptrVolume+indexP))*(*(Drag + indexP))/(*(zone->Flowsol->RhoDisperse+indexP*flagrhod)); 
					}
					indexP++;
				}
			//}
			
	//	if(base->flagTurb==1&&base->flagParticle==1)//Added by chandrasekhar --- Turbulent Gas-Particle
		//	{
		//		indexP=0;
	//			for(k=0;k<NK-1;k++)
	//			for(j=0;j<NJ-1;j++)
	//			for(i=0;i<NI-1;i++)
//				{
//					index[0]= k ; index[1]=j; index[2]=i;
//					*(variable+indexP) = *(variable+indexP) + (*(ptrVolume+indexP))*Findgradient(base,zone,index,variableTag);
//					indexP++;
	//			}
//			}//
			//if(base->flagTurb==1&&base->flagParticle==1)//Added by chandrasekhar --- Turbulent Gas-Particle
		//	{
	//			indexP=0;
//				for(k=0;k<NK-1;k++)
	//			for(j=0;j<NJ-1;j++)
//				for(i=0;i<NI-1;i++)
//				{
	//				index[0]= k ; index[1]=j; index[2]=i;
	//				*(variable+indexP) = *(variable+indexP) + (*(ptrVolume+indexP))*Findgradient_Disperse(base,zone,index,variableTag);
	//				indexP++;
//				}
//			}
			if(base->flagTurb==1&&base->flagParticle==1)
			{ printf("\n xxxxxxxxxxxxxxxxxxx\n ");
				DiffFlux=zone->Flowsol->Diff_Disperse;
				flagrhod=base->BaseParameters->flagRhoDisperse;
				if(variableTag==62)
				{
					indexP=0;
					for(k=0;k<NK-1;k++)
					for(j=0;j<NJ-1;j++)
					for(i=0;i<NI-1;i++)
					{
						if((*(zone->PrevFlowsol->VarthetaDisperse+indexP))>EPSILON)
							Ud= *(zone->PrevFlowsol->MomentumDispersex+indexP)/(*(zone->Flowsol->RhoDisperse+indexP*flagrhod)* *(zone->PrevFlowsol->VarthetaDisperse+indexP));
						else
							Ud=*(zone->PrevFlowsol->u+indexP);
						
						*(variable+indexP) = *(variable+indexP)- *(DiffFlux+indexP)* *(zone->Flowsol->RhoDisperse)* *(zone->Flowsol->Tcoeff->DiffVolfrac_disperse+indexP)* Ud;
						indexP++;
					}
				}		
				if(variableTag==63)
				{
					indexP=0;
					for(k=0;k<NK-1;k++)
					for(j=0;j<NJ-1;j++)
					for(i=0;i<NI-1;i++)
					{
						if((*(zone->PrevFlowsol->VarthetaDisperse+indexP))>EPSILON)
							Vd= *(zone->PrevFlowsol->MomentumDispersey+indexP)/(*(zone->Flowsol->RhoDisperse+indexP*flagrhod)* *(zone->PrevFlowsol->VarthetaDisperse+indexP));
						else
							Vd=*(zone->PrevFlowsol->v+indexP);
						
						*(variable+indexP) = *(variable+indexP)- *(DiffFlux+indexP)* *(zone->Flowsol->RhoDisperse)* *(zone->Flowsol->Tcoeff->DiffVolfrac_disperse+indexP)* Vd;
						indexP++;
					}
						
				}
				if(variableTag==64)
				{
					indexP=0;
					for(k=0;k<NK-1;k++)
					for(j=0;j<NJ-1;j++)
					for(i=0;i<NI-1;i++)
					{
						if((*(zone->PrevFlowsol->VarthetaDisperse+indexP))>EPSILON)
							Wd= *(zone->PrevFlowsol->MomentumDispersez+indexP)/(*(zone->Flowsol->RhoDisperse+indexP*flagrhod)* *(zone->PrevFlowsol->VarthetaDisperse+indexP));
						else
							Wd=*(zone->PrevFlowsol->w+indexP);
						
						*(variable+indexP) = *(variable+indexP)- *(DiffFlux+indexP)* *(zone->Flowsol->RhoDisperse)* *(zone->Flowsol->Tcoeff->DiffVolfrac_disperse+indexP)* Wd;
						indexP++;
					}
						
				}
			}
			if(base->BaseParameters->flagGravity==1)
			{
				indexP=0;
				for(k=0;k<NK-1;k++)
				for(j=0;j<NJ-1;j++)
				for(i=0;i<NI-1;i++)
				{ 
                    			if((*(zone->PrevFlowsol->VarthetaDisperse+indexP))<=EPSILON)
		   			{
						*(variable+indexP)=*(variable+indexP) + 0.0;
					}
		    			else
		        		{ //changed by Mrunalini
						*(variable+indexP)=*(variable+indexP)+(*(ptrVolume+indexP))*gravity; 
		        		}
					indexP++;
				}
				
			}
		}
		else if(variableTag==65)
		{
			indexP=0;
			
			if((base->BaseParameters->flagVariableDia==1)&&(base->flagDroplet==1))
			{//changed by Mrunalini
				indexP=0;
				for(k=0;k<NK-1;k++)
				for(j=0;j<NJ-1;j++)
				for(i=0;i<NI-1;i++)
				{
					*(variable+indexP)=0.0;
					*(variable_DIAG+indexP)=0.0;
									
					if((*(zone->PrevFlowsol->VarthetaDisperse+indexP))<=EPSILON)
		   			{
						*(variable+indexP)=*(variable+indexP) + 0.0;
					}
		  		  	else
				   	{
						VolumeFrac = *(zone->PrevFlowsol->VarthetaDisperse+indexP)*base->BaseParameters->InitialVolumeFraction;
						*(variable_DIAG+indexP)=*(variable_DIAG+indexP)+ (*(ptrVolume+indexP)/(*(zone->Flowsol->RhoDisperse+indexP*flagrhod)*VolumeFrac))*(*(zone->PrevFlowsol->Evaporate+indexP));
					}
					indexP++;
				}	
			}
			else
			{
				for(k=0;k<NK-1;k++)
				for(j=0;j<NJ-1;j++)
				for(i=0;i<NI-1;i++)
				{
					*(variable+indexP)=0.0;
					*(variable_DIAG+indexP)=0.0;
					indexP++;
				}
			}
		}

		else if(variableTag==710)
		{
			if((base->BaseParameters->flagVariableDia==1)&&(base->flagDroplet==1))
			{//changed by Mrunalini
				indexP=0;

				for(k=0;k<NK-1;k++)
				for(j=0;j<NJ-1;j++)
				for(i=0;i<NI-1;i++)
				{
					*(variable+indexP)=0.0;
					*(variable_DIAG+indexP)=0.0;

					if((*(zone->PrevFlowsol->VarthetaDisperse+indexP))<=EPSILON)
		   			{
						*(variable+indexP)=*(variable+indexP) + 0.0;
					}
		  			else
					{
						VolumeFrac = *(zone->Flowsol->VarthetaDisperse+indexP)*base->BaseParameters->InitialVolumeFraction;
						*(variable+indexP)=*(variable+indexP) + (*(ptrVolume+indexP))*(*(zone->PrevFlowsol->Evaporate+indexP))/(1-VolumeFrac);
					}
					indexP++;
				}
				
			}
			else
			{
				indexP=0;

				for(k=0;k<NK-1;k++)
				for(j=0;j<NJ-1;j++)
				for(i=0;i<NI-1;i++)
				{
					*(variable+indexP)=0.0;
					*(variable_DIAG+indexP)=0.0;
					indexP++;
				}
			}
		}
		else if(variableTag==68)
		{
			indexP=0;
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				*(variable+indexP)=0.0;
				*(variable_DIAG+indexP)=0.0;
		             
				massFlux = zone->Flowsol->Massflux+indexP*6;
				velP = *(domain_prev+indexP);
								
				sumMassFlux = *(massFlux)+*(massFlux+1)+*(massFlux+2)+*(massFlux+3)+*(massFlux+4)+*(massFlux+5);
				
				*(variable+indexP) = *(variable+indexP) + velP*sumMassFlux;

				if((*(zone->PrevFlowsol->VarthetaDisperse+indexP))<=EPSILON)
		   		{
					*(variable+indexP)=*(variable+indexP) + 0.0;
				}
		  		else
				{	   
					VolumeFrac = *(zone->Flowsol->VarthetaDisperse+indexP)*base->BaseParameters->InitialVolumeFraction;
					*(variable_DIAG+indexP) = *(variable_DIAG+indexP)+ (1.0/(*(zone->Flowsol->Rho+indexP*flagrho)*(1-VolumeFrac)))*(*(ptrVolume+indexP))*(*(zone->PrevFlowsol->Evaporate+indexP)); //*(zone->Flowsol->Rho+indexP*flagrho)
					*(variable+indexP)=*(variable+indexP)+ ((*(ptrVolume+indexP))/(*(zone->Flowsol->Rho+indexP*flagrho)*(1-VolumeFrac)))*(*(zone->PrevFlowsol->Evaporate+indexP)); //
				}
		 		
				indexP++;
		        }
		}
		else if(variableTag==69)
		{
			if((base->BaseParameters->flagVariableDia==1)&&(base->flagDroplet==1))
		  	{
				indexP=0;

				for(k=0;k<NK-1;k++)
				for(j=0;j<NJ-1;j++)
				for(i=0;i<NI-1;i++)
				{
					*(variable+indexP)=0.0;
					*(variable_DIAG+indexP)=0.0;

					massFlux = zone->Flowsol->MassfluxDisperse+indexP*6;
					velP = *(domain_prev+indexP)*(*(zone->PrevFlowsol->CPdrop+indexP*flagCPDrop));
					sumMassFlux = *(massFlux)+*(massFlux+1)+*(massFlux+2)+*(massFlux+3)+*(massFlux+4)+*(massFlux+5);
				
					*(variable+indexP) = *(variable+indexP) + velP*sumMassFlux;
				
         	  		  	if((*(zone->PrevFlowsol->VarthetaDisperse+indexP))<=EPSILON)
		   			{
						*(variable+indexP)=*(variable+indexP) + 0.0;
					}
		  		  	else
				   	{	
						VolumeFrac = *(zone->Flowsol->VarthetaDisperse+indexP)*base->BaseParameters->InitialVolumeFraction;
										
				   		*(variable+indexP) = *(variable+indexP) + (*(ptrVolume+indexP)/(*(zone->PrevFlowsol->RhoDisperse)*VolumeFrac))*((*(zone->PrevFlowsol->QDroptoGas+indexP))-(*(zone->PrevFlowsol->QLatent+indexP)));
				   	}
					 indexP++;	
				  }	
			}
			else
			{
				indexP =0;
				for(k=0;k<NK-1;k++)
				for(j=0;j<NJ-1;j++)
				for(i=0;i<NI-1;i++)
				{
					*(variable+indexP)=0.0;
					*(variable_DIAG+indexP)=0.0;
			
					indexP++;
				}
			}			
		}
//End of source terms for disperse phase adhi+shash #############################
		else
		{
			printf("Your module yet to be implemented in CalcSource_opti.c\n\n");
			printf("Remember that this is a newly added function and your variables has to be incorporated here \n");
			exit(123);
		}
	}

// 	if(base->BaseParameters->Turb->TProductionMaxLimitCount != 0)//This should be done only for k-epsilon model
// 		printf("######### Max. Turbulent Production is limited in %d Cell(s) ########\n", base->BaseParameters->Turb->TProductionMaxLimitCount);

	return (diagonalCoef);
}

