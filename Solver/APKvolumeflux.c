#include "APKvolumeflux.h"
/**
This function is used only once for calculating the volumeflux just before the first timestep using
the guess values of mass velocities. This is used for 3D geometries.
*/
int volumeflux(Base *base)
{
	int nZones,numZones;
	int k,j,i,index,indexB,indexT,indexE,indexW,indexS,indexN,indexP;
	int numiface,numBc,face,iface,zBc,dir,Index;
	double *SurfArea,*F,*FD;
	double *domainVolume,*domainUmass,*domainVmass,*domainWmass,*domainRho;
	double *slabVolume[6],*slabUmass[6],*slabVmass[6],*slabWmass[6],*slabRho[6];
	
	double *domainRhoDisperse,*domainVartheta,*domainUDisperse,*domainVDisperse,*domainWDisperse;
	double *slabVartheta[6],*slabUDisperse[6],*slabVDisperse[6],*slabWDisperse[6],*slabRhoDisperse[6];
	Zone* zone;
	
	int NI,NJ,NK,*ptrSize;
	double uw,vw,ww,ue,ve,we,un,vn,wn,us,vs,ws,ut,vt,wt,ub,vb,wb;
	double u_P,u_E,u_W,u_S,u_N,u_T,u_B,v_P,v_E,v_W,v_S,v_N,v_T,v_B,w_P,w_E,w_W,w_S,w_N,w_T,w_B;
	//added by adhiraj for disperse phase
	double uwD,vwD,wwD,ueD,veD,weD,unD,vnD,wnD,usD,vsD,wsD,utD,vtD,wtD,ubD,vbD,wbD;
	double u_PD,u_ED,u_WD,u_SD,u_ND,u_TD,u_BD,v_PD,v_ED,v_WD,v_SD,v_ND,v_TD,v_BD,w_PD,w_ED,w_WD,w_SD,w_ND,w_TD,w_BD;
	double vol_P,vol_E,vol_W,vol_N,vol_S,vol_T,vol_B,rho_P,rho_E,rho_W,rho_S,rho_N,rho_T,rho_B;
	double rho_PD,rho_ED,rho_WD,rho_SD,rho_ND,rho_TD,rho_BD;

	double *slabUDmass[6], *slabVDmass[6], *slabWDmass[6], *domainUDmass, *domainVDmass, *domainWDmass;
	
	double EPSILON;
	EPSILON = base->BaseParameters->EPSILONDrop;
	int che;	
	
// 	FILE *fp1;
// 	fp1 = fopen("volflux.dat","w");
	printf("\nEntering volume flux\n");

	nZones=base->Number_of_zones;	
	for(numZones=0;numZones<nZones;numZones++)
	{
		/*Copying down u,v,w into umass,vmass,wmass ....domain as well as slabs*/
		ptrSize=((base)->CgnsZone+numZones)->Size;
		zone=(base->CgnsZone+numZones);
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		
		for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{
			indexP=(NJ-1)*(NI-1)*k+(NI-1)*j+i;
			if(base->flagMPFV==0)
			{
				*((zone->Flowsol)->umass+indexP)=*((zone->Flowsol)->u+indexP);
				*((zone->Flowsol)->vmass+indexP)=*((zone->Flowsol)->v+indexP);
				*((zone->Flowsol)->wmass+indexP)=*((zone->Flowsol)->w+indexP);

				if(base->flagDroplet==1 && base->BaseParameters->flagVariableDia==1)
				{
					*((zone->Flowsol)->U_Dmass+indexP)=*((zone->Flowsol)->U_Disperse+indexP);
					*((zone->Flowsol)->V_Dmass+indexP)=*((zone->Flowsol)->V_Disperse+indexP);
					*((zone->Flowsol)->W_Dmass+indexP)=*((zone->Flowsol)->W_Disperse+indexP);
				}

			}
			else if(base->flagMPFV==1)
			{
				*((zone->Flowsol)->Momentumxstar+indexP)=*((zone->Flowsol)->Momentumx+indexP);
				*((zone->Flowsol)->Momentumystar+indexP)=*((zone->Flowsol)->Momentumy+indexP);
				*((zone->Flowsol)->Momentumzstar+indexP)=*((zone->Flowsol)->Momentumz+indexP);
			}
		}
		numiface=(zone)->Number_of_Interfaces;
		
		numBc=(zone)->Number_of_BC_Patches;
		
		for(dir=1;dir<=6;dir++)
		{
			for(iface=0;iface<numiface;iface++)
			{
				if(((zone)->Iface+iface)->Direction==dir)
				{
					for(k=0;k<NK-1;k++)
					for(j=0;j<NJ-1;j++)
					for(i=0;i<NI-1;i++)
					{
						if(dir==1||dir==2)
							index= k*(NJ-1)+j;
						if(dir==3||dir==4)
							index= k*(NI-1)+i;
						if(dir==5||dir==6)
							index= j*(NI-1)+i;
						if(base->flagMPFV==0)
						{
							*((((zone)->Iface+iface)->ZSlab)->umass+index) = *((((zone)->Iface+iface)->ZSlab)->u+index);
							*((((zone)->Iface+iface)->ZSlab)->vmass+index) = *((((zone)->Iface+iface)->ZSlab)->v+index);
							*((((zone)->Iface+iface)->ZSlab)->wmass+index) = *((((zone)->Iface+iface)->ZSlab)->w+index);

							if(base->flagDroplet==1 && base->BaseParameters->flagVariableDia==1)
							{
								*((((zone)->Iface+iface)->ZSlab)->U_Dmass+index) = *((((zone)->Iface+iface)->ZSlab)->U_Disperse+index);
								*((((zone)->Iface+iface)->ZSlab)->V_Dmass+index) = *((((zone)->Iface+iface)->ZSlab)->V_Disperse+index);
								*((((zone)->Iface+iface)->ZSlab)->W_Dmass+index) = *((((zone)->Iface+iface)->ZSlab)->W_Disperse+index);
							}
						}
						else if(base->flagMPFV==1)
						{
							*((((zone)->Iface+iface)->ZSlab)->Momentumxstar+index) = *((((zone)->Iface+iface)->ZSlab)->Momentumx+index);
							*((((zone)->Iface+iface)->ZSlab)->Momentumystar+index) = *((((zone)->Iface+iface)->ZSlab)->Momentumy+index);
							*((((zone)->Iface+iface)->ZSlab)->Momentumzstar+index) = *((((zone)->Iface+iface)->ZSlab)->Momentumz+index);
						}
					}
				
				}
			}
			for(zBc=0;zBc<numBc;zBc++)
			{
				if(((zone)->ZoneBc+zBc)->Direction==dir)
				{
					for(k=0;k<NK-1;k++)
					for(j=0;j<NJ-1;j++)
					for(i=0;i<NI-1;i++)
					{
						if(dir==1||dir==2)
							index= k*(NJ-1)+j;
						if(dir==3||dir==4)
							index= k*(NI-1)+i; 
						if(dir==5||dir==6)
							index= j*(NI-1)+i;
						if(base->flagMPFV==0)
						{
							*((((zone)->ZoneBc+zBc)->ZSlab)->umass+index) = *((((zone)->ZoneBc+zBc)->ZSlab)->u+index);
							*((((zone)->ZoneBc+zBc)->ZSlab)->vmass+index) = *((((zone)->ZoneBc+zBc)->ZSlab)->v+index);
							*((((zone)->ZoneBc+zBc)->ZSlab)->wmass+index) = *((((zone)->ZoneBc+zBc)->ZSlab)->w+index);

							if(base->flagDroplet==1 && base->BaseParameters->flagVariableDia==1)
							{	
								*((((zone)->ZoneBc+zBc)->ZSlab)->U_Dmass+index) = *((((zone)->ZoneBc+zBc)->ZSlab)->U_Disperse+index);
								*((((zone)->ZoneBc+zBc)->ZSlab)->V_Dmass+index) = *((((zone)->ZoneBc+zBc)->ZSlab)->V_Disperse+index);
								*((((zone)->ZoneBc+zBc)->ZSlab)->W_Dmass+index) = *((((zone)->ZoneBc+zBc)->ZSlab)->W_Disperse+index);	
							}
						}
						else if(base->flagMPFV==1)
						{
							*((((zone)->ZoneBc+zBc)->ZSlab)->Momentumxstar+index) = *((((zone)->ZoneBc+zBc)->ZSlab)->Momentumx+index);
							*((((zone)->ZoneBc+zBc)->ZSlab)->Momentumystar+index) = *((((zone)->ZoneBc+zBc)->ZSlab)->Momentumy+index);
							*((((zone)->ZoneBc+zBc)->ZSlab)->Momentumzstar+index) = *((((zone)->ZoneBc+zBc)->ZSlab)->Momentumz+index);
						}
					}
				}
			}
		}/*End of direction Loop*/
	}/*End of Zone loop*/

	nZones=base->Number_of_zones;
	for(numZones=0;numZones<nZones;numZones++)
	{
		zone=(base->CgnsZone+numZones);
		
		ptrSize=zone->Size;

		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);

		F = zone->Flowsol->Massflux;

		SurfArea=zone->ZGeom->SurfaceArea;
		if(base->flagMPFV==0)
		{
			domainUmass=MakeDomainArray(zone,100);
			domainVmass=MakeDomainArray(zone,200);
			domainWmass=MakeDomainArray(zone,300);
			domainVolume=MakeDomainArray(zone,5);
	
			MakeSlabPtrArray(zone,slabUmass,100);
			MakeSlabPtrArray(zone,slabVmass,200);
			MakeSlabPtrArray(zone,slabWmass,300);
			MakeSlabPtrArray(zone,slabVolume,5);

			if(base->flagDroplet==1 && base->BaseParameters->flagVariableDia==1)
			{
				domainUDmass=MakeDomainArray(zone,6200);
				domainVDmass=MakeDomainArray(zone,6300);
				domainWDmass=MakeDomainArray(zone,6400);
					
				MakeSlabPtrArray(zone,slabUDmass,6200);
				MakeSlabPtrArray(zone,slabVDmass,6300);
				MakeSlabPtrArray(zone,slabWDmass,6400);
			}
		}
		else if(base->flagMPFV==1)
		{
		
			domainUmass=MakeDomainArray(zone,3100);
			domainVmass=MakeDomainArray(zone,3200);
			domainWmass=MakeDomainArray(zone,3300);
			domainVolume=MakeDomainArray(zone,5);
			domainRho=MakeDomainArray(zone,710);
			
			MakeSlabPtrArray(zone,slabUmass,3100);
			MakeSlabPtrArray(zone,slabVmass,3200);
			MakeSlabPtrArray(zone,slabWmass,3300);
			MakeSlabPtrArray(zone,slabVolume,5);
			MakeSlabPtrArray(zone,slabRho,710);
		
		}
		if ((base->flagParticle==1)|| (base->flagDroplet==1))
		{
			FD=zone->Flowsol->MassfluxDisperse;
			
			domainUDisperse=MakeDomainArray(zone,62);
			domainVDisperse=MakeDomainArray(zone,63);
			domainWDisperse=MakeDomainArray(zone,64);
			domainVolume=MakeDomainArray(zone,5);
			//domainRhoDisperse=MakeDomainArray(zone,66);
			//domainVartheta=MakeDomainArray(zone,65);
		
			MakeSlabPtrArray(zone,slabUDisperse,62);
			MakeSlabPtrArray(zone,slabVDisperse,63);
			MakeSlabPtrArray(zone,slabWDisperse,64);
			MakeSlabPtrArray(zone,slabVolume,5);
			//MakeSlabPtrArray(zone,slabRhoDisperse,66);
			//MakeSlabPtrArray(zone,slabVartheta,65);
		}

		for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{
			indexP=(NJ-1)*(NI-1)*k+(NI-1)*j+i;
			indexE=(NJ-1)*(NI-1)*k+(NI-1)*j+(i+1);
			indexW=(NJ-1)*(NI-1)*k+(NI-1)*j+(i-1);
			indexN=(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i;
			indexS=(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i;
			indexT=(NJ-1)*(NI-1)*(k+1)+(NI-1)*j+i;
			indexB=(NJ-1)*(NI-1)*(k-1)+(NI-1)*j+i;

			u_P=*(domainUmass+indexP);
			v_P=*(domainVmass+indexP);
			w_P=*(domainWmass+indexP);
			vol_P=*(domainVolume+indexP);
			
			if ((base->flagParticle==1)|| (base->flagDroplet==1 &&(base->BaseParameters->flagVariableDia==0)))
			{
				u_PD=*(domainUDisperse+indexP);
				v_PD=*(domainVDisperse+indexP);
				w_PD=*(domainWDisperse+indexP);
			}
			if ((base->flagParticle==1)|| (base->flagDroplet==1 &&(base->BaseParameters->flagVariableDia==1)))
			{
				u_PD=*(domainUDmass+indexP);
				v_PD=*(domainVDmass+indexP);
				w_PD=*(domainWDmass+indexP);
			}
			
		/*	if((base->flagParticle==1)|| (base->flagDroplet==1))
			{
				rho_PD=1.0/(*(domainRhoDisperse+indexP)*(*(domainVartheta+indexP)));
				u_PD = u_PD*rho_PD;
				v_PD = v_PD*rho_PD;
				w_PD = w_PD*rho_PD;
				
			}	*/
			if(base->flagMPFV==1)
			{
				rho_P=1/(*(domainRho+indexP));
				u_P = u_P*rho_P;
				v_P = v_P*rho_P;
				w_P = w_P*rho_P;
				
			}
			if(i==(NI-2))
			{
				u_E=*(slabUmass[0]+(NJ-1)*k+j);
				v_E=*(slabVmass[0]+(NJ-1)*k+j);
				w_E=*(slabWmass[0]+(NJ-1)*k+j);
				vol_E=*(slabVolume[0]+(NJ-1)*k+j);
				
				if ((base->flagParticle==1)|| (base->flagDroplet==1 && (base->BaseParameters->flagVariableDia==0)))
				{
					u_ED=*(slabUDisperse[0]+(NJ-1)*k+j);
					v_ED=*(slabVDisperse[0]+(NJ-1)*k+j);
					w_ED=*(slabWDisperse[0]+(NJ-1)*k+j);
				/*	rho_ED = 1.0/(*(slabRhoDisperse[0]+(NJ-1)*k+j)*(*(slabVartheta[0]+(NJ-1)*k+j)));
					u_ED = u_ED*rho_ED;
					v_ED = v_ED*rho_ED;
					w_ED = w_ED*rho_ED;	*/
				}
				if ((base->flagParticle==1)|| (base->flagDroplet==1 && (base->BaseParameters->flagVariableDia==1)))
				{
					u_ED=*(slabUDmass[0]+(NJ-1)*k+j);
					v_ED=*(slabVDmass[0]+(NJ-1)*k+j);
					w_ED=*(slabWDmass[0]+(NJ-1)*k+j);
				}
				if(base->flagMPFV==1)
				{
					rho_E = 1/(*(slabRho[0]+(NJ-1)*k+j));
					u_E = u_E*rho_E;
					v_E = v_E*rho_E;
					w_E = w_E*rho_E;
				}
				
			}
			else
			{
				u_E=*(domainUmass+indexE);
				v_E=*(domainVmass+indexE);
				w_E=*(domainWmass+indexE);
				vol_E=*(domainVolume+indexE);
				if ((base->flagParticle==1)|| (base->flagDroplet==1 && (base->BaseParameters->flagVariableDia==0)))
				{	
					u_ED=*(domainUDisperse+indexE);
					v_ED=*(domainVDisperse+indexE);
					w_ED=*(domainWDisperse+indexE);
				/*	rho_ED = 1.0/(*(domainRhoDisperse+indexE)*(*(domainVartheta+indexE)));
					u_ED = u_ED*rho_ED;
					v_ED = v_ED*rho_ED;
					w_ED = w_ED*rho_ED;	*/
					
				}
				if ((base->flagParticle==1)|| (base->flagDroplet==1 && (base->BaseParameters->flagVariableDia==1)))
				{
					u_ED=*(domainUDmass+indexE);
					v_ED=*(domainVDmass+indexE);
					w_ED=*(domainWDmass+indexE);
				}
				if(base->flagMPFV==1)
				{
					rho_E = 1/(*(domainRho+indexE));
					u_E = u_E*rho_E;
					v_E = v_E*rho_E;
					w_E = w_E*rho_E;
					
				}
			}

			if(i==0)
			{
				
				u_W =*(slabUmass[1]+(NJ-1)*k+j);
				v_W =*(slabVmass[1]+(NJ-1)*k+j);
				w_W =*(slabWmass[1]+(NJ-1)*k+j);
				vol_W =*(slabVolume[1]+(NJ-1)*k+j);
				if ((base->flagParticle==1)|| (base->flagDroplet==1 && (base->BaseParameters->flagVariableDia==0)))
				{
					u_WD =*(slabUDisperse[1]+(NJ-1)*k+j);
					v_WD =*(slabVDisperse[1]+(NJ-1)*k+j);
					w_WD =*(slabWDisperse[1]+(NJ-1)*k+j);
				/*	rho_WD = 1.0/(*(slabRhoDisperse[1]+(NJ-1)*k+j)*(*(slabVartheta[1]+(NJ-1)*k+j)));
					u_WD = u_WD*rho_WD;
					v_WD = v_WD*rho_WD;
					w_WD = w_WD*rho_WD;	*/
				}
				if ((base->flagParticle==1)|| (base->flagDroplet==1 && (base->BaseParameters->flagVariableDia==1)))
				{
					u_WD =*(slabUDmass[1]+(NJ-1)*k+j);
					v_WD =*(slabVDmass[1]+(NJ-1)*k+j);
					w_WD =*(slabWDmass[1]+(NJ-1)*k+j);
				}
				if(base->flagMPFV==1)
				{
					rho_W = 1/(*(slabRho[1]+(NJ-1)*k+j));
					u_W = u_W*rho_W;
					v_W = v_W*rho_W;
					w_W = w_W*rho_W;
				}
				
			}
			
			else
			{
				u_W=*(domainUmass+indexW);
				v_W=*(domainVmass+indexW);
				w_W=*(domainWmass+indexW);
				vol_W=*(domainVolume+indexW);
				if ((base->flagParticle==1)|| (base->flagDroplet==1 && (base->BaseParameters->flagVariableDia==0)))
				{	
					u_WD=*(domainUDisperse+indexW);
					v_WD=*(domainVDisperse+indexW);
					w_WD=*(domainWDisperse+indexW);
				/*	rho_WD = 1.0/(*(domainRhoDisperse+indexW)*(*(domainVartheta+indexW)));
					u_WD = u_WD*rho_WD;
					v_WD = v_WD*rho_WD;
					w_WD = w_WD*rho_WD;	*/
				}
				if ((base->flagParticle==1)|| (base->flagDroplet==1 && (base->BaseParameters->flagVariableDia==1)))
				{
					u_WD =*(slabUDmass[1]+(NJ-1)*k+j);
					v_WD =*(slabVDmass[1]+(NJ-1)*k+j);
					w_WD =*(slabWDmass[1]+(NJ-1)*k+j);
				}
				if(base->flagMPFV==1)
				{
					rho_W = 1/(*(domainRho+indexW));
					u_W = u_W*rho_W;
					v_W = v_W*rho_W;
					w_W = w_W*rho_W;
				}
			}

			if(j==(NJ-2))
			{
				
				u_N=*(slabUmass[2]+(NI-1)*k+i);
				v_N=*(slabVmass[2]+(NI-1)*k+i);
				w_N=*(slabWmass[2]+(NI-1)*k+i);
				vol_N=*(slabVolume[2]+(NI-1)*k+i);
				if ((base->flagParticle==1)|| (base->flagDroplet==1 && (base->BaseParameters->flagVariableDia==0)))
				{
					u_ND=*(slabUDisperse[2]+(NI-1)*k+i);
					v_ND=*(slabVDisperse[2]+(NI-1)*k+i);
					w_ND=*(slabWDisperse[2]+(NI-1)*k+i);
				/*	rho_ND = 1.0/(*(slabRhoDisperse[2]+(NI-1)*k+i)*(*(slabVartheta[2]+(NI-1)*k+i)));
					u_ND = u_ND*rho_ND;
					v_ND = v_ND*rho_ND;
					w_ND = w_ND*rho_ND;	*/
				}
				if ((base->flagParticle==1)|| (base->flagDroplet==1 && (base->BaseParameters->flagVariableDia==1)))
				{
					u_ND=*(slabUDmass[2]+(NI-1)*k+i);
					v_ND=*(slabVDmass[2]+(NI-1)*k+i);
					w_ND=*(slabWDmass[2]+(NI-1)*k+i);
				}
				if(base->flagMPFV==1)
				{
					rho_N = 1/(*(slabRho[2]+(NI-1)*k+i));
					u_N = u_N*rho_N;
					v_N = v_N*rho_N;
					w_N = w_N*rho_N;
				}
				
			}
			else
			{
				u_N=*(domainUmass+indexN);
				v_N=*(domainVmass+indexN);
				w_N=*(domainWmass+indexN);
				vol_N=*(domainVolume+indexN);
				if ((base->flagParticle==1)|| (base->flagDroplet==1 && (base->BaseParameters->flagVariableDia==0)))
				{
					u_ND=*(domainUDisperse+indexN);
					v_ND=*(domainVDisperse+indexN);
					w_ND=*(domainWDisperse+indexN);
				/*	rho_ND = 1/(*(domainRhoDisperse+indexN)*(*(domainVartheta+indexN)));
					u_ND = u_ND*rho_ND;
					v_ND = v_ND*rho_ND;
					w_ND = w_ND*rho_ND; */
				}
				if ((base->flagParticle==1)|| (base->flagDroplet==1 && (base->BaseParameters->flagVariableDia==1)))
				{
					u_ND=*(domainUDmass+indexN);
					v_ND=*(domainVDmass+indexN);
					w_ND=*(domainWDmass+indexN);
				}
				if(base->flagMPFV==1)
				{
					rho_N = 1/(*(domainRho+indexN));
					u_N = u_N*rho_N;
					v_N = v_N*rho_N;
					w_N = w_N*rho_N;
				}
			}

			if(j==0)
			{	
				
				u_S=*(slabUmass[3]+(NI-1)*k+i);
				v_S=*(slabVmass[3]+(NI-1)*k+i);
				w_S=*(slabWmass[3]+(NI-1)*k+i);
				vol_S=*(slabVolume[3]+(NI-1)*k+i);
				if ((base->flagParticle==1)|| (base->flagDroplet==1 && (base->BaseParameters->flagVariableDia==0)))
				{
					u_SD=*(slabUDisperse[3]+(NI-1)*k+i);
					v_SD=*(slabVDisperse[3]+(NI-1)*k+i);
					w_SD=*(slabWDisperse[3]+(NI-1)*k+i);
				/*	rho_SD = 1.0/(*(slabRhoDisperse[3]+(NI-1)*k+i)*(*(slabVartheta[3]+(NI-1)*k+i)));
					u_SD = u_SD*rho_SD;
					v_SD = v_SD*rho_SD;
					w_SD = w_SD*rho_SD;	*/
				}
				if ((base->flagParticle==1)|| (base->flagDroplet==1 && (base->BaseParameters->flagVariableDia==1)))
				{
					u_SD=*(slabUDmass[3]+(NI-1)*k+i);
					v_SD=*(slabVDmass[3]+(NI-1)*k+i);
					w_SD=*(slabWDmass[3]+(NI-1)*k+i);
				}
				if(base->flagMPFV==1)
				{
					rho_S = 1/(*(slabRho[3]+(NI-1)*k+i));
					u_S = u_S*rho_S;
					v_S = v_S*rho_S;
					w_S = w_S*rho_S;
				}
				
			}
			else
			{
				u_S=*(domainUmass+indexS);
				v_S=*(domainVmass+indexS);
				w_S=*(domainWmass+indexS);
				vol_S=*(domainVolume+indexS);
				if ((base->flagParticle==1)|| (base->flagDroplet==1 && (base->BaseParameters->flagVariableDia==0)))
				{
					u_SD=*(domainUDisperse+indexS);
					v_SD=*(domainVDisperse+indexS);
					w_SD=*(domainWDisperse+indexS);
				/*	rho_SD = 1.0/(*(domainRhoDisperse+indexS)*(*(domainVartheta+indexS)));
					u_SD = u_SD*rho_SD;
					v_SD = v_SD*rho_SD;
					w_SD = w_SD*rho_SD; */
				}
				if ((base->flagParticle==1)|| (base->flagDroplet==1 && (base->BaseParameters->flagVariableDia==1)))
				{
					u_SD=*(domainUDmass+indexS);
					v_SD=*(domainVDmass+indexS);
					w_SD=*(domainWDmass+indexS);
				}
				if(base->flagMPFV==1)
				{
					rho_S = 1/(*(domainRho+indexS));
					u_S = u_S*rho_S;
					v_S = v_S*rho_S;
					w_S = w_S*rho_S;
				}
			}

			if(k==(NK-2))
			{	
				
				u_T=*(slabUmass[4]+(NI-1)*j+i);
				v_T=*(slabVmass[4]+(NI-1)*j+i);
				w_T=*(slabWmass[4]+(NI-1)*j+i);
				vol_T=*(slabVolume[4]+(NI-1)*j+i);
				if ((base->flagParticle==1)|| (base->flagDroplet==1 && (base->BaseParameters->flagVariableDia==0)))
				{
					u_TD=*(slabUDisperse[4]+(NI-1)*j+i);
					v_TD=*(slabVDisperse[4]+(NI-1)*j+i);
					w_TD=*(slabWDisperse[4]+(NI-1)*j+i);
				/*	rho_TD = 1.0/(*(slabRhoDisperse[4]+(NI-1)*j+i)*(*(slabRho[4]+(NI-1)*j+i)));
					u_TD = u_TD*rho_TD;
					v_TD = v_TD*rho_TD;
					w_TD = w_TD*rho_TD; */
				}
				vol_T=*(slabVolume[4]+(NI-1)*j+i);
				if ((base->flagParticle==1)|| (base->flagDroplet==1 && (base->BaseParameters->flagVariableDia==1)))
				{
					u_TD=*(slabUDmass[4]+(NI-1)*j+i);
					v_TD=*(slabVDmass[4]+(NI-1)*j+i);
					w_TD=*(slabWDmass[4]+(NI-1)*j+i);
				}
				if(base->flagMPFV==1)
				{
					rho_T = 1/(*(slabRho[4]+(NI-1)*j+i));
					u_T = u_T*rho_T;
					v_T = v_T*rho_T;
					w_T = w_T*rho_T;
				}
			}
			else
			{
				u_T=*(domainUmass+indexT);
				v_T=*(domainVmass+indexT);
				w_T=*(domainWmass+indexT);
				vol_T=*(domainVolume+indexT);
				if ((base->flagParticle==1)|| (base->flagDroplet==1 && (base->BaseParameters->flagVariableDia==0)))
				{
					u_TD=*(domainUDisperse+indexT);
					v_TD=*(domainVDisperse+indexT);
					w_TD=*(domainWDisperse+indexT);
				/*	rho_TD = 1.0/(*(domainRhoDisperse+indexT)*(*(domainVartheta+indexT)));
					u_TD = u_TD*rho_TD;
					v_TD = v_TD*rho_TD;
					w_TD = w_TD*rho_TD;*/
				}
				if ((base->flagParticle==1)|| (base->flagDroplet==1 && (base->BaseParameters->flagVariableDia==1)))
				{
					u_TD=*(domainUDmass+indexT);
					v_TD=*(domainVDmass+indexT);
					w_TD=*(domainWDmass+indexT);
				}
				if(base->flagMPFV==1)
				{
					rho_T = 1/(*(domainRho+indexT));
					u_T = u_T*rho_T;
					v_T = v_T*rho_T;
					w_T = w_T*rho_T;
				}
			}

			if(k==0)
			{
				u_B=*(slabUmass[5]+(NI-1)*j+i);
				v_B=*(slabVmass[5]+(NI-1)*j+i);
				w_B=*(slabWmass[5]+(NI-1)*j+i);
				vol_B=*(slabVolume[5]+(NI-1)*j+i);
				if ((base->flagParticle==1)|| (base->flagDroplet==1 && (base->BaseParameters->flagVariableDia==0)))
				{
					u_BD=*(slabUDisperse[5]+(NI-1)*j+i);
					v_BD=*(slabVDisperse[5]+(NI-1)*j+i);
					w_BD=*(slabWDisperse[5]+(NI-1)*j+i);
				/*	rho_BD = 1.0/(*(slabRhoDisperse[5]+(NI-1)*j+i)*(*(slabVartheta[5]+(NI-1)*j+i)));
					u_BD = u_BD*rho_BD;
					v_BD = v_BD*rho_BD;
					w_BD = w_BD*rho_BD; */
				}
				if ((base->flagParticle==1)|| (base->flagDroplet==1 && (base->BaseParameters->flagVariableDia==1)))
				{
					u_BD=*(slabUDmass[5]+(NI-1)*j+i);
					v_BD=*(slabVDmass[5]+(NI-1)*j+i);
					w_BD=*(slabWDmass[5]+(NI-1)*j+i);
				}
				if(base->flagMPFV==1)
				{
					rho_B = 1/(*(slabRho[5]+(NI-1)*j+i));
					u_B = u_B*rho_B;
					v_B = v_B*rho_B;
					w_B = w_B*rho_B;
				}
			}
			else
			{
				u_B=*(domainUmass+indexB);
				v_B=*(domainVmass+indexB);
				w_B=*(domainWmass+indexB);
				vol_B=*(domainVolume+indexB);
				if ((base->flagParticle==1)|| (base->flagDroplet==1 && (base->BaseParameters->flagVariableDia==0)))
				{
					u_BD=*(domainUDisperse+indexB);
					v_BD=*(domainVDisperse+indexB);
					w_BD=*(domainWDisperse+indexB);
				/*	rho_BD = 1.0/(*(domainRhoDisperse+indexB)*(*(domainVartheta+indexB)));
					u_BD = u_BD*rho_BD;
					v_BD = v_BD*rho_BD;
					w_BD = w_BD*rho_BD; */
				}
				if ((base->flagParticle==1)|| (base->flagDroplet==1 && (base->BaseParameters->flagVariableDia==1)))
				{
					u_BD=*(domainUDmass+indexB);
					v_BD=*(domainVDmass+indexB);
					w_BD=*(domainWDmass+indexB);
				}
				if(base->flagMPFV==1)
				{
					rho_B = 1/(*(domainRho+indexB));
					u_B = u_B*rho_B;
					v_B = v_B*rho_B;
					w_B = w_B*rho_B;
				}
			}

			ue=((u_E*vol_P+u_P*vol_E)/(vol_P+vol_E));
			uw=((u_W*vol_P+u_P*vol_W)/(vol_P+vol_W));
			us=((u_S*vol_P+u_P*vol_S)/(vol_P+vol_S));
			un=((u_N*vol_P+u_P*vol_N)/(vol_P+vol_N));
			ut=((u_T*vol_P+u_P*vol_T)/(vol_P+vol_T));
			ub=((u_B*vol_P+u_P*vol_B)/(vol_P+vol_B));

			ve=((v_E*vol_P+v_P*vol_E)/(vol_P+vol_E));
			vw=((v_W*vol_P+v_P*vol_W)/(vol_P+vol_W));
			vs=((v_S*vol_P+v_P*vol_S)/(vol_P+vol_S));
			vn=((v_N*vol_P+v_P*vol_N)/(vol_P+vol_N));
			vt=((v_T*vol_P+v_P*vol_T)/(vol_P+vol_T));
			vb=((v_B*vol_P+v_P*vol_B)/(vol_P+vol_B));

			we=((w_E*vol_P+w_P*vol_E)/(vol_P+vol_E));
			ww=((w_W*vol_P+w_P*vol_W)/(vol_P+vol_W));
			ws=((w_S*vol_P+w_P*vol_S)/(vol_P+vol_S));
			wn=((w_N*vol_P+w_P*vol_N)/(vol_P+vol_N));
			wt=((w_T*vol_P+w_P*vol_T)/(vol_P+vol_T));
			wb=((w_B*vol_P+w_P*vol_B)/(vol_P+vol_B));

			Index=((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6*3;
			face = 0;  // west face
				*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = (uw*(*(SurfArea+Index+3*0+0))+vw*(*(SurfArea+Index+3*0+1))+ww*(*(SurfArea+Index+3*0+2)));
				//fprintf(fp1,"west=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
			face = 1;  // east face
				*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = (ue*(*(SurfArea+Index+3*1+0))+ve*(*(SurfArea+Index+3*1+1))+we*(*(SurfArea+Index+3*1+2)));
				//fprintf(fp1,"east=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
			face = 2;  // south face
				*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = (us*(*(SurfArea+Index+3*2+0))+vs*(*(SurfArea+Index+3*2+1))+ws*(*(SurfArea+Index+3*2+2)));
				//fprintf(fp1,"south=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
			face = 3;  // north face
				*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = (un*(*(SurfArea+Index+3*3+0))+vn*(*(SurfArea+Index+3*3+1))+wn*(*(SurfArea+Index+3*3+2)));
				//fprintf(fp1,"north=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
			face = 4;  // bottom face
				*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = (ub*(*(SurfArea+Index+3*4+0))+vb*(*(SurfArea+Index+3*4+1))+wb*(*(SurfArea+Index+3*4+2)));
				//fprintf(fp1,"bottom=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
			face = 5;  // top face
				*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = (ut*(*(SurfArea+Index+3*5+0))+vt*(*(SurfArea+Index+3*5+1))+wt*(*(SurfArea+Index+3*5+2)));
				//fprintf(fp1,"top=%lf\n",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
			
			//******************added by adhiraj for disperse phase************************
			if ((base->flagParticle==1)|| (base->flagDroplet==1))
			{
			ueD=((u_ED*vol_P+u_PD*vol_E)/(vol_P+vol_E));
			uwD=((u_WD*vol_P+u_PD*vol_W)/(vol_P+vol_W));
			usD=((u_SD*vol_P+u_PD*vol_S)/(vol_P+vol_S));
			unD=((u_ND*vol_P+u_PD*vol_N)/(vol_P+vol_N));
			utD=((u_TD*vol_P+u_PD*vol_T)/(vol_P+vol_T));
			ubD=((u_BD*vol_P+u_PD*vol_B)/(vol_P+vol_B));
					
			veD=((v_ED*vol_P+v_PD*vol_E)/(vol_P+vol_E));
			vwD=((v_WD*vol_P+v_PD*vol_W)/(vol_P+vol_W));
			vsD=((v_SD*vol_P+v_PD*vol_S)/(vol_P+vol_S));
			vnD=((v_ND*vol_P+v_PD*vol_N)/(vol_P+vol_N));
			vtD=((v_TD*vol_P+v_PD*vol_T)/(vol_P+vol_T));
			vbD=((v_B*vol_P+v_P*vol_B)/(vol_P+vol_B));

			weD=((w_ED*vol_P+w_PD*vol_E)/(vol_P+vol_E));
			wwD=((w_WD*vol_P+w_PD*vol_W)/(vol_P+vol_W));
			wsD=((w_SD*vol_P+w_PD*vol_S)/(vol_P+vol_S));
			wnD=((w_ND*vol_P+w_PD*vol_N)/(vol_P+vol_N));
			wtD=((w_TD*vol_P+w_PD*vol_T)/(vol_P+vol_T));
			wbD=((w_BD*vol_P+w_PD*vol_B)/(vol_P+vol_B));

			Index=((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6*3;
			face = 0;  // west face
				*(FD+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = (uwD*(*(SurfArea+Index+3*0+0))+vwD*(*(SurfArea+Index+3*0+1))+wwD*(*(SurfArea+Index+3*0+2)));
				//fprintf(fp1,"west=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
			face = 1;  // east face
				*(FD+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = (ueD*(*(SurfArea+Index+3*1+0))+veD*(*(SurfArea+Index+3*1+1))+weD*(*(SurfArea+Index+3*1+2)));
				//fprintf(fp1,"east=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
			face = 2;  // south face
				*(FD+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = (usD*(*(SurfArea+Index+3*2+0))+vsD*(*(SurfArea+Index+3*2+1))+wsD*(*(SurfArea+Index+3*2+2)));
				//fprintf(fp1,"south=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
			face = 3;  // north face
				*(FD+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = (unD*(*(SurfArea+Index+3*3+0))+vnD*(*(SurfArea+Index+3*3+1))+wnD*(*(SurfArea+Index+3*3+2)));
				//fprintf(fp1,"north=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
			face = 4;  // bottom face
				*(FD+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = (ubD*(*(SurfArea+Index+3*4+0))+vbD*(*(SurfArea+Index+3*4+1))+wbD*(*(SurfArea+Index+3*4+2)));
				//fprintf(fp1,"bottom=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
			face = 5;  // top face
				*(FD+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = (utD*(*(SurfArea+Index+3*5+0))+vtD*(*(SurfArea+Index+3*5+1))+wtD*(*(SurfArea+Index+3*5+2)));
			
				}
		
		
		}/*End of looping over NI,NJ,NK*/
	}/*End of looping over the zones */
	return 0;
}

/**
This function is used only once for calculating the volumeflux just before the first timestep using
the guess values of mass velocities. This is used for 2D geometries.  This also calculates the volume flux of disperse phase that is either particle or droplet phase if disperse phase flow is choosen by the user.
*/
int volumeflux2d(Base *base)
{
	int nZones,numZones;
	int k,j,i,index,indexB,indexT,indexE,indexW,indexS,indexN,indexP;
	int numiface,numBc,face,iface,zBc,dir,Index;
	double *SurfArea,*F,*FD;
	double *domainVolume,*domainUmass,*domainVmass,*domainRho;
	double *slabVolume[6],*slabUmass[6],*slabVmass[6],*slabRho[6];
	Zone* zone;
	
	double *domainRhoDisperse,*domainVartheta,*domainUDisperse,*domainVDisperse;
	double *slabVartheta[6],*slabUDisperse[6],*slabVDisperse[6],*slabRhoDisperse[6];
	double *slabUDmass[6], *slabVDmass[6], *domainUDmass, *domainVDmass;
	
	int NI,NJ,NK,*ptrSize;
	double uw,vw,ue,ve,un,vn,us,vs,ut,vt,ub,vb;
	double u_P,u_E,u_W,u_S,u_N,u_T,u_B,v_P,v_E,v_W,v_S,v_N,v_T,v_B;
	double uwD,vwD,ueD,veD,unD,vnD,usD,vsD,utD,vtD,ubD,vbD;
	double u_PD,u_ED,u_WD,u_SD,u_ND,u_TD,u_BD,v_PD,v_ED,v_WD,v_SD,v_ND,v_TD,v_BD;
	double vol_P,vol_E,vol_W,vol_N,vol_S,vol_T,vol_B,rho_P,rho_E,rho_W,rho_S,rho_N,rho_T,rho_B;
	
	double rho_PD,rho_ED,rho_WD,rho_SD,rho_ND,rho_TD,rho_BD;
	int flagrhodisperse;
	
	//shash
	double EPSILON;
	EPSILON = base->BaseParameters->EPSILONDrop;
	int che;
	//FILE *fp;
	//fp = fopen("volflux.dat","w");	
	nZones=base->Number_of_zones;	
	

	for(numZones=0;numZones<nZones;numZones++)
	{
		/*Copying down u,v,w into umass,vmass,wmass ....domain as well as slabs*/
		ptrSize=((base)->CgnsZone+numZones)->Size;
		zone=(base->CgnsZone+numZones);
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		
		for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{
			indexP=(NJ-1)*(NI-1)*k+(NI-1)*j+i;
			if(base->flagMPFV==0)
			{
				*((zone->Flowsol)->umass+indexP)=*((zone->Flowsol)->u+indexP);
				*((zone->Flowsol)->vmass+indexP)=*((zone->Flowsol)->v+indexP);

				if(base->flagDroplet==1 && base->BaseParameters->flagVariableDia==1)
				{
					*((zone->Flowsol)->U_Dmass+indexP)=*((zone->Flowsol)->U_Disperse+indexP);
					*((zone->Flowsol)->V_Dmass+indexP)=*((zone->Flowsol)->V_Disperse+indexP);
				}

			}
			else if(base->flagMPFV==1)
			{
				*((zone->Flowsol)->Momentumxstar+indexP)=*((zone->Flowsol)->Momentumx+indexP);
				*((zone->Flowsol)->Momentumystar+indexP)=*((zone->Flowsol)->Momentumy+indexP);
			}
		}
		numiface=(zone)->Number_of_Interfaces;
		
		numBc=(zone)->Number_of_BC_Patches;
		
		for(dir=1;dir<=6;dir++)
		{
			for(iface=0;iface<numiface;iface++)
			{
				if(((zone)->Iface+iface)->Direction==dir)
				{
					for(k=0;k<NK-1;k++)
					for(j=0;j<NJ-1;j++)
					for(i=0;i<NI-1;i++)
					{
						if(dir==1||dir==2)
							index= k*(NJ-1)+j;
						if(dir==3||dir==4)
							index= k*(NI-1)+i;
						if(dir==5||dir==6)
							index= j*(NI-1)+i;
						if(base->flagMPFV==0)
						{
							*((((zone)->Iface+iface)->ZSlab)->umass+index) = *((((zone)->Iface+iface)->ZSlab)->u+index);
							*((((zone)->Iface+iface)->ZSlab)->vmass+index) = *((((zone)->Iface+iface)->ZSlab)->v+index);

							if(base->flagDroplet==1 && base->BaseParameters->flagVariableDia==1)
							{
								*((((zone)->Iface+iface)->ZSlab)->U_Dmass+index) = *((((zone)->Iface+iface)->ZSlab)->U_Disperse+index);
								*((((zone)->Iface+iface)->ZSlab)->V_Dmass+index) = *((((zone)->Iface+iface)->ZSlab)->V_Disperse+index);
							}
				
						}
						else if(base->flagMPFV==1)
						{
							*((((zone)->Iface+iface)->ZSlab)->Momentumxstar+index) = *((((zone)->Iface+iface)->ZSlab)->Momentumx+index);
							*((((zone)->Iface+iface)->ZSlab)->Momentumystar+index) = *((((zone)->Iface+iface)->ZSlab)->Momentumy+index);
						}
					}
				
				}
			}
			for(zBc=0;zBc<numBc;zBc++)
			{
				if(((zone)->ZoneBc+zBc)->Direction==dir)
				{
					for(k=0;k<NK-1;k++)
					for(j=0;j<NJ-1;j++)
					for(i=0;i<NI-1;i++)
					{
						if(dir==1||dir==2)
							index= k*(NJ-1)+j;
						if(dir==3||dir==4)
							index= k*(NI-1)+i; 
						if(dir==5||dir==6)
							index= j*(NI-1)+i;
						if(base->flagMPFV==0)
						{
							*((((zone)->ZoneBc+zBc)->ZSlab)->umass+index) = *((((zone)->ZoneBc+zBc)->ZSlab)->u+index);
							*((((zone)->ZoneBc+zBc)->ZSlab)->vmass+index) = *((((zone)->ZoneBc+zBc)->ZSlab)->v+index);

							if(base->flagDroplet==1 && base->BaseParameters->flagVariableDia==1)
							{	
								*((((zone)->ZoneBc+zBc)->ZSlab)->U_Dmass+index) = *((((zone)->ZoneBc+zBc)->ZSlab)->U_Disperse+index);
								*((((zone)->ZoneBc+zBc)->ZSlab)->V_Dmass+index) = *((((zone)->ZoneBc+zBc)->ZSlab)->V_Disperse+index);
							}	
						}
						else if(base->flagMPFV==1)
						{
						
							*((((zone)->ZoneBc+zBc)->ZSlab)->Momentumxstar+index) = *((((zone)->ZoneBc+zBc)->ZSlab)->Momentumx+index);
							*((((zone)->ZoneBc+zBc)->ZSlab)->Momentumystar+index) = *((((zone)->ZoneBc+zBc)->ZSlab)->Momentumy+index);
						}
					}
				}
			}
		}/*End of direction Loop*/
	}/*End of Zone loop*/
	

	nZones=base->Number_of_zones;
	for(numZones=0;numZones<nZones;numZones++)
	{
		zone=(base->CgnsZone+numZones);
		
		ptrSize=zone->Size;

		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);

		F = zone->Flowsol->Massflux;

		SurfArea=zone->ZGeom->SurfaceArea;
		if(base->flagMPFV==0)
		{
			domainUmass=MakeDomainArray(zone,100);
			domainVmass=MakeDomainArray(zone,200);
			domainVolume=MakeDomainArray(zone,5);
	
			MakeSlabPtrArray(zone,slabUmass,100);
			MakeSlabPtrArray(zone,slabVmass,200);
			MakeSlabPtrArray(zone,slabVolume,5);

			if(base->flagDroplet==1 && base->BaseParameters->flagVariableDia==1)
			{
				FD = zone->Flowsol->MassfluxDisperse;
				domainUDmass=MakeDomainArray(zone,6200);
				domainVDmass=MakeDomainArray(zone,6300);
					
				MakeSlabPtrArray(zone,slabUDmass,6200);
				MakeSlabPtrArray(zone,slabVDmass,6300);
			}
		}
		
		else if(base->flagMPFV==1)
		{
			domainUmass=MakeDomainArray(zone,3100);
			domainVmass=MakeDomainArray(zone,3200);
			domainVolume=MakeDomainArray(zone,5);
			domainRho = MakeDomainArray(zone,710);
	
			MakeSlabPtrArray(zone,slabUmass,3100);
			MakeSlabPtrArray(zone,slabVmass,3200);
			MakeSlabPtrArray(zone,slabVolume,5);
			MakeSlabPtrArray(zone,slabRho,710);
		}
	//NOTE: This function was earlier written to deal with specific momenta from which velocity components was calculated. Since in the new gas=particle and gas-droplet flows solver, we calculate velocity components from the momentum equations of the dispersed-phase, parts of the code has been commented. Refer the thesis on "Development of Finite Volume Solver for Dilute Gas-Droplet flows with evaporation" for the governing equations. 
		
		if ((base->flagParticle==1) || (base->flagDroplet==1 && base->BaseParameters->flagVariableDia==0))
		{
			flagrhodisperse=base->BaseParameters->flagRhoDisperse;
			FD = zone->Flowsol->MassfluxDisperse;
			domainUDisperse=MakeDomainArray(zone,62);
			domainVDisperse=MakeDomainArray(zone,63);
			domainVolume=MakeDomainArray(zone,5);
			//domainRhoDisperse=MakeDomainArray(zone,66);//RhoDisperse
			//domainVartheta=MakeDomainArray(zone,65);
		
			MakeSlabPtrArray(zone,slabUDisperse,62);
			MakeSlabPtrArray(zone,slabVDisperse,63);
			MakeSlabPtrArray(zone,slabVolume,5);
		//	MakeSlabPtrArray(zone,slabRhoDisperse,66);//RhoDisperse
		//	MakeSlabPtrArray(zone,slabVartheta,65);	
		}	
		for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{
			indexP=(NJ-1)*(NI-1)*k+(NI-1)*j+i;
			indexE=(NJ-1)*(NI-1)*k+(NI-1)*j+(i+1);
			indexW=(NJ-1)*(NI-1)*k+(NI-1)*j+(i-1);
			indexN=(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i;
			indexS=(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i;
			indexT=(NJ-1)*(NI-1)*(k+1)+(NI-1)*j+i;
			indexB=(NJ-1)*(NI-1)*(k-1)+(NI-1)*j+i;

			u_P=*(domainUmass+indexP);
			v_P=*(domainVmass+indexP);
			vol_P=*(domainVolume+indexP);
			
			if ((base->flagParticle==1)|| (base->flagDroplet==1 &&(base->BaseParameters->flagVariableDia==0)))
			{
				u_PD=*(domainUDisperse+indexP);
				v_PD=*(domainVDisperse+indexP);
			}
			if (base->flagDroplet==1 &&(base->BaseParameters->flagVariableDia==1))
			{
				u_PD=*(domainUDmass+indexP);
				v_PD=*(domainVDmass+indexP);
			}
//The following lines have been commented as this function is called only once before entering the time stepping loop for the first time-step. At this point, the normalized-volume fraction of dispersed phase will not be less than epsilon. Hence this if condition is not required
			
		/*	if ((base->flagParticle==1)|| (base->flagDroplet==1))
			{
				if(!(*(domainVartheta+indexP) <= EPSILON))
				{
					rho_PD = (*(domainVartheta+indexP) <= EPSILON) ? 0.0 :  (1.0/(*(domainRhoDisperse+indexP*flagrhodisperse)*(*(domainVartheta+indexP)))); 
					
					u_PD = u_PD*rho_PD;
					v_PD = v_PD*rho_PD;
				}
				else
				{
					u_PD = u_P;
					v_PD = v_P;
				}  
			}*/
			if(base->flagMPFV==1)
			{
				rho_P=1/(*(domainRho+indexP));
				u_P = u_P*rho_P;
				v_P = v_P*rho_P;
				
				
			}
			
			if(i==(NI-2))
			{
				u_E=*(slabUmass[0]+(NJ-1)*k+j);
				v_E=*(slabVmass[0]+(NJ-1)*k+j);
				vol_E=*(slabVolume[0]+(NJ-1)*k+j);

				if ((base->flagParticle==1)|| (base->flagDroplet==1 && (base->BaseParameters->flagVariableDia==0)))
				{
					u_ED=*(slabUDisperse[0]+(NJ-1)*k+j);
					v_ED=*(slabVDisperse[0]+(NJ-1)*k+j);
				/*	if(!(*(slabVartheta[0]+(NJ-1)*k+j) <= EPSILON))
					{
						rho_ED = (*(slabVartheta[0]+(NJ-1)*k+j) <= EPSILON) ? 0.0 :  (1.0/(*(slabRhoDisperse[0]+((NJ-1)*k+j)*flagrhodisperse)*(*(slabVartheta[0]+(NJ-1)*k+j))));
						
						u_ED= u_ED*rho_ED;
						v_ED = v_ED*rho_ED;
					}
					else				
					{
						u_ED= u_E;
						v_ED = v_E;
					} */
				}
				if (base->flagDroplet==1 && (base->BaseParameters->flagVariableDia==1))
				{
					u_ED=*(slabUDmass[0]+(NJ-1)*k+j);
					v_ED=*(slabVDmass[0]+(NJ-1)*k+j);
				}
				if(base->flagMPFV==1)
				{
					rho_E = 1/(*(slabRho[0]+(NJ-1)*k+j));
					u_E = u_E*rho_E;
					v_E = v_E*rho_E;
					
					
				}
				
			}
			else
			{
				u_E=*(domainUmass+indexE);
				v_E=*(domainVmass+indexE);
				vol_E=*(domainVolume+indexE);
				if ((base->flagParticle==1)|| (base->flagDroplet==1 && (base->BaseParameters->flagVariableDia==0)))
				{
					u_ED=*(domainUDisperse+indexE);
					v_ED=*(domainVDisperse+indexE);
				/*	if(!(*(domainVartheta+indexE)<=EPSILON))
					{
						rho_ED = (*(domainVartheta+indexE)<=EPSILON)? 0.0 : (1.0/(*(domainRhoDisperse+indexE*flagrhodisperse)*(*(domainVartheta+indexE))));
						u_ED = u_ED*rho_ED;
						v_ED = v_ED*rho_ED;
					}
					else
					{
						u_ED = u_E;
						v_ED = v_E;
					}	*/		
				}
				if ((base->flagDroplet==1 && (base->BaseParameters->flagVariableDia==1)))
				{
					u_ED=*(domainUDmass+indexE);
					v_ED=*(domainVDmass+indexE);
				}
				if(base->flagMPFV==1)
				{
					rho_E = 1/(*(domainRho+indexE));
					u_E = u_E*rho_E;
					v_E = v_E*rho_E;

				}
			}
			
			if(i==0)
			{
				u_W =*(slabUmass[1]+(NJ-1)*k+j);
				v_W =*(slabVmass[1]+(NJ-1)*k+j);
				vol_W =*(slabVolume[1]+(NJ-1)*k+j);
				if ((base->flagParticle==1)|| (base->flagDroplet==1 && (base->BaseParameters->flagVariableDia==0)))
				{
					u_WD =*(slabUDisperse[1]+(NJ-1)*k+j);
					v_WD =*(slabVDisperse[1]+(NJ-1)*k+j);
				/*	if(!(*(slabVartheta[1]+(NJ-1)*k+j)<=EPSILON))
					{
						rho_WD = (*(slabVartheta[1]+(NJ-1)*k+j)<=EPSILON)? 0.0 : (1.0/(*(slabRhoDisperse[1]+((NJ-1)*k+j)*flagrhodisperse)*(*(slabVartheta[1]+(NJ-1)*k+j))));
						u_WD = u_WD*rho_WD;
						v_WD = v_WD*rho_WD;
					}
					else
					{
						u_WD = u_W;
						v_WD = v_W;
					}				
				*/	 					
				}
				if ((base->flagDroplet==1 && (base->BaseParameters->flagVariableDia==1)))
				{
					u_WD =*(slabUDmass[1]+(NJ-1)*k+j);
					v_WD =*(slabVDmass[1]+(NJ-1)*k+j);
				}
				if(base->flagMPFV==1)
				{
					rho_W = 1/(*(slabRho[1]+(NJ-1)*k+j));
					u_W = u_W*rho_W;
					v_W = v_W*rho_W;
	
				}
				
			}
			else
			{
				u_W=*(domainUmass+indexW);
				v_W=*(domainVmass+indexW);
				vol_W=*(domainVolume+indexW);
				if ((base->flagParticle==1)|| (base->flagDroplet==1 && (base->BaseParameters->flagVariableDia==0)))
				{
					u_WD=*(domainUDisperse+indexW);
					v_WD=*(domainVDisperse+indexW);
				/*	if(!(*(domainVartheta+indexW) <=EPSILON))
					{
						rho_WD = (*(domainVartheta+indexW) <=EPSILON)? 0.0 : (1.0/(*(domainRhoDisperse+indexW*flagrhodisperse)*(*(domainVartheta+indexW))));
						u_WD = u_WD*rho_WD;
						v_WD = v_WD*rho_WD;
					}
					else
					{
						u_WD = u_W;
						v_WD = v_W;
					}	
				*/
				}
				if ((base->flagDroplet==1 && (base->BaseParameters->flagVariableDia==1)))
				{
					u_WD=*(domainUDmass+indexW);
					v_WD=*(domainVDmass+indexW);
				}
				if(base->flagMPFV==1)
				{
					rho_W = 1/(*(domainRho+indexW));
					u_W = u_W*rho_W;
					v_W = v_W*rho_W;

				}
			}
			
			if(j==(NJ-2))
			{
				u_N=*(slabUmass[2]+(NI-1)*k+i);
				v_N=*(slabVmass[2]+(NI-1)*k+i);
				vol_N=*(slabVolume[2]+(NI-1)*k+i);
				if ((base->flagParticle==1)|| (base->flagDroplet==1 && (base->BaseParameters->flagVariableDia==0)))
				{
					u_ND=*(slabUDisperse[2]+(NI-1)*k+i);
					v_ND=*(slabVDisperse[2]+(NI-1)*k+i);
				/*	if(!(*(slabVartheta[2]+(NI-1)*k+i)<= EPSILON))
					{
						rho_ND = (*(slabVartheta[2]+(NI-1)*k+i)<= EPSILON) ? 0.0 : (1.0/(*(slabRhoDisperse[2]+((NI-1)*k+i)*flagrhodisperse)*(*(slabVartheta[2]+(NI-1)*k+i))));
						u_ND = u_ND*rho_ND;
						v_ND = v_ND*rho_ND;
					}
					else
					{
						u_ND = u_N;
						v_ND = v_N;
					}
				*/
	
				}
				if ((base->flagDroplet==1 && (base->BaseParameters->flagVariableDia==1)))
				{
					u_ND=*(slabUDmass[2]+(NI-1)*k+i);
					v_ND=*(slabVDmass[2]+(NI-1)*k+i);
				}
				if(base->flagMPFV==1)
				{
					rho_N = 1/(*(slabRho[2]+(NI-1)*k+i));
					u_N = u_N*rho_N;
					v_N = v_N*rho_N;
					
				}
			}
			else
			{
				u_N=*(domainUmass+indexN);
				v_N=*(domainVmass+indexN);
				vol_N=*(domainVolume+indexN);
				if ((base->flagParticle==1)|| (base->flagDroplet==1 && (base->BaseParameters->flagVariableDia==0)))
				{
					u_ND=*(domainUDisperse+indexN);
					v_ND=*(domainVDisperse+indexN);
				/*	if(!(*(domainVartheta+indexN)<=EPSILON))
					{
						rho_ND = (*(domainVartheta+indexN)<=EPSILON) ? 0.0 :(1.0/(*(domainRhoDisperse+indexN*flagrhodisperse)*(*(domainVartheta+indexN))));
						u_ND = u_ND*rho_ND;
						v_ND = v_ND*rho_ND;
					}
					else
					{
						u_ND = u_N;
						v_ND = v_N;
					}
				*/				
				}
				if ((base->flagDroplet==1 && (base->BaseParameters->flagVariableDia==1)))
				{
					u_ND=*(domainUDmass+indexN);
					v_ND=*(domainVDmass+indexN);
				}
				if(base->flagMPFV==1)
				{
					rho_N = 1/(*(domainRho+indexN));
					u_N = u_N*rho_N;
					v_N = v_N*rho_N;
					
				}
			}
			
			if(j==0)
			{
				u_S=*(slabUmass[3]+(NI-1)*k+i);
				v_S=*(slabVmass[3]+(NI-1)*k+i);
				vol_S=*(slabVolume[3]+(NI-1)*k+i);
				if ((base->flagParticle==1)|| (base->flagDroplet==1 && (base->BaseParameters->flagVariableDia==0)))
				{
					u_SD=*(slabUDisperse[3]+(NI-1)*k+i);
					v_SD=*(slabVDisperse[3]+(NI-1)*k+i);
				/*	if(!(*(slabVartheta[3]+(NI-1)*k+i)<=EPSILON))
					{
						rho_SD = (*(slabVartheta[3]+(NI-1)*k+i)<=EPSILON) ? 0.0:  (1.0/(*(slabRhoDisperse[3]+((NI-1)*k+i)*flagrhodisperse)*(*(slabVartheta[3]+(NI-1)*k+i))));
					
						u_SD = u_SD*rho_SD;
						v_SD = v_SD*rho_SD;
					}
					else
					{
						u_SD = u_S;
						v_SD = v_S;
					}
				*/				
				}
				if((base->flagDroplet==1 && (base->BaseParameters->flagVariableDia==1)))
				{
					u_SD=*(slabUDmass[3]+(NI-1)*k+i);
					v_SD=*(slabVDmass[3]+(NI-1)*k+i);
				}
				if(base->flagMPFV==1)
				{
					rho_S = 1/(*(slabRho[3]+(NI-1)*k+i));
					u_S = u_S*rho_S;
					v_S = v_S*rho_S;
					
				}
			}
			else
			{
				u_S=*(domainUmass+indexS);
				v_S=*(domainVmass+indexS);
				vol_S=*(domainVolume+indexS);
				if ((base->flagParticle==1)|| (base->flagDroplet==1 && (base->BaseParameters->flagVariableDia==0)))
				{
					u_SD=*(domainUDisperse+indexS);
					v_SD=*(domainVDisperse+indexS);
				/*	if(!(*(domainVartheta+indexS)<=EPSILON) )
					{
						rho_SD = (*(domainVartheta+indexS)<=EPSILON) ? 0.0 : (1.0/(*(domainRhoDisperse+indexS*flagrhodisperse)*(*(domainVartheta+indexS))));
					
						u_SD = u_SD*rho_SD;
						v_SD = v_SD*rho_SD;
					}
					else
					{
						u_SD = u_S;
						v_SD = v_S;
					}
				*/
				}
				if((base->flagDroplet==1 && (base->BaseParameters->flagVariableDia==1)))
				{
					u_SD=*(domainUDmass+indexS);
					v_SD=*(domainVDmass+indexS);
				}
				if(base->flagMPFV==1)
				{
					rho_S = 1/(*(domainRho+indexS));
					u_S = u_S*rho_S;
					v_S = v_S*rho_S;
					
				}
			}
			
			ue=((u_E*vol_P+u_P*vol_E)/(vol_P+vol_E));
			uw=((u_W*vol_P+u_P*vol_W)/(vol_P+vol_W));
			us=((u_S*vol_P+u_P*vol_S)/(vol_P+vol_S));
			un=((u_N*vol_P+u_P*vol_N)/(vol_P+vol_N));
			
			ve=((v_E*vol_P+v_P*vol_E)/(vol_P+vol_E));
			vw=((v_W*vol_P+v_P*vol_W)/(vol_P+vol_W));
			vs=((v_S*vol_P+v_P*vol_S)/(vol_P+vol_S));
			vn=((v_N*vol_P+v_P*vol_N)/(vol_P+vol_N));
			
			Index=((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6*3;
			face = 0;  // west face
				*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = (uw*(*(SurfArea+Index+3*0+0))+vw*(*(SurfArea+Index+3*0+1)));
			face = 1;  // east face
				*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = (ue*(*(SurfArea+Index+3*1+0))+ve*(*(SurfArea+Index+3*1+1)));
			face = 2;  // south face
				*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = (us*(*(SurfArea+Index+3*2+0))+vs*(*(SurfArea+Index+3*2+1)));
			face = 3;  // north face
				*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = (un*(*(SurfArea+Index+3*3+0))+vn*(*(SurfArea+Index+3*3+1)));
			face = 4;  // bottom face
				*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = 0.0;
			face = 5;  // top face
				*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = 0.0;

		//	fprintf(fp,"%e\t%e\t%e\t%e\t%e\t%e\n",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+0),*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+1),*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+2),*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+3),*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+4),*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+5));
			
			//####################added by adhiraj for disperse phase flow#################################################
			
				if ((base->flagParticle==1)|| (base->flagDroplet==1))
				{
					ueD=((u_ED*vol_P+u_PD*vol_E)/(vol_P+vol_E));
					uwD=((u_WD*vol_P+u_PD*vol_W)/(vol_P+vol_W));
					usD=((u_SD*vol_P+u_PD*vol_S)/(vol_P+vol_S));
					unD=((u_ND*vol_P+u_PD*vol_N)/(vol_P+vol_N));
					
					veD=((v_ED*vol_P+v_PD*vol_E)/(vol_P+vol_E));
					vwD=((v_WD*vol_P+v_PD*vol_W)/(vol_P+vol_W));
					vsD=((v_SD*vol_P+v_PD*vol_S)/(vol_P+vol_S));
					vnD=((v_ND*vol_P+v_PD*vol_N)/(vol_P+vol_N));
					
			
					Index=((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6*3;
					face = 0;  // west face
					*(FD+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = (uwD*(*(SurfArea+Index+3*0+0))+vwD*(*(SurfArea+Index+3*0+1)));
					face = 1;  // east face
					*(FD+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = (ueD*(*(SurfArea+Index+3*1+0))+veD*(*(SurfArea+Index+3*1+1)));
					face = 2;  // south face
					*(FD+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = (usD*(*(SurfArea+Index+3*2+0))+vsD*(*(SurfArea+Index+3*2+1)));
					face = 3;  // north face
					*(FD+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = (unD*(*(SurfArea+Index+3*3+0))+vnD*(*(SurfArea+Index+3*3+1)));
					face = 4;  // bottom face
					*(FD+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = 0.0;
					face = 5;  // top face
					*(FD+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = 0.0;

							
					//fprintf(fp,"%e\t%e\t%e\t%e\t%e\t%e\n",*(FD+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+0),*(FD+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+1),*(FD+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+2),*(FD+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+3),*(FD+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+4),*(FD+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+5));
				}

				
//####################end of flux calculation for disperse phase###############################################
		}/*End of looping over NI,NJ,NK*/
	}/*End of looping over the zones */
	//fclose(fp);

	return 0;
}
