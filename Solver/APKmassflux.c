#include "APKmassflux.h"
/**
This function is used only once for calculating the massflux just before the first timestep using
the guess values of mass velocities. This is used for 3D geometries.
*/

int massflux(Base *base)
{
	int nZones,numZones;
	int k,j,i,index,indexB,indexT,indexE,indexW,indexS,indexN,indexP;
	int numiface,numBc,face,iface,zBc,dir,Index;
	double *SurfArea,*F,*rho;
	double *domainVolume,*domainUmass,*domainVmass,*domainWmass;
	double *slabVolume[6],*slabUmass[6],*slabVmass[6],*slabWmass[6];
	Zone* zone;
	double *slabRho[6],*domainRho;

	int NI,NJ,NK,*ptrSize;
	double uw,vw,ww,ue,ve,we,un,vn,wn,us,vs,ws,ut,vt,wt,ub,vb,wb;
	double u_P,u_E,u_W,u_S,u_N,u_T,u_B,v_P,v_E,v_W,v_S,v_N,v_T,v_B,w_P,w_E,w_W,w_S,w_N,w_T,w_B;
	double vol_P,vol_E,vol_W,vol_N,vol_S,vol_T,vol_B;
	double rho_P,rho_E,rho_W,rho_N,rho_S,rho_T,rho_B;
	double rhoue,rhove,rhowe,rhouw,rhovw,rhoww,rhoun,rhovn,rhown,rhous,rhovs,rhows,rhout,rhovt,rhowt,rhoub,rhovb,rhowb;
	
	nZones=base->Number_of_zones;	
	for(numZones=0;numZones<nZones;numZones++)
	{
		/*Copying down u,v,w into umass,vmass,wmass ....domain as well as slabs*/
		ptrSize=((base)->CgnsZone+numZones)->Size;
		zone=(base->CgnsZone+numZones);
		if((zone->flagSolidFluid==0))
		{
			NI=*(ptrSize+0);
			NJ=*(ptrSize+1);
			NK=*(ptrSize+2);

			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				indexP=(NJ-1)*(NI-1)*k+(NI-1)*j+i;
				*((zone->Flowsol)->umass+indexP)=*((zone->Flowsol)->u+indexP);
				*((zone->Flowsol)->vmass+indexP)=*((zone->Flowsol)->v+indexP);
				*((zone->Flowsol)->wmass+indexP)=*((zone->Flowsol)->w+indexP);
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
							*((((zone)->Iface+iface)->ZSlab)->umass+index) = *((((zone)->Iface+iface)->ZSlab)->u+index);
							*((((zone)->Iface+iface)->ZSlab)->vmass+index) = *((((zone)->Iface+iface)->ZSlab)->v+index);
							*((((zone)->Iface+iface)->ZSlab)->wmass+index) = *((((zone)->Iface+iface)->ZSlab)->w+index);
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
							*((((zone)->ZoneBc+zBc)->ZSlab)->umass+index) = *((((zone)->ZoneBc+zBc)->ZSlab)->u+index);
							*((((zone)->ZoneBc+zBc)->ZSlab)->vmass+index) = *((((zone)->ZoneBc+zBc)->ZSlab)->v+index);
							*((((zone)->ZoneBc+zBc)->ZSlab)->wmass+index) = *((((zone)->ZoneBc+zBc)->ZSlab)->w+index);
						}
					}
				}
			}/*End of direction Loop*/
		}
	}/*End of Zone loop*/

	nZones=base->Number_of_zones;
	for(numZones=0;numZones<nZones;numZones++)
	{
		zone=(base->CgnsZone+numZones);
			
		if( (zone->flagSolidFluid==0) && ((zone->Flowsol)->flagRho==1) )
		{/*Done By Ashish for Conjugate Heat Transfer Case*/	
			ptrSize=zone->Size;
			NI=*(ptrSize+0);
			NJ=*(ptrSize+1);
			NK=*(ptrSize+2);

			F = zone->Flowsol->Massflux;

			SurfArea=zone->ZGeom->SurfaceArea;

			rho=zone->Flowsol->Rho;

			domainUmass=MakeDomainArray(zone,100);
			domainVmass=MakeDomainArray(zone,200);
			domainWmass=MakeDomainArray(zone,300);
			domainVolume=MakeDomainArray(zone,5);
			domainRho=MakeDomainArray(zone,710);/*Density*/

			MakeSlabPtrArray(zone,slabUmass,100);
			MakeSlabPtrArray(zone,slabVmass,200);
			MakeSlabPtrArray(zone,slabWmass,300);
			MakeSlabPtrArray(zone,slabVolume,5);
			MakeSlabPtrArray(zone,slabRho,710);

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
				rho_P=*(domainRho+indexP);
				vol_P=*(domainVolume+indexP);

				if(i==(NI-2))
				{
					u_E=*(slabUmass[0]+(NJ-1)*k+j);
					v_E=*(slabVmass[0]+(NJ-1)*k+j);
					w_E=*(slabWmass[0]+(NJ-1)*k+j);
					rho_E=*(slabRho[0]+(NJ-1)*k+j);
					vol_E=*(slabVolume[0]+(NJ-1)*k+j);
				}
				else
				{
					u_E=*(domainUmass+indexE);
					v_E=*(domainVmass+indexE);
					w_E=*(domainWmass+indexE);
					rho_E=*(domainRho+indexE);
					vol_E=*(domainVolume+indexE);
				}

				if(i==0)
				{
					u_W =*(slabUmass[1]+(NJ-1)*k+j);
					v_W =*(slabVmass[1]+(NJ-1)*k+j);
					w_W =*(slabWmass[1]+(NJ-1)*k+j);
					rho_W =*(slabRho[1]+(NJ-1)*k+j);
					vol_W =*(slabVolume[1]+(NJ-1)*k+j);
				}
				else
				{
					u_W=*(domainUmass+indexW);
					v_W=*(domainVmass+indexW);
					w_W=*(domainWmass+indexW);
					rho_W=*(domainRho+indexW);
					vol_W=*(domainVolume+indexW);
				}

				if(j==(NJ-2))
				{
					u_N=*(slabUmass[2]+(NI-1)*k+i);
					v_N=*(slabVmass[2]+(NI-1)*k+i);
					w_N=*(slabWmass[2]+(NI-1)*k+i);
					rho_N=*(slabRho[2]+(NI-1)*k+i);
					vol_N=*(slabVolume[2]+(NI-1)*k+i);
				}
				else
				{
					u_N=*(domainUmass+indexN);
					v_N=*(domainVmass+indexN);
					w_N=*(domainWmass+indexN);
					rho_N=*(domainRho+indexN);
					vol_N=*(domainVolume+indexN);
				}

				if(j==0)
				{
					u_S=*(slabUmass[3]+(NI-1)*k+i);
					v_S=*(slabVmass[3]+(NI-1)*k+i);
					w_S=*(slabWmass[3]+(NI-1)*k+i);
					rho_S=*(slabRho[3]+(NI-1)*k+i);
					vol_S=*(slabVolume[3]+(NI-1)*k+i);
				}
				else
				{
					u_S=*(domainUmass+indexS);
					v_S=*(domainVmass+indexS);
					w_S=*(domainWmass+indexS);
					rho_S=*(domainRho+indexS);
					vol_S=*(domainVolume+indexS);
				}

				if(k==(NK-2))
				{
					u_T=*(slabUmass[4]+(NI-1)*j+i);
					v_T=*(slabVmass[4]+(NI-1)*j+i);
					w_T=*(slabWmass[4]+(NI-1)*j+i);
					rho_T=*(slabRho[4]+(NI-1)*j+i);
					vol_T=*(slabVolume[4]+(NI-1)*j+i);
				}
				else
				{
					u_T=*(domainUmass+indexT);
					v_T=*(domainVmass+indexT);
					w_T=*(domainWmass+indexT);
					rho_T=*(domainRho+indexT);
					vol_T=*(domainVolume+indexT);
				}

				if(k==0)
				{
					u_B=*(slabUmass[5]+(NI-1)*j+i);
					v_B=*(slabVmass[5]+(NI-1)*j+i);
					w_B=*(slabWmass[5]+(NI-1)*j+i);
					rho_B=*(slabRho[5]+(NI-1)*j+i);
					vol_B=*(slabVolume[5]+(NI-1)*j+i);
				}
				else
				{
					u_B=*(domainUmass+indexB);
					v_B=*(domainVmass+indexB);
					w_B=*(domainWmass+indexB);
					rho_B=*(domainRho+indexB);
					vol_B=*(domainVolume+indexB);
				}


				rhoue=((rho_E*u_E*vol_P+rho_P*u_P*vol_E)/(vol_P+vol_E));
				rhove=((rho_E*v_E*vol_P+rho_P*v_P*vol_E)/(vol_P+vol_E));
				rhowe=((rho_E*w_E*vol_P+rho_P*w_P*vol_E)/(vol_P+vol_E));

				rhouw=((rho_W*u_W*vol_P+rho_P*u_P*vol_W)/(vol_P+vol_W));
				rhovw=((rho_W*v_W*vol_P+rho_P*v_P*vol_W)/(vol_P+vol_W));
				rhoww=((rho_W*w_W*vol_P+rho_P*w_P*vol_W)/(vol_P+vol_W));

				rhoun=((rho_N*u_N*vol_P+rho_P*u_P*vol_N)/(vol_P+vol_N));
				rhovn=((rho_N*v_N*vol_P+rho_P*v_P*vol_N)/(vol_P+vol_N));
				rhown=((rho_N*w_N*vol_P+rho_P*w_P*vol_N)/(vol_P+vol_N));

				rhous=((rho_S*u_S*vol_P+rho_P*u_P*vol_S)/(vol_P+vol_S));
				rhovs=((rho_S*v_S*vol_P+rho_P*v_P*vol_S)/(vol_P+vol_S));
				rhows=((rho_S*w_S*vol_P+rho_P*w_P*vol_S)/(vol_P+vol_S));

				rhout=((rho_T*u_T*vol_P+rho_P*u_P*vol_T)/(vol_P+vol_T));
				rhovt=((rho_T*v_T*vol_P+rho_P*v_P*vol_T)/(vol_P+vol_T));
				rhowt=((rho_T*w_T*vol_P+rho_P*w_P*vol_T)/(vol_P+vol_T));

				rhoub=((rho_B*u_B*vol_P+rho_P*u_P*vol_B)/(vol_P+vol_B));
				rhovb=((rho_B*v_B*vol_P+rho_P*v_P*vol_B)/(vol_P+vol_B));
				rhowb=((rho_B*w_B*vol_P+rho_P*w_P*vol_B)/(vol_P+vol_B));

				Index=((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6*3;
				face = 0;  // west face
					*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = (rhouw*(*(SurfArea+Index+3*0+0))+rhovw*(*(SurfArea+Index+3*0+1))+rhoww*(*(SurfArea+Index+3*0+2)));
					//fprintf(fp1,"west=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
				face = 1;  // east face
					*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = (rhoue*(*(SurfArea+Index+3*1+0))+rhove*(*(SurfArea+Index+3*1+1))+rhowe*(*(SurfArea+Index+3*1+2)));
					//fprintf(fp1,"east=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
				face = 2;  // south face
					*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = (rhous*(*(SurfArea+Index+3*2+0))+rhovs*(*(SurfArea+Index+3*2+1))+rhows*(*(SurfArea+Index+3*2+2)));
					//fprintf(fp1,"south=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
				face = 3;  // north face
					*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = (rhoun*(*(SurfArea+Index+3*3+0))+rhovn*(*(SurfArea+Index+3*3+1))+rhown*(*(SurfArea+Index+3*3+2)));
					//fprintf(fp1,"north=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
				face = 4;  // bottom face
					*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = (rhoub*(*(SurfArea+Index+3*4+0))+rhovb*(*(SurfArea+Index+3*4+1))+rhowb*(*(SurfArea+Index+3*4+2)));
					//fprintf(fp1,"bottom=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
				face = 5;  // top face
					*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = (rhout*(*(SurfArea+Index+3*5+0))+rhovt*(*(SurfArea+Index+3*5+1))+rhowt*(*(SurfArea+Index+3*5+2)));
					//fprintf(fp1,"top=%lf\n",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
			}/*End of looping over NI,NJ,NK*/
		}/*end of if flagRho==1*/
		if( (zone->flagSolidFluid==0) && ((zone->Flowsol)->flagRho==0) )
		{/*Done By Ashish for Conjugate Heat Transfer Case*/	
			ptrSize=zone->Size;

			NI=*(ptrSize+0);
			NJ=*(ptrSize+1);
			NK=*(ptrSize+2);

			F = zone->Flowsol->Massflux;

			SurfArea=zone->ZGeom->SurfaceArea;

			rho=zone->Flowsol->Rho;

			domainUmass=MakeDomainArray(zone,100);
			domainVmass=MakeDomainArray(zone,200);
			domainWmass=MakeDomainArray(zone,300);
			domainVolume=MakeDomainArray(zone,5);
			
			MakeSlabPtrArray(zone,slabUmass,100);
			MakeSlabPtrArray(zone,slabVmass,200);
			MakeSlabPtrArray(zone,slabWmass,300);
			MakeSlabPtrArray(zone,slabVolume,5);
			
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

				if(i==(NI-2))
				{
					u_E=*(slabUmass[0]+(NJ-1)*k+j);
					v_E=*(slabVmass[0]+(NJ-1)*k+j);
					w_E=*(slabWmass[0]+(NJ-1)*k+j);
					vol_E=*(slabVolume[0]+(NJ-1)*k+j);
				}
				else
				{
					u_E=*(domainUmass+indexE);
					v_E=*(domainVmass+indexE);
					w_E=*(domainWmass+indexE);
					vol_E=*(domainVolume+indexE);
				}

				if(i==0)
				{
					u_W =*(slabUmass[1]+(NJ-1)*k+j);
					v_W =*(slabVmass[1]+(NJ-1)*k+j);
					w_W =*(slabWmass[1]+(NJ-1)*k+j);
					vol_W =*(slabVolume[1]+(NJ-1)*k+j);
				}
				else
				{
					u_W=*(domainUmass+indexW);
					v_W=*(domainVmass+indexW);
					w_W=*(domainWmass+indexW);
					vol_W=*(domainVolume+indexW);
				}

				if(j==(NJ-2))
				{
					u_N=*(slabUmass[2]+(NI-1)*k+i);
					v_N=*(slabVmass[2]+(NI-1)*k+i);
					w_N=*(slabWmass[2]+(NI-1)*k+i);
					vol_N=*(slabVolume[2]+(NI-1)*k+i);
				}
				else
				{
					u_N=*(domainUmass+indexN);
					v_N=*(domainVmass+indexN);
					w_N=*(domainWmass+indexN);
					vol_N=*(domainVolume+indexN);
				}

				if(j==0)
				{
					u_S=*(slabUmass[3]+(NI-1)*k+i);
					v_S=*(slabVmass[3]+(NI-1)*k+i);
					w_S=*(slabWmass[3]+(NI-1)*k+i);
					vol_S=*(slabVolume[3]+(NI-1)*k+i);
				}
				else
				{
					u_S=*(domainUmass+indexS);
					v_S=*(domainVmass+indexS);
					w_S=*(domainWmass+indexS);
					vol_S=*(domainVolume+indexS);
				}

				if(k==(NK-2))
				{
					u_T=*(slabUmass[4]+(NI-1)*j+i);
					v_T=*(slabVmass[4]+(NI-1)*j+i);
					w_T=*(slabWmass[4]+(NI-1)*j+i);
					vol_T=*(slabVolume[4]+(NI-1)*j+i);
				}
				else
				{
					u_T=*(domainUmass+indexT);
					v_T=*(domainVmass+indexT);
					w_T=*(domainWmass+indexT);
					vol_T=*(domainVolume+indexT);
				}

				if(k==0)
				{
					u_B=*(slabUmass[5]+(NI-1)*j+i);
					v_B=*(slabVmass[5]+(NI-1)*j+i);
					w_B=*(slabWmass[5]+(NI-1)*j+i);
					vol_B=*(slabVolume[5]+(NI-1)*j+i);
				}
				else
				{
					u_B=*(domainUmass+indexB);
					v_B=*(domainVmass+indexB);
					w_B=*(domainWmass+indexB);
					vol_B=*(domainVolume+indexB);
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
					*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = *rho*(uw*(*(SurfArea+Index+3*0+0))+vw*(*(SurfArea+Index+3*0+1))+ww*(*(SurfArea+Index+3*0+2)));
					//fprintf(fp1,"west=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
				face = 1;  // east face
					*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = *rho*(ue*(*(SurfArea+Index+3*1+0))+ve*(*(SurfArea+Index+3*1+1))+we*(*(SurfArea+Index+3*1+2)));
					//fprintf(fp1,"east=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
				face = 2;  // south face
					*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = *rho*(us*(*(SurfArea+Index+3*2+0))+vs*(*(SurfArea+Index+3*2+1))+ws*(*(SurfArea+Index+3*2+2)));
					//fprintf(fp1,"south=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
				face = 3;  // north face
					*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = *rho*(un*(*(SurfArea+Index+3*3+0))+vn*(*(SurfArea+Index+3*3+1))+wn*(*(SurfArea+Index+3*3+2)));
					//fprintf(fp1,"north=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
				face = 4;  // bottom face
					*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = *rho*(ub*(*(SurfArea+Index+3*4+0))+vb*(*(SurfArea+Index+3*4+1))+wb*(*(SurfArea+Index+3*4+2)));
					//fprintf(fp1,"bottom=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
				face = 5;  // top face
					*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = *rho*(ut*(*(SurfArea+Index+3*5+0))+vt*(*(SurfArea+Index+3*5+1))+wt*(*(SurfArea+Index+3*5+2)));
					//fprintf(fp1,"top=%lf\n",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
			}/*End of looping over NI,NJ,NK*/
		}/*end of else */
	}/*End of looping over the zones */
	return 0;
}
	
		
