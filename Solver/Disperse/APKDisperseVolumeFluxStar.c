#include "APKDisperseVolumeFluxStar.h"
/**
This function is used for calculating the volumeflux from the current converged
values of massvelocities of Dispersed phase. This function is used for 3D domains.
*/
int DisperseVolumeFluxStar(Base *base)
{
	int nZones,numZones;
	int k,j,i;
	int indexB,indexT,indexE,indexW,indexS,indexN,indexP;

	double *domainVolume,*domainUDmass,*domainVDmass,*domainWDmass;
	double *slabVolume[6],*slabUDmass[6],*slabVDmass[6],*slabWDmass[6],*slabRho[6],*domainRho;	
	double *SurfArea,*F;
	Zone* zone;
	nZones=(base->Number_of_zones);
	
	int NI,NJ,NK,*ptrSize;
	int K_Increment,J_Increment,I_Increment,face_index,area_index;
	
	double uw,vw,ww,ue,ve,we,un,vn,wn,us,vs,ws,ut,vt,wt,ub,vb,wb;
	double u_P,u_E,u_W,u_S,u_N,u_T,u_B,v_P,v_E,v_W,v_S,v_N,v_T,v_B,w_P,w_E,w_W,w_S,w_N,w_T,w_B;
	double vol_P,vol_E,vol_W,vol_N,vol_S,vol_T,vol_B,rho_P,rho_E,rho_W,rho_S,rho_N,rho_T,rho_B;
	
	double temp;

	for(numZones=0;numZones<nZones;numZones++)
	{
		zone=(base->CgnsZone+numZones);
		
		ptrSize=zone->Size;

		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);

		F = zone->Flowsol->MassfluxDisperse;

		SurfArea=zone->ZGeom->SurfaceArea;
		
		domainUDmass=MakeDomainArray(zone,6200);
		domainVDmass=MakeDomainArray(zone,6300);
		domainWDmass=MakeDomainArray(zone,6400);
		domainVolume=MakeDomainArray(zone,5);
	
		MakeSlabPtrArray(zone,slabUDmass,6200);
		MakeSlabPtrArray(zone,slabVDmass,6300);
		MakeSlabPtrArray(zone,slabWDmass,6400);
		MakeSlabPtrArray(zone,slabVolume,5);
		indexP = 0;
		face_index=0;
		area_index=0;
		
		K_Increment = (NI-1)*(NJ-1);
		J_Increment = (NI-1);
		I_Increment = 1;
		
		for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{
					
			indexE=indexP+1;
			indexW=indexP-1;
			indexN=indexP+J_Increment;
			indexS=indexP-J_Increment;
			indexT=indexP+K_Increment;
			indexB=indexP-K_Increment;

			u_P=*(domainUDmass+indexP);
			v_P=*(domainVDmass+indexP);
			w_P=*(domainWDmass+indexP);
			vol_P=*(domainVolume+indexP);
				
			if(i==(NI-2))
			{
				u_E=*(slabUDmass[0]+(NJ-1)*k+j);
				v_E=*(slabVDmass[0]+(NJ-1)*k+j);
				w_E=*(slabWDmass[0]+(NJ-1)*k+j);
				vol_E=*(slabVolume[0]+(NJ-1)*k+j);
			}
			else
			{
				u_E=*(domainUDmass+indexE);
				v_E=*(domainVDmass+indexE);
				w_E=*(domainWDmass+indexE);
				vol_E=*(domainVolume+indexE);	
			}
			if(i==0)
			{
				u_W =*(slabUDmass[1]+(NJ-1)*k+j);
				v_W =*(slabVDmass[1]+(NJ-1)*k+j);
				w_W =*(slabWDmass[1]+(NJ-1)*k+j);
				vol_W =*(slabVolume[1]+(NJ-1)*k+j);
			}
			else
			{
				u_W=*(domainUDmass+indexW);
				v_W=*(domainVDmass+indexW);
				w_W=*(domainWDmass+indexW);
				vol_W=*(domainVolume+indexW);
			}
			if(j==(NJ-2))
			{
				u_N=*(slabUDmass[2]+(NI-1)*k+i);
				v_N=*(slabVDmass[2]+(NI-1)*k+i);
				w_N=*(slabWDmass[2]+(NI-1)*k+i);
				vol_N=*(slabVolume[2]+(NI-1)*k+i);
			}
			else
			{
				u_N=*(domainUDmass+indexN);
				v_N=*(domainVDmass+indexN);
				w_N=*(domainWDmass+indexN);
				vol_N=*(domainVolume+indexN);
			}
			if(j==0)
			{
				u_S=*(slabUDmass[3]+(NI-1)*k+i);
				v_S=*(slabVDmass[3]+(NI-1)*k+i);
				w_S=*(slabWDmass[3]+(NI-1)*k+i);
				vol_S=*(slabVolume[3]+(NI-1)*k+i);
			}
			else
			{
				u_S=*(domainUDmass+indexS);
				v_S=*(domainVDmass+indexS);
				w_S=*(domainWDmass+indexS);
				vol_S=*(domainVolume+indexS);
			}
			if(k==(NK-2))
			{
				u_T=*(slabUDmass[4]+(NI-1)*j+i);
				v_T=*(slabVDmass[4]+(NI-1)*j+i);
				w_T=*(slabWDmass[4]+(NI-1)*j+i);
				vol_T=*(slabVolume[4]+(NI-1)*j+i);
			}
			else
			{
				u_T=*(domainUDmass+indexT);
				v_T=*(domainVDmass+indexT);
				w_T=*(domainWDmass+indexT);
				vol_T=*(domainVolume+indexT);
			}
			if(k==0)
			{
				u_B=*(slabUDmass[5]+(NI-1)*j+i);
				v_B=*(slabVDmass[5]+(NI-1)*j+i);
				w_B=*(slabWDmass[5]+(NI-1)*j+i);
				vol_B=*(slabVolume[5]+(NI-1)*j+i);
			}
			else
			{
				u_B=*(domainUDmass+indexB);
				v_B=*(domainVDmass+indexB);
				w_B=*(domainWDmass+indexB);
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

			*(F+face_index) = (uw*(*(SurfArea+area_index+0))+vw*(*(SurfArea+area_index+1))+ww*(*(SurfArea+area_index+2)));
			*(F+face_index+1) = (ue*(*(SurfArea+area_index+3))+ve*(*(SurfArea+area_index+4))+we*(*(SurfArea+area_index+5)));
			*(F+face_index+2) = (us*(*(SurfArea+area_index+6))+vs*(*(SurfArea+area_index+7))+ws*(*(SurfArea+area_index+8)));
			*(F+face_index+3) = (un*(*(SurfArea+area_index+9))+vn*(*(SurfArea+area_index+10))+wn*(*(SurfArea+area_index+11)));
			*(F+face_index+4) = (ub*(*(SurfArea+area_index+12))+vb*(*(SurfArea+area_index+13))+wb*(*(SurfArea+area_index+14)));
			*(F+face_index+5) = (ut*(*(SurfArea+area_index+15))+vt*(*(SurfArea+area_index+16))+wt*(*(SurfArea+area_index+17)));
			
	/*		if(i==(int)NI/2 && j==(int)NJ/2)
			{
				temp = *(F+face_index)+*(F+face_index+1)+*(F+face_index+2)+*(F+face_index+3)+*(F+face_index+4)+*(F+face_index+5);
				printf("total volume flux of dispersed phase at %d,%d= %e",i,j,temp);
			}
		*/
			indexP++;
			face_index = face_index + 6;
			area_index = area_index + 18;
			
		}	/*End of looping over NI,NJ,NK*/
		
	}/*End of looping over the zones */
	return 0;
}

/**
This function is used for calculating the volumeflux from the current converged
values of massvelocities of Disperse Phase. This function is used for 2D domains.
*/
int DisperseVolumeFluxStar2d(Base *base)
{
	int nZones,numZones;
	int k,j,i;
	int indexE,indexW,indexS,indexN,indexP;

	double *domainVolume,*domainUDmass,*domainVDmass;
	double *slabVolume[6],*slabUDmass[6],*slabVDmass[6],*slabRho[6],*domainRho;	
	double *SurfArea,*F;
	Zone* zone;
	nZones=(base->Number_of_zones);
	
	int NI,NJ,NK,*ptrSize;
	int K_Increment,J_Increment,I_Increment,face_index,area_index;
	
	double uw,vw,ue,ve,un,vn,us,vs;
	double u_P,u_E,u_W,u_S,u_N,v_P,v_E,v_W,v_S,v_N;
	double vol_P,vol_E,vol_W,vol_N,vol_S,rho_P,rho_E,rho_W,rho_S,rho_N;
	
	double temp;

	for(numZones=0;numZones<nZones;numZones++)
	{
		zone=(base->CgnsZone+numZones);
		
		ptrSize=zone->Size;

		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);

		F = zone->Flowsol->MassfluxDisperse;

		SurfArea=zone->ZGeom->SurfaceArea;
		
		domainUDmass=MakeDomainArray(zone,6200);
		domainVDmass=MakeDomainArray(zone,6300);
		domainVolume=MakeDomainArray(zone,5);
	
		MakeSlabPtrArray(zone,slabUDmass,6200);
		MakeSlabPtrArray(zone,slabVDmass,6300);
		MakeSlabPtrArray(zone,slabVolume,5);
		indexP = 0;
		face_index=0;
		area_index=0;
		
		K_Increment = (NI-1)*(NJ-1);
		J_Increment = (NI-1);
		I_Increment = 1;
		
		for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{
					
			indexE=indexP+1;
			indexW=indexP-1;
			indexN=indexP+J_Increment;
			indexS=indexP-J_Increment;
			
			u_P=*(domainUDmass+indexP);
			v_P=*(domainVDmass+indexP);
			vol_P=*(domainVolume+indexP);
				
			if(i==(NI-2))
			{
				u_E=*(slabUDmass[0]+(NJ-1)*k+j);
				v_E=*(slabVDmass[0]+(NJ-1)*k+j);
				vol_E=*(slabVolume[0]+(NJ-1)*k+j);
			}
			else
			{
				u_E=*(domainUDmass+indexE);
				v_E=*(domainVDmass+indexE);
				vol_E=*(domainVolume+indexE);	
			}
			if(i==0)
			{
				u_W =*(slabUDmass[1]+(NJ-1)*k+j);
				v_W =*(slabVDmass[1]+(NJ-1)*k+j);
				vol_W =*(slabVolume[1]+(NJ-1)*k+j);
			}
			else
			{
				u_W=*(domainUDmass+indexW);
				v_W=*(domainVDmass+indexW);
				vol_W=*(domainVolume+indexW);
			}
			if(j==(NJ-2))
			{
				u_N=*(slabUDmass[2]+(NI-1)*k+i);
				v_N=*(slabVDmass[2]+(NI-1)*k+i);
				vol_N=*(slabVolume[2]+(NI-1)*k+i);
			}
			else
			{
				u_N=*(domainUDmass+indexN);
				v_N=*(domainVDmass+indexN);
				vol_N=*(domainVolume+indexN);
			}
			if(j==0)
			{
				u_S=*(slabUDmass[3]+(NI-1)*k+i);
				v_S=*(slabVDmass[3]+(NI-1)*k+i);
				vol_S=*(slabVolume[3]+(NI-1)*k+i);
			}
			else
			{
				u_S=*(domainUDmass+indexS);
				v_S=*(domainVDmass+indexS);
				vol_S=*(domainVolume+indexS);
			}
			
			ue=((u_E*vol_P+u_P*vol_E)/(vol_P+vol_E));
			uw=((u_W*vol_P+u_P*vol_W)/(vol_P+vol_W));
			us=((u_S*vol_P+u_P*vol_S)/(vol_P+vol_S));
			un=((u_N*vol_P+u_P*vol_N)/(vol_P+vol_N));
			
			ve=((v_E*vol_P+v_P*vol_E)/(vol_P+vol_E));
			vw=((v_W*vol_P+v_P*vol_W)/(vol_P+vol_W));
			vs=((v_S*vol_P+v_P*vol_S)/(vol_P+vol_S));
			vn=((v_N*vol_P+v_P*vol_N)/(vol_P+vol_N));
			
			*(F+face_index) = uw*(*(SurfArea+area_index+0))+vw*(*(SurfArea+area_index+1));
			*(F+face_index+1) = ue*(*(SurfArea+area_index+3))+ve*(*(SurfArea+area_index+4));
			*(F+face_index+2) = us*(*(SurfArea+area_index+6))+vs*(*(SurfArea+area_index+7));
			*(F+face_index+3) = un*(*(SurfArea+area_index+9))+vn*(*(SurfArea+area_index+10));
			*(F+face_index+4) = 0.0;
			*(F+face_index+5) = 0.0;

/*			if(i==(int)NI/2 && j==(int)NJ/2)
			{
				printf("\nFor%d,%d:\nue= %f\tuw = %f\tus=%f\tun=%f\nve=%f\tvw=%f\tvs=%f\tvn=%f\n",i,j,ue,uw,us,un,ve,vw,vs,vn);
				temp = *(F+face_index)+*(F+face_index+1)+*(F+face_index+2)+*(F+face_index+3)+*(F+face_index+4)+*(F+face_index+5);
				printf("total volume flux of dispersed phase at %d,%d= %e\n\n",i,j,temp);
			}
*/
			indexP++;
			face_index = face_index + 6;
			area_index = area_index + 18;
			
		}	/*End of looping over NI,NJ,NK*/
		
	}/*End of looping over the zones */
	return 0;
}

