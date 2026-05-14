#include "APKvolumefluxstar.h"
/**
This function is used for calculating the volumeflux from the current converged
values of massvelocities. This function is used for 3D domains.
See Eqn No(3.36) Project Report No. 1, Dec 2005.
*/
int volumefluxstar(Base *base)
{
	int nZones,numZones;
	int k,j,i;
	int indexB,indexT,indexE,indexW,indexS,indexN,indexP;

	double *domainVolume,*domainUmass,*domainVmass,*domainWmass;
	double *slabVolume[6],*slabUmass[6],*slabVmass[6],*slabWmass[6],*slabRho[6],*domainRho;	
	double *SurfArea,*F;
	Zone* zone;
	nZones=(base->Number_of_zones);
	
	int NI,NJ,NK,*ptrSize;
	int K_Increment,J_Increment,I_Increment,face_index,area_index;
	
	double uw,vw,ww,ue,ve,we,un,vn,wn,us,vs,ws,ut,vt,wt,ub,vb,wb;
	double u_P,u_E,u_W,u_S,u_N,u_T,u_B,v_P,v_E,v_W,v_S,v_N,v_T,v_B,w_P,w_E,w_W,w_S,w_N,w_T,w_B;
	double vol_P,vol_E,vol_W,vol_N,vol_S,vol_T,vol_B,rho_P,rho_E,rho_W,rho_S,rho_N,rho_T,rho_B;
	
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
			
// 				indexE=(NJ-1)*(NI-1)*k+(NI-1)*j+(i+1);
// 				indexW=(NJ-1)*(NI-1)*k+(NI-1)*j+(i-1);
// 				indexN=(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i;
// 				indexS=(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i;
// 				indexT=(NJ-1)*(NI-1)*(k+1)+(NI-1)*j+i;
// 				indexB=(NJ-1)*(NI-1)*(k-1)+(NI-1)*j+i;
			
				indexE=indexP+1;
				indexW=indexP-1;
				indexN=indexP+J_Increment;
				indexS=indexP-J_Increment;
				indexT=indexP+K_Increment;
				indexB=indexP-K_Increment;

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

			
				// west face
				*(F+face_index) = (uw*(*(SurfArea+area_index+0))+vw*(*(SurfArea+area_index+1))+ww*(*(SurfArea+area_index+2)));
				
				// east face
				*(F+face_index+1) = (ue*(*(SurfArea+area_index+3))+ve*(*(SurfArea+area_index+4))+we*(*(SurfArea+area_index+5)));
				
				// south face
				*(F+face_index+2) = (us*(*(SurfArea+area_index+6))+vs*(*(SurfArea+area_index+7))+ws*(*(SurfArea+area_index+8)));
				
				// north face
				*(F+face_index+3) = (un*(*(SurfArea+area_index+9))+vn*(*(SurfArea+area_index+10))+wn*(*(SurfArea+area_index+11)));
				
				// bottom face
				*(F+face_index+4) = (ub*(*(SurfArea+area_index+12))+vb*(*(SurfArea+area_index+13))+wb*(*(SurfArea+area_index+14)));
				
				// top face
				*(F+face_index+5) = (ut*(*(SurfArea+area_index+15))+vt*(*(SurfArea+area_index+16))+wt*(*(SurfArea+area_index+17)));
				
				indexP++;
				face_index = face_index + 6;
				area_index = area_index + 18;
				//fprintf(fp1,"west=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
				//fprintf(fp1,"east=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
				//fprintf(fp1,"south=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
				//fprintf(fp1,"north=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
				//fprintf(fp1,"bottom=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
				//fprintf(fp1,"top=%lf\n",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
			}	/*End of looping over NI,NJ,NK*/
		}
		
		else if(base->flagMPFV==1)//Adde by sushil For multiphase flow problem 
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
				
// 				indexE=(NJ-1)*(NI-1)*k+(NI-1)*j+(i+1);
// 				indexW=(NJ-1)*(NI-1)*k+(NI-1)*j+(i-1);
// 				indexN=(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i;
// 				indexS=(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i;
// 				indexT=(NJ-1)*(NI-1)*(k+1)+(NI-1)*j+i;
// 				indexB=(NJ-1)*(NI-1)*(k-1)+(NI-1)*j+i;
				
				indexE=indexP+1;
				indexW=indexP-1;
				indexN=indexP+J_Increment;
				indexS=indexP-J_Increment;
				indexT=indexP+K_Increment;
				indexB=indexP-K_Increment;
	
				u_P=*(domainUmass+indexP);
				v_P=*(domainVmass+indexP);
				w_P=*(domainWmass+indexP);
				vol_P=*(domainVolume+indexP);
				
				
				rho_P=1/(*(domainRho+indexP));
				u_P = u_P*rho_P;
				v_P = v_P*rho_P;
				w_P = w_P*rho_P;
					
				
				if(i==(NI-2))
				{
					u_E=*(slabUmass[0]+(NJ-1)*k+j);
					v_E=*(slabVmass[0]+(NJ-1)*k+j);
					w_E=*(slabWmass[0]+(NJ-1)*k+j);
					vol_E=*(slabVolume[0]+(NJ-1)*k+j);
					
					
					rho_E = 1/(*(slabRho[0]+(NJ-1)*k+j));
					u_E = u_E*rho_E;
					v_E = v_E*rho_E;
					w_E = w_E*rho_E;
					
				}
				else
				{
					u_E=*(domainUmass+indexE);
					v_E=*(domainVmass+indexE);
					w_E=*(domainWmass+indexE);
					vol_E=*(domainVolume+indexE);
					
					
					rho_E = 1/(*(domainRho+indexE));
					u_E = u_E*rho_E;
					v_E = v_E*rho_E;
					w_E = w_E*rho_E;
					
					
				}
	
				if(i==0)
				{
					u_W =*(slabUmass[1]+(NJ-1)*k+j);
					v_W =*(slabVmass[1]+(NJ-1)*k+j);
					w_W =*(slabWmass[1]+(NJ-1)*k+j);
					vol_W =*(slabVolume[1]+(NJ-1)*k+j);
					
					
					rho_W = 1/(*(slabRho[1]+(NJ-1)*k+j));
					u_W = u_W*rho_W;
					v_W = v_W*rho_W;
					w_W = w_W*rho_W;
					
				}
				else
				{
					u_W=*(domainUmass+indexW);
					v_W=*(domainVmass+indexW);
					w_W=*(domainWmass+indexW);
					vol_W=*(domainVolume+indexW);
					
					
					rho_W = 1/(*(domainRho+indexW));
					u_W = u_W*rho_W;
					v_W = v_W*rho_W;
					w_W = w_W*rho_W;
					
				}
	
				if(j==(NJ-2))
				{
					u_N=*(slabUmass[2]+(NI-1)*k+i);
					v_N=*(slabVmass[2]+(NI-1)*k+i);
					w_N=*(slabWmass[2]+(NI-1)*k+i);
					vol_N=*(slabVolume[2]+(NI-1)*k+i);
					
					
					rho_N = 1/(*(slabRho[2]+(NI-1)*k+i));
					u_N = u_N*rho_N;
					v_N = v_N*rho_N;
					w_N = w_N*rho_N;
					
				}
				else
				{
					u_N=*(domainUmass+indexN);
					v_N=*(domainVmass+indexN);
					w_N=*(domainWmass+indexN);
					vol_N=*(domainVolume+indexN);
					
					
					rho_N = 1/(*(domainRho+indexN));
					u_N = u_N*rho_N;
					v_N = v_N*rho_N;
					w_N = w_N*rho_N;
					
				}
	
				if(j==0)
				{
					u_S=*(slabUmass[3]+(NI-1)*k+i);
					v_S=*(slabVmass[3]+(NI-1)*k+i);
					w_S=*(slabWmass[3]+(NI-1)*k+i);
					vol_S=*(slabVolume[3]+(NI-1)*k+i);
					
					
					rho_S = 1/(*(slabRho[3]+(NI-1)*k+i));
					u_S = u_S*rho_S;
					v_S = v_S*rho_S;
					w_S = w_S*rho_S;
					
				}
				else
				{
					u_S=*(domainUmass+indexS);
					v_S=*(domainVmass+indexS);
					w_S=*(domainWmass+indexS);
					vol_S=*(domainVolume+indexS);
					
					
					rho_S = 1/(*(domainRho+indexS));
					u_S = u_S*rho_S;
					v_S = v_S*rho_S;
					w_S = w_S*rho_S;
					
				}
	
				if(k==(NK-2))
				{
					u_T=*(slabUmass[4]+(NI-1)*j+i);
					v_T=*(slabVmass[4]+(NI-1)*j+i);
					w_T=*(slabWmass[4]+(NI-1)*j+i);
					vol_T=*(slabVolume[4]+(NI-1)*j+i);
					
					
					rho_T = 1/(*(slabRho[4]+(NI-1)*j+i));
					u_T = u_T*rho_T;
					v_T = v_T*rho_T;
					w_T = w_T*rho_T;
				
				}
				else
				{
					u_T=*(domainUmass+indexT);
					v_T=*(domainVmass+indexT);
					w_T=*(domainWmass+indexT);
					vol_T=*(domainVolume+indexT);
					
					
					rho_T = 1/(*(domainRho+indexT));
					u_T = u_T*rho_T;
					v_T = v_T*rho_T;
					w_T = w_T*rho_T;
				
				}
	
				if(k==0)
				{
					u_B=*(slabUmass[5]+(NI-1)*j+i);
					v_B=*(slabVmass[5]+(NI-1)*j+i);
					w_B=*(slabWmass[5]+(NI-1)*j+i);
					vol_B=*(slabVolume[5]+(NI-1)*j+i);
					
					
					rho_B = 1/(*(slabRho[5]+(NI-1)*j+i));
					u_B = u_B*rho_B;
					v_B = v_B*rho_B;
					w_B = w_B*rho_B;
					
				}
				else
				{
					u_B=*(domainUmass+indexB);
					v_B=*(domainVmass+indexB);
					w_B=*(domainWmass+indexB);
					vol_B=*(domainVolume+indexB);
					
					
					rho_B = 1/(*(domainRho+indexB));
					u_B = u_B*rho_B;
					v_B = v_B*rho_B;
					w_B = w_B*rho_B;
					
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
	
				
				// west face
				*(F+face_index) = (uw*(*(SurfArea+area_index+0))+vw*(*(SurfArea+area_index+1))+ww*(*(SurfArea+area_index+2)));
					
				// east face
				*(F+face_index+1) = (ue*(*(SurfArea+area_index+3))+ve*(*(SurfArea+area_index+4))+we*(*(SurfArea+area_index+5)));
					
				// south face
				*(F+face_index+2) = (us*(*(SurfArea+area_index+6))+vs*(*(SurfArea+area_index+7))+ws*(*(SurfArea+area_index+8)));
					
				// north face
				*(F+face_index+3) = (un*(*(SurfArea+area_index+9))+vn*(*(SurfArea+area_index+10))+wn*(*(SurfArea+area_index+11)));
					
				// bottom face
				*(F+face_index+4) = (ub*(*(SurfArea+area_index+12))+vb*(*(SurfArea+area_index+13))+wb*(*(SurfArea+area_index+14)));
					
				// top face
				*(F+face_index+5) = (ut*(*(SurfArea+area_index+15))+vt*(*(SurfArea+area_index+16))+wt*(*(SurfArea+area_index+17)));
					
				indexP++;
				face_index = face_index + 6;
				area_index = area_index + 18;
				//fprintf(fp1,"west=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
				//fprintf(fp1,"east=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
				//fprintf(fp1,"south=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
				//fprintf(fp1,"north=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
				//fprintf(fp1,"bottom=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
				//fprintf(fp1,"top=%lf\n",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
			}/*End of looping over NI,NJ,NK*/
	
		}	
	}/*End of looping over the zones */
	return 0;
}

/**
This function is used for calculating the volumeflux from the current converged
values of massvelocities. This function is used for 2D domains.
See Eqn No(3.36) Project Report No. 1, Dec 2005.
*/
int volumefluxstar2d(Base *base)
{
	int nZones,numZones;
	int k,j,i;
	int indexE,indexW,indexS,indexN,indexP;//,indexB,indexT;
	int face_index,area_index;

	double *domainVolume,*domainUmass,*domainVmass;
	double *slabVolume[6],*slabUmass[6],*slabVmass[6],*slabRho[6],*domainRho;
	double *SurfArea,*F;
	Zone* zone;
	nZones=(base->Number_of_zones);
	
	int NI,NJ,NK,*ptrSize;
	double uw,vw,ue,ve,un,vn,us,vs;//,ut,vt,ub,vb;
	double u_P,u_E,u_W,u_S,u_N;//,u_T,u_B;
	double v_P,v_E,v_W,v_S,v_N;//,v_T,v_B;
	double vol_P,vol_E,vol_W,vol_N,vol_S,rho_P,rho_E,rho_W,rho_S,rho_N;//,rho_T,rho_B;//,vol_T,vol_B;
// 	FILE *fp;
// 	double sum;
// 	fp = fopen("volumeflux.dat","w");
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
			//domainWmass=MakeDomainArray(zone,300);
			domainVolume=MakeDomainArray(zone,5);
	
			MakeSlabPtrArray(zone,slabUmass,100);
			MakeSlabPtrArray(zone,slabVmass,200);
			//MakeSlabPtrArray(zone,slabWmass,300);
			MakeSlabPtrArray(zone,slabVolume,5);
			
			indexP=0;
			face_index=0;
			area_index=0;
		
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				indexE=indexP+1;
				indexW=indexP-1;
				indexN=indexP+(NI-1);
				indexS=indexP-(NI-1);
			
// 				indexN=(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i;
// 				indexS=(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i;
// 				indexT=(NJ-1)*(NI-1)*(k+1)+(NI-1)*j+i;
// 				indexB=(NJ-1)*(NI-1)*(k-1)+(NI-1)*j+i;

				u_P=*(domainUmass+indexP);
				v_P=*(domainVmass+indexP);
				//w_P=*(domainWmass+indexP);
				vol_P=*(domainVolume+indexP);
				
				
				if(i==(NI-2))
				{
					u_E=*(slabUmass[0]+j);
					v_E=*(slabVmass[0]+j);
				//w_E=*(slabWmass[0]+(NJ-1)*k+j);
					vol_E=*(slabVolume[0]+j);
				
				}
				else
				{
					u_E=*(domainUmass+indexE);
					v_E=*(domainVmass+indexE);
				//w_E=*(domainWmass+indexE);
					vol_E=*(domainVolume+indexE);
				
				}

				if(i==0)
				{
					u_W =*(slabUmass[1]+j);
					v_W =*(slabVmass[1]+j);
				//w_W =*(slabWmass[1]+(NJ-1)*k+j);
					vol_W =*(slabVolume[1]+j);
				
				}
				else
				{
					u_W=*(domainUmass+indexW);
					v_W=*(domainVmass+indexW);
				//w_W=*(domainWmass+indexW);
					vol_W=*(domainVolume+indexW);
				
				}

				if(j==(NJ-2))
				{
					u_N=*(slabUmass[2]+i);
					v_N=*(slabVmass[2]+i);
				//w_N=*(slabWmass[2]+(NI-1)*k+i);
					vol_N=*(slabVolume[2]+i);
				
				}
				else
				{
					u_N=*(domainUmass+indexN);
					v_N=*(domainVmass+indexN);
				//w_N=*(domainWmass+indexN);
					vol_N=*(domainVolume+indexN);
				
				}

				if(j==0)
				{
					u_S=*(slabUmass[3]+i);
					v_S=*(slabVmass[3]+i);
				//w_S=*(slabWmass[3]+(NI-1)*k+i);
					vol_S=*(slabVolume[3]+i);
				
				}
				else
				{
					u_S=*(domainUmass+indexS);
					v_S=*(domainVmass+indexS);
				//w_S=*(domainWmass+indexS);
					vol_S=*(domainVolume+indexS);
				
				}

// 				if(k==(NK-2))
// 				{
// 				u_T=*(slabUmass[4]+(NI-1)*j+i);
// 				v_T=*(slabVmass[4]+(NI-1)*j+i);
// 				//w_T=*(slabWmass[4]+(NI-1)*j+i);
// 				vol_T=*(slabVolume[4]+(NI-1)*j+i);
// 				}
// 				else
// 				{
// 				u_T=*(domainUmass+indexT);
// 				v_T=*(domainVmass+indexT);
// 				//w_T=*(domainWmass+indexT);
// 				vol_T=*(domainVolume+indexT);
// 				}
				// 
// 				if(k==0)
// 				{
// 				u_B=*(slabUmass[5]+(NI-1)*j+i);
// 				v_B=*(slabVmass[5]+(NI-1)*j+i);
// 				//w_B=*(slabWmass[5]+(NI-1)*j+i);
// 				vol_B=*(slabVolume[5]+(NI-1)*j+i);
// 				}
// 				else
// 				{
// 				u_B=*(domainUmass+indexB);
// 				v_B=*(domainVmass+indexB);
// 				//w_B=*(domainWmass+indexB);
// 				vol_B=*(domainVolume+indexB);
// 				}


				ue=((u_E*vol_P+u_P*vol_E)/(vol_P+vol_E));
				uw=((u_W*vol_P+u_P*vol_W)/(vol_P+vol_W));
				us=((u_S*vol_P+u_P*vol_S)/(vol_P+vol_S));
				un=((u_N*vol_P+u_P*vol_N)/(vol_P+vol_N));
// 				ut=((u_T*vol_P+u_P*vol_T)/(vol_P+vol_T));
// 				ub=((u_B*vol_P+u_P*vol_B)/(vol_P+vol_B));

				ve=((v_E*vol_P+v_P*vol_E)/(vol_P+vol_E));
				vw=((v_W*vol_P+v_P*vol_W)/(vol_P+vol_W));
				vs=((v_S*vol_P+v_P*vol_S)/(vol_P+vol_S));
				vn=((v_N*vol_P+v_P*vol_N)/(vol_P+vol_N));
// 				vt=((v_T*vol_P+v_P*vol_T)/(vol_P+vol_T));
// 				vb=((v_B*vol_P+v_P*vol_B)/(vol_P+vol_B));
				//printf("\nus=%e\tun=%e",us,un);
			
				*(F+face_index)   = (uw*(*(SurfArea+area_index))  +vw*(*(SurfArea+area_index+1 )));
				//+ww*(*(SurfArea+Index+3*0+2)));
				*(F+face_index+1) = (ue*(*(SurfArea+area_index+3))+ve*(*(SurfArea+area_index+4 )));
				//+we*(*(SurfArea+Index+3*1+2)));
				*(F+face_index+2) = (us*(*(SurfArea+area_index+6))+vs*(*(SurfArea+area_index+7 )));
				//+ws*(*(SurfArea+Index+3*2+2)));
				*(F+face_index+3) = (un*(*(SurfArea+area_index+9))+vn*(*(SurfArea+area_index+10)));
				//+wn*(*(SurfArea+Index+3*3+2));
				*(F+face_index+4) = 0.0;
				//(ub*(*(SurfArea+Index+3*4+0))+vb*(*(SurfArea+Index+3*4+1))+wb*(*(SurfArea+Index+3*4+2)));
				*(F+face_index+5) =0.0;
				//(ut*(*(SurfArea+Index+3*5+0))+vt*(*(SurfArea+Index+3*5+1))+wt*(*(SurfArea+Index+3*5+2)));
				//sum = 0 ;
				//sum = sum + *(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+0)+ *(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+1)+*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+2)+*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+3)+*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+4)+*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+5);
				//fprintf(fp,"%e\t%e\t%e\t%e\t%e\t%e\t%e\n",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+0),*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+1),*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+2),*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+3),*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+4),*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+5),sum);
			
			
				indexP++;
				face_index = face_index + 6;
				area_index = area_index + 18;
			}/*End of looping over NI,NJ,NK*/
		}
		else if(base->flagMPFV==1)
		{
			domainUmass=MakeDomainArray(zone,3100);
			domainVmass=MakeDomainArray(zone,3200);
			domainRho=MakeDomainArray(zone,710);
			domainVolume=MakeDomainArray(zone,5);
	
			MakeSlabPtrArray(zone,slabUmass,3100);
			MakeSlabPtrArray(zone,slabVmass,3200);
			MakeSlabPtrArray(zone,slabRho,710);
			MakeSlabPtrArray(zone,slabVolume,5);
			
			indexP=0;
			face_index=0;
			area_index=0;
		
			for(k=0;k<NK-1;k++)
			for(j=0;j<NJ-1;j++)
			for(i=0;i<NI-1;i++)
			{
				indexE=indexP+1;
				indexW=indexP-1;
				indexN=indexP+(NI-1);
				indexS=indexP-(NI-1);
			
	
				u_P=*(domainUmass+indexP);
				v_P=*(domainVmass+indexP);
				//w_P=*(domainWmass+indexP);
				vol_P=*(domainVolume+indexP);
				
				rho_P=1/(*(domainRho+indexP));
				u_P = u_P*rho_P;
				v_P = v_P*rho_P;

				if(i==(NI-2))
				{
					u_E=*(slabUmass[0]+j);
					v_E=*(slabVmass[0]+j);
				//w_E=*(slabWmass[0]+(NJ-1)*k+j);
					vol_E=*(slabVolume[0]+j);
					
					rho_E = 1/(*(slabRho[0]+(NJ-1)*k+j));
					u_E = u_E*rho_E;
					v_E = v_E*rho_E;
					
					
					
				}
				else
				{
					u_E=*(domainUmass+indexE);
					v_E=*(domainVmass+indexE);
					//w_E=*(domainWmass+indexE);
					vol_E=*(domainVolume+indexE);
					
					
					rho_E = 1/(*(domainRho+indexE));
					u_E = u_E*rho_E;
					v_E = v_E*rho_E;
				
					
				}

				if(i==0)
				{
					u_W =*(slabUmass[1]+j);
					v_W =*(slabVmass[1]+j);
					//w_W =*(slabWmass[1]+(NJ-1)*k+j);
					vol_W =*(slabVolume[1]+j);

					rho_W = 1/(*(slabRho[1]+(NJ-1)*k+j));
					u_W = u_W*rho_W;
					v_W = v_W*rho_W;

				}
				else
				{
					u_W=*(domainUmass+indexW);
					v_W=*(domainVmass+indexW);
					//w_W=*(domainWmass+indexW);
					vol_W=*(domainVolume+indexW);
					
					
					rho_W = 1/(*(domainRho+indexW));
					u_W = u_W*rho_W;
					v_W = v_W*rho_W;
					
				}

				if(j==(NJ-2))
				{
					u_N=*(slabUmass[2]+i);
					v_N=*(slabVmass[2]+i);
					//w_N=*(slabWmass[2]+(NI-1)*k+i);
					vol_N=*(slabVolume[2]+i);
					
					
					rho_N = 1/(*(slabRho[2]+(NI-1)*k+i));
					u_N = u_N*rho_N;
					v_N = v_N*rho_N;
					
					
				}
				else
				{
					u_N=*(domainUmass+indexN);
					v_N=*(domainVmass+indexN);
					//w_N=*(domainWmass+indexN);
					vol_N=*(domainVolume+indexN);
					
					
					rho_N = 1/(*(domainRho+indexN));
					u_N = u_N*rho_N;
					v_N = v_N*rho_N;
					
					
				}

				if(j==0)
				{
					u_S=*(slabUmass[3]+i);
					v_S=*(slabVmass[3]+i);
					//w_S=*(slabWmass[3]+(NI-1)*k+i);
					vol_S=*(slabVolume[3]+i);
					
					
					rho_S = 1/(*(slabRho[3]+(NI-1)*k+i));
					u_S = u_S*rho_S;
					v_S = v_S*rho_S;
					
					
				}
				else
				{
					u_S=*(domainUmass+indexS);
					v_S=*(domainVmass+indexS);
					//w_S=*(domainWmass+indexS);
					vol_S=*(domainVolume+indexS);
					
					
					rho_S = 1/(*(domainRho+indexS));
					u_S = u_S*rho_S;
					v_S = v_S*rho_S;
				
					
				}

	

				ue=((u_E*vol_P+u_P*vol_E)/(vol_P+vol_E));
				uw=((u_W*vol_P+u_P*vol_W)/(vol_P+vol_W));
				us=((u_S*vol_P+u_P*vol_S)/(vol_P+vol_S));
				un=((u_N*vol_P+u_P*vol_N)/(vol_P+vol_N));


				ve=((v_E*vol_P+v_P*vol_E)/(vol_P+vol_E));
				vw=((v_W*vol_P+v_P*vol_W)/(vol_P+vol_W));
				vs=((v_S*vol_P+v_P*vol_S)/(vol_P+vol_S));
				vn=((v_N*vol_P+v_P*vol_N)/(vol_P+vol_N));

			
				*(F+face_index)   = (uw*(*(SurfArea+area_index))  +vw*(*(SurfArea+area_index+1 )));
				//+ww*(*(SurfArea+Index+3*0+2)));
				*(F+face_index+1) = (ue*(*(SurfArea+area_index+3))+ve*(*(SurfArea+area_index+4 )));
				//+we*(*(SurfArea+Index+3*1+2)));
				*(F+face_index+2) = (us*(*(SurfArea+area_index+6))+vs*(*(SurfArea+area_index+7 )));
				//+ws*(*(SurfArea+Index+3*2+2)));
				*(F+face_index+3) = (un*(*(SurfArea+area_index+9))+vn*(*(SurfArea+area_index+10)));
				//+wn*(*(SurfArea+Index+3*3+2));
				*(F+face_index+4) = 0.0;
				//(ub*(*(SurfArea+Index+3*4+0))+vb*(*(SurfArea+Index+3*4+1))+wb*(*(SurfArea+Index+3*4+2)));
				*(F+face_index+5) =0.0;
				//(ut*(*(SurfArea+Index+3*5+0))+vt*(*(SurfArea+Index+3*5+1))+wt*(*(SurfArea+Index+3*5+2)));
				//sum = 0 ;
				//sum = sum + *(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+0)+ *(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+1)+*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+2)+*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+3)+*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+4)+*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+5);
				//fprintf(fp,"%e\t%e\t%e\t%e\t%e\t%e\t%e\n",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+0),*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+1),*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+2),*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+3),*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+4),*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+5),sum);
			
			
				indexP++;
				face_index = face_index + 6;
				area_index = area_index + 18;
			}/*End of looping over NI,NJ,NK*/
		}
		
	}/*End of looping over the zones */
	//fclose(fp);
	return 0;
}
