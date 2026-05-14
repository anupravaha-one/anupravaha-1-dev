#include "APKvolumefluxDisperse.h"
/**
These functions are used to compute the volume flux for the disperse phase variables. It is called only if disperse phase is being solved.
*/
int volumefluxDisperse2d(Base *base)
{
	int nZones,numZones;
	int k,j,i;
	int indexB,indexT,indexE,indexW,indexS,indexN,indexP,flagrho;

	double *domainVolume,*domainu,*domainv,*domainVartheta;
	double *slabVolume[6],*slabu[6],*slabv[6];	
	double *SurfArea,*F,A,B,*FG, total;
	Zone* zone;
	nZones=(base->Number_of_zones);
	
	
	int K_Increment,J_Increment,I_Increment,face_index,area_index;
	
	int NI,NJ,NK,*ptrSize;
	double uw,vw,ue,ve,un,vn,us,vs;
	
	double U_P,V_P,U_E,U_W,U_N,U_S,U_T,U_B,V_E,V_W,V_N,V_S,V_T,V_B;

	double vol_P,vol_E,vol_W,vol_N,vol_S,vol_T,vol_B;
	
	//shash
	double EPSILONU;
	EPSILONU = base->BaseParameters->EPSILONDrop;
	
	for(numZones=0;numZones<nZones;numZones++)
	{
		zone=(base->CgnsZone+numZones);
		
		ptrSize=zone->Size;

		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);

		F = zone->Flowsol->MassfluxDisperse;
		FG = zone->Flowsol->Massflux;

		SurfArea=zone->ZGeom->SurfaceArea;
                 
		
		domainu=MakeDomainArray(zone,62);
		domainv=MakeDomainArray(zone,63);
		
		domainVolume=MakeDomainArray(zone,5);
		domainVartheta=MakeDomainArray(zone,650);

		MakeSlabPtrArray(zone,slabu,62);
		MakeSlabPtrArray(zone,slabv,63);
		
		
		
		MakeSlabPtrArray(zone,slabVolume,5);
		

//#######################################################################
		K_Increment = (NI-1)*(NJ-1);
		J_Increment = (NI-1);
		I_Increment = 1;

		indexP=0;
		face_index=0;
		area_index=0;

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

			U_P=*(domainu+indexP);
			V_P=*(domainv+indexP);
			
			vol_P=*(domainVolume+indexP);

			if(i==(NI-2))
			{
				U_E=*(slabu[0]+(NJ-1)*k+j);
				V_E=*(slabv[0]+(NJ-1)*k+j);
				
				vol_E=*(slabVolume[0]+(NJ-1)*k+j);
			}
			else
			{
				U_E=*(domainu+indexE);
				V_E=*(domainv+indexE);
				
				vol_E=*(domainVolume+indexE);
			}

			if(i==0)
			{
				U_W =*(slabu[1]+(NJ-1)*k+j);
				V_W =*(slabv[1]+(NJ-1)*k+j);
				
				vol_W =*(slabVolume[1]+(NJ-1)*k+j);
			}
			else
			{
				U_W=*(domainu+indexW);
				V_W=*(domainv+indexW);
				
				vol_W=*(domainVolume+indexW);
			}

			if(j==(NJ-2))
			{
				U_N=*(slabu[2]+(NI-1)*k+i);
				V_N=*(slabv[2]+(NI-1)*k+i);
				
				vol_N=*(slabVolume[2]+(NI-1)*k+i);
			}
			else
			{
				U_N=*(domainu+indexN);
				V_N=*(domainv+indexN);
				
				vol_N=*(domainVolume+indexN);
			}

			if(j==0)
			{
				U_S=*(slabu[3]+(NI-1)*k+i);
				V_S=*(slabv[3]+(NI-1)*k+i);
				
				vol_S=*(slabVolume[3]+(NI-1)*k+i);
			}
			else
			{
				U_S=*(domainu+indexS);
				V_S=*(domainv+indexS);
				
				vol_S=*(domainVolume+indexS);
				
			}

			if(k==(NK-2))
			{
				U_T=*(slabu[4]+(NI-1)*j+i);
				V_T=*(slabv[4]+(NI-1)*j+i);
				
				vol_T=*(slabVolume[4]+(NI-1)*j+i);
			}
			else
			{
				U_T=*(domainu+indexT);
				V_T=*(domainv+indexT);
				
				vol_T=*(domainVolume+indexT);
			}

			if(k==0)
			{
				U_B=*(slabu[5]+(NI-1)*j+i);
				V_B=*(slabv[5]+(NI-1)*j+i);
				
				vol_B=*(slabVolume[5]+(NI-1)*j+i);
			}
			else
			{
				U_B=*(domainu+indexB);
				V_B=*(domainv+indexB);
				
				vol_B=*(domainVolume+indexB);
			}

			ue=((U_E*vol_P+U_P*vol_E)/(vol_P+vol_E));
			uw=((U_W*vol_P+U_P*vol_W)/(vol_P+vol_W));
			us=((U_S*vol_P+U_P*vol_S)/(vol_P+vol_S));
			un=((U_N*vol_P+U_P*vol_N)/(vol_P+vol_N));
			

			ve=((V_E*vol_P+V_P*vol_E)/(vol_P+vol_E));
			vw=((V_W*vol_P+V_P*vol_W)/(vol_P+vol_W));
			vs=((V_S*vol_P+V_P*vol_S)/(vol_P+vol_S));
			vn=((V_N*vol_P+V_P*vol_N)/(vol_P+vol_N));
			
			
			if(*(domainVartheta+indexP)>EPSILONU)
			{
				*(F+face_index)   = (uw*(*(SurfArea+area_index))  +vw*(*(SurfArea+area_index+1 )));
		
				*(F+face_index+1) = (ue*(*(SurfArea+area_index+3))+ve*(*(SurfArea+area_index+4 )));
			
				*(F+face_index+2) = (us*(*(SurfArea+area_index+6))+vs*(*(SurfArea+area_index+7 )));
			
				*(F+face_index+3) = (un*(*(SurfArea+area_index+9))+vn*(*(SurfArea+area_index+10)));
		
				*(F+face_index+4) = 0.0;
			
				*(F+face_index+5) =0.0;

			/*	if((i==1)&&(j==(int)NJ/4 || j==(int)NJ/2 || j==(int)3*NJ/4))
				{	
					total = *(F+face_index) + *(F+face_index+1) + *(F+face_index+2) + *(F+face_index+3) + *(F+face_index+4) + *(F+face_index+5);
					printf("\n For i=%d, j=%d \t Sum of Fdf = %e\n",i,j,total);
				} */

			}
			else
			{
				
				*(F+face_index)=*(FG+face_index);
				*(F+face_index+1)=*(FG+face_index+1);
				*(F+face_index+2)=*(FG+face_index+2);
				*(F+face_index+3)=*(FG+face_index+3);
				*(F+face_index+4)= 0.0;
				*(F+face_index+5)=0.0;

				*(zone->Flowsol->U_Disperse+indexP)=*(zone->Flowsol->u+indexP);
				*(zone->Flowsol->V_Disperse+indexP)=*(zone->Flowsol->v+indexP);
			}
			
			if((base->flagParticle==1)&&base->flagTurb==1)
			{
				//if((Var_P> EPSILONU))
				//{
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
				//}
				/*else if((Var_P< EPSILONU))
				{
				//Equating gas flux to particle flux when volume fraction is less than the limiting value
				*(F+face_index)=*(FG+face_index);
				*(F+face_index+1)=*(FG+face_index+1);
				*(F+face_index+2)=*(FG+face_index+2);
				*(F+face_index+3)=*(FG+face_index+3);
				*(F+face_index+4)= 0.0;
				*(F+face_index+5)=0.0;
				}*/
			}
			indexP++;
			face_index = face_index + 6;
			area_index = area_index + 18;		
		}/*End of looping over NI,NJ,NK*/

	
	}/*End of looping over the zones */
	return 0;
}

/**
This function computes the volume flux for disperse phase in 3d geometry
*/
int volumefluxDisperse(Base *base)
{
	
	int nZones,numZones;
	int k,j,i;
	int indexB,indexT,indexE,indexW,indexS,indexN,indexP;

	double *domainVolume,*domainu,*domainv, *domainw;
	double *slabVolume[6],*slabu[6],*slabv[6], *slabw[6];	
	double *SurfArea,*F,A,B,*FG;
	Zone* zone;
	nZones=(base->Number_of_zones);
	
	
	int K_Increment,J_Increment,I_Increment,face_index,area_index;
	
	int NI,NJ,NK,*ptrSize;
	
	double uw,vw,ue,ve,un,vn,us,vs, ut,vt,ub,vb, ww, we, wn, ws, wt, wb;

	double U_P,V_P,U_E,U_W,U_N,U_S,U_T,U_B,V_E,V_W,V_N,V_S,V_T,V_B,W_P,W_W,W_E,W_N,W_S,W_T,W_B;

	double vol_P,vol_E,vol_W,vol_N,vol_S,vol_T,vol_B;
	
	//shash
	double EPSILONU;
	EPSILONU = base->BaseParameters->EPSILONDrop;
	
	for(numZones=0;numZones<nZones;numZones++)
	{
		zone=(base->CgnsZone+numZones);
		
		ptrSize=zone->Size;

		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);

		F = zone->Flowsol->MassfluxDisperse;
		FG = zone->Flowsol->Massflux;

		SurfArea=zone->ZGeom->SurfaceArea;
                 
		
		domainu=MakeDomainArray(zone,62);
		domainv=MakeDomainArray(zone,63);
		domainw=MakeDomainArray(zone,64);
		
		domainVolume=MakeDomainArray(zone,5);
		
		MakeSlabPtrArray(zone,slabu,62);
		MakeSlabPtrArray(zone,slabv,63);
		MakeSlabPtrArray(zone,slabw,64);
		
		MakeSlabPtrArray(zone,slabVolume,5);
		

//#######################################################################
		K_Increment = (NI-1)*(NJ-1);
		J_Increment = (NI-1);
		I_Increment = 1;

		indexP=0;
		face_index=0;
		area_index=0;
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

			U_P=*(domainu+indexP);
			V_P=*(domainv+indexP);
			W_P=*(domainw+indexP);
		
			vol_P=*(domainVolume+indexP);

			if(i==(NI-2))
			{
				U_E=*(slabu[0]+(NJ-1)*k+j);
				V_E=*(slabv[0]+(NJ-1)*k+j);
				W_E=*(slabw[0]+(NJ-1)*k+j);
				
				vol_E=*(slabVolume[0]+(NJ-1)*k+j);
			}
			else
			{
				U_E=*(domainu+indexE);
				V_E=*(domainv+indexE);
				W_E=*(domainw+indexE);
				
				vol_E=*(domainVolume+indexE);
			}

			if(i==0)
			{
				U_W =*(slabu[1]+(NJ-1)*k+j);
				V_W =*(slabv[1]+(NJ-1)*k+j);
				W_W =*(slabw[1]+(NJ-1)*k+j);
				
				vol_W =*(slabVolume[1]+(NJ-1)*k+j);
			}
			else
			{
				U_W=*(domainu+indexW);
				V_W=*(domainv+indexW);
				W_W=*(domainw+indexW);

				vol_W=*(domainVolume+indexW);
			}

			if(j==(NJ-2))
			{
				U_N=*(slabu[2]+(NI-1)*k+i);
				V_N=*(slabv[2]+(NI-1)*k+i);
				W_N=*(slabw[2]+(NI-1)*k+i);

				vol_N=*(slabVolume[2]+(NI-1)*k+i);
			}
			else
			{
				U_N=*(domainu+indexN);
				V_N=*(domainv+indexN);
				W_N=*(domainw+indexN);

				vol_N=*(domainVolume+indexN);
			}

			if(j==0)
			{
				U_S=*(slabu[3]+(NI-1)*k+i);
				V_S=*(slabv[3]+(NI-1)*k+i);
				W_S=*(slabw[3]+(NI-1)*k+i);
				
				vol_S=*(slabVolume[3]+(NI-1)*k+i);
			}
			else
			{
				U_S=*(domainu+indexS);
				V_S=*(domainv+indexS);
				W_S=*(domainw+indexS);

				vol_S=*(domainVolume+indexS);
				
			}

			if(k==(NK-2))
			{
				U_T=*(slabu[4]+(NI-1)*j+i);
				V_T=*(slabv[4]+(NI-1)*j+i);
				W_T=*(slabw[4]+(NI-1)*j+i);
				
				vol_T=*(slabVolume[4]+(NI-1)*j+i);
			}
			else
			{
				U_T=*(domainu+indexT);
				V_T=*(domainv+indexT);
				W_T=*(domainw+indexT);
				
				vol_T=*(domainVolume+indexT);
			}

			if(k==0)
			{
				U_B=*(slabu[5]+(NI-1)*j+i);
				V_B=*(slabv[5]+(NI-1)*j+i);
				W_B=*(slabw[5]+(NI-1)*j+i);
				
				vol_B=*(slabVolume[5]+(NI-1)*j+i);
			}
			else
			{
				U_B=*(domainu+indexB);
				V_B=*(domainv+indexB);
				W_B=*(domainw+indexB);

				vol_B=*(domainVolume+indexB);
				
			}


			ue=((U_E*vol_P+U_P*vol_E)/(vol_P+vol_E));
			uw=((U_W*vol_P+U_P*vol_W)/(vol_P+vol_W));
			us=((U_S*vol_P+U_P*vol_S)/(vol_P+vol_S));
			un=((U_N*vol_P+U_P*vol_N)/(vol_P+vol_N));
			ut=((U_T*vol_P+U_P*vol_T)/(vol_P+vol_T));
			ub=((U_B*vol_P+U_P*vol_B)/(vol_P+vol_B));

			ve=((V_E*vol_P+V_P*vol_E)/(vol_P+vol_E));
			vw=((V_W*vol_P+V_P*vol_W)/(vol_P+vol_W));
			vs=((V_S*vol_P+V_P*vol_S)/(vol_P+vol_S));
			vn=((V_N*vol_P+V_P*vol_N)/(vol_P+vol_N));
			vt=((V_T*vol_P+V_P*vol_T)/(vol_P+vol_T));
			vb=((V_B*vol_P+V_P*vol_B)/(vol_P+vol_B));

			we=((W_E*vol_P+W_P*vol_E)/(vol_P+vol_E));
			ww=((W_W*vol_P+W_P*vol_W)/(vol_P+vol_W));
			ws=((W_S*vol_P+W_P*vol_S)/(vol_P+vol_S));
			wn=((W_N*vol_P+W_P*vol_N)/(vol_P+vol_N));
			wt=((W_T*vol_P+W_P*vol_T)/(vol_P+vol_T));
			wb=((W_B*vol_P+W_P*vol_B)/(vol_P+vol_B));

			*(F+face_index)   = (uw*(*(SurfArea+area_index))  +vw*(*(SurfArea+area_index+1 ))+ww*(*(SurfArea+area_index+3*0+2)));
			*(F+face_index+1) = (ue*(*(SurfArea+area_index+3))+ve*(*(SurfArea+area_index+4 ))+we*(*(SurfArea+area_index+3*1+2)));
			*(F+face_index+2) = (us*(*(SurfArea+area_index+6))+vs*(*(SurfArea+area_index+7 ))+ws*(*(SurfArea+area_index+3*2+2)));
			*(F+face_index+3) = (un*(*(SurfArea+area_index+9))+vn*(*(SurfArea+area_index+10))+wn*(*(SurfArea+area_index+3*3+2)));
			*(F+face_index+4) = (ub*(*(SurfArea+area_index+3*4+0))+vb*(*(SurfArea+area_index+3*4+1))+wb*(*(SurfArea+area_index+3*4+2)));
			*(F+face_index+5) = (ut*(*(SurfArea+area_index+3*5+0))+vt*(*(SurfArea+area_index+3*5+1))+wt*(*(SurfArea+area_index+3*5+2)));
			
			/*if((base->flagParticle==1)&&base->flagTurb==1)
			{
				if((Var_P> EPSILONU))
				{
				*(F+face_index)   = (uw*(*(SurfArea+area_index))  +vw*(*(SurfArea+area_index+1 ))+ww*(*(SurfArea+area_index+3*0+2)));
				*(F+face_index+1) = (ue*(*(SurfArea+area_index+3))+ve*(*(SurfArea+area_index+4 ))+we*(*(SurfArea+area_index+3*1+2)));
				*(F+face_index+2) = (us*(*(SurfArea+area_index+6))+vs*(*(SurfArea+area_index+7 ))+ws*(*(SurfArea+area_index+3*2+2)));
				*(F+face_index+3) = (un*(*(SurfArea+area_index+9))+vn*(*(SurfArea+area_index+10))+wn*(*(SurfArea+area_index+3*3+2));
				*(F+face_index+4) = (ub*(*(SurfArea+area_index+3*4+0))+vb*(*(SurfArea+area_index+3*4+1))+wb*(*(SurfArea+area_index+3*4+2)));
				*(F+face_index+5) = (ut*(*(SurfArea+area_index+3*5+0))+vt*(*(SurfArea+area_index+3*5+1))+wt*(*(SurfArea+area_index+3*5+2)));
				}
				else if((Var_P< EPSILONU))
				{
				//Equating gas flux to particle flux when volume fraction is less than the limiting value
				*(F+face_index)=*(FG+face_index);
				*(F+face_index+1)=*(FG+face_index+1);
				*(F+face_index+2)=*(FG+face_index+2);
				*(F+face_index+3)=*(FG+face_index+3);
				*(F+face_index+4)= *(FG+face_index+4);
				*(F+face_index+5)=*(FG+face_index+5);
				}
			}*/
			indexP++;
			face_index = face_index + 6;
			area_index = area_index + 18;		
		}/*End of looping over NI,NJ,NK*/

		
	}/*End of looping over the zones */
	return 0;
}


