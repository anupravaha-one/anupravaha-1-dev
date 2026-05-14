#include "APKvolumefluxMPFV.h"
/**
These functions are used to compute the volume flux for the disperse phase variables. It is called only if disperse phase is being solved.
*/
int volumefluxMPFV2d(Base *base)
{
	int nZones,numZones;
	int k,j,i;
	int indexB,indexT,indexE,indexW,indexS,indexN,indexP,flagrho;

	double *domainVolume,*domainMomentumx,*domainMomentumy,*domainVartheta,*rhoPhase;
	double *slabVolume[6],*slabMomentumx[6],*slabMomentumy[6],*slabVartheta[6];	
	double *SurfArea,*F;
	Zone* zone;
	nZones=(base->Number_of_zones);
	
	
	int K_Increment,J_Increment,I_Increment,face_index,area_index;
	
	int NI,NJ,NK,*ptrSize;
	double uw,vw,ue,ve,un,vn,us,vs;//,ut,vt,ub,vb;
	double Mxw,Myw,Mxe,Mye,Mxn,Myn,Mxs,Mys,Mxt,Myt,Mxb,Myb;
	double Varw,Vare,Varn,Vars,Vart,Varb;
	
	double Mx_P,Mx_E,Mx_W,Mx_S,Mx_N,Mx_T,Mx_B,My_P,My_E,My_W,My_S,My_N,My_T,My_B;
	double Var_P,Var_E,Var_W,Var_S,Var_N,Var_T,Var_B;
	double vol_P,vol_E,vol_W,vol_N,vol_S,vol_T,vol_B;
	

	
	for(numZones=0;numZones<nZones;numZones++)
	{
		zone=(base->CgnsZone+numZones);
		
		ptrSize=zone->Size;

		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);

		F = zone->PrevFlowsol->MassfluxDisperse;
		
		flagrho=base->BaseParameters->flagRhoPhase[1];
		rhoPhase=zone->Flowsol->RhoofPhase[1];//disperse phase density
		
		SurfArea=zone->ZGeom->SurfaceArea;
                
		
		domainMomentumx=MakeDomainArray(zone,44);//V_dx
		domainMomentumy=MakeDomainArray(zone,45);//V_dy
		
		domainVolume=MakeDomainArray(zone,5);
		
		domainVartheta=MakeDomainArray(zone,35);//prevflowsol--vartheta_d...

	
		
		
		MakeSlabPtrArray(zone,slabMomentumx,44);//V_dx
		MakeSlabPtrArray(zone,slabMomentumy,45);//V_dy
	
		MakeSlabPtrArray(zone,slabVolume,5);
	
		MakeSlabPtrArray(zone,slabVartheta,35);//vartheta_d...
		

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

			Mx_P=*(domainMomentumx+indexP);
			My_P=*(domainMomentumy+indexP);
			Var_P=*(domainVartheta+indexP);
			vol_P=*(domainVolume+indexP);

			if(i==(NI-2))
			{
				Mx_E=*(slabMomentumx[0]+(NJ-1)*k+j);
				My_E=*(slabMomentumy[0]+(NJ-1)*k+j);
				Var_E=*(slabVartheta[0]+(NJ-1)*k+j);
				vol_E=*(slabVolume[0]+(NJ-1)*k+j);
			}
			else
			{
				Mx_E=*(domainMomentumx+indexE);
				My_E=*(domainMomentumy+indexE);
				Var_E=*(domainVartheta+indexE);
				vol_E=*(domainVolume+indexE);
			}

			if(i==0)
			{
				Mx_W =*(slabMomentumx[1]+(NJ-1)*k+j);
				My_W =*(slabMomentumy[1]+(NJ-1)*k+j);
				Var_W =*(slabVartheta[1]+(NJ-1)*k+j);
				vol_W =*(slabVolume[1]+(NJ-1)*k+j);
			}
			else
			{
				Mx_W=*(domainMomentumx+indexW);
				My_W=*(domainMomentumy+indexW);
				Var_W=*(domainVartheta+indexW);
				vol_W=*(domainVolume+indexW);
			}

			if(j==(NJ-2))
			{
				Mx_N=*(slabMomentumx[2]+(NI-1)*k+i);
				My_N=*(slabMomentumy[2]+(NI-1)*k+i);
				Var_N=*(slabVartheta[2]+(NI-1)*k+i);
				vol_N=*(slabVolume[2]+(NI-1)*k+i);
			}
			else
			{
				Mx_N=*(domainMomentumx+indexN);
				My_N=*(domainMomentumy+indexN);
				Var_N=*(domainVartheta+indexN);
				vol_N=*(domainVolume+indexN);
			}

			if(j==0)
			{
				Mx_S=*(slabMomentumx[3]+(NI-1)*k+i);
				My_S=*(slabMomentumy[3]+(NI-1)*k+i);
				Var_S=*(slabVartheta[3]+(NI-1)*k+i);
				vol_S=*(slabVolume[3]+(NI-1)*k+i);
			}
			else
			{
				Mx_S=*(domainMomentumx+indexS);
				My_S=*(domainMomentumy+indexS);
				Var_S=*(domainVartheta+indexS);
				vol_S=*(domainVolume+indexS);
			}

			if(k==(NK-2))
			{
				Mx_T=*(slabMomentumx[4]+(NI-1)*j+i);
				My_T=*(slabMomentumy[4]+(NI-1)*j+i);
				Var_T=*(slabVartheta[4]+(NI-1)*j+i);
				vol_T=*(slabVolume[4]+(NI-1)*j+i);
			}
			else
			{
				Mx_T=*(domainMomentumx+indexT);
				My_T=*(domainMomentumy+indexT);
				Var_T=*(domainVartheta+indexT);
				vol_T=*(domainVolume+indexT);
			}

			if(k==0)
			{
				Mx_B=*(slabMomentumx[5]+(NI-1)*j+i);
				My_B=*(slabMomentumy[5]+(NI-1)*j+i);
				Var_B=*(slabVartheta[5]+(NI-1)*j+i);
				vol_B=*(slabVolume[5]+(NI-1)*j+i);
			}
			else
			{
				Mx_B=*(domainMomentumx+indexB);
				My_B=*(domainMomentumy+indexB);
				Var_B=*(domainVartheta+indexB);
				vol_B=*(domainVolume+indexB);
				
			}

			
			Mxe=((Mx_E*vol_P+Mx_P*vol_E)/(vol_P+vol_E));
			Mxw=((Mx_W*vol_P+Mx_P*vol_W)/(vol_P+vol_W));
			Mxs=((Mx_S*vol_P+Mx_P*vol_S)/(vol_P+vol_S));
			Mxn=((Mx_N*vol_P+Mx_P*vol_N)/(vol_P+vol_N));
			Mxt=((Mx_T*vol_P+Mx_P*vol_T)/(vol_P+vol_T));
			Mxb=((Mx_B*vol_P+Mx_P*vol_B)/(vol_P+vol_B));

			Mye=((My_E*vol_P+My_P*vol_E)/(vol_P+vol_E));
			Myw=((My_W*vol_P+My_P*vol_W)/(vol_P+vol_W));
			Mys=((My_S*vol_P+My_P*vol_S)/(vol_P+vol_S));
			Myn=((My_N*vol_P+My_P*vol_N)/(vol_P+vol_N));
			Myt=((My_T*vol_P+My_P*vol_T)/(vol_P+vol_T));
			Myb=((My_B*vol_P+My_P*vol_B)/(vol_P+vol_B));

			Vare=((Var_E*vol_P+Var_P*vol_E)/(vol_P+vol_E));//No problem with this as there is no volfrac in the denominator
			Varw=((Var_W*vol_P+Var_P*vol_W)/(vol_P+vol_W));
			Vars=((Var_S*vol_P+Var_P*vol_S)/(vol_P+vol_S));
			Varn=((Var_N*vol_P+Var_P*vol_N)/(vol_P+vol_N));
			Vart=((Var_T*vol_P+Var_P*vol_T)/(vol_P+vol_T));
			Varb=((Var_B*vol_P+Var_P*vol_B)/(vol_P+vol_B));
//#######################################################################		
            //Computation of velocities of the disperse phase                                                               
			ue=Mxe/(*(rhoPhase+indexP*flagrho)*Vare);
			ve=Mye/(*(rhoPhase+indexP*flagrho)*Vare);
				
			uw=Mxw/(*(rhoPhase+indexP*flagrho)*Varw);
			vw=Myw/(*(rhoPhase+indexP*flagrho)*Varw);
			
			us=Mxs/(*(rhoPhase+indexP*flagrho)*Vars);
			vs=Mys/(*(rhoPhase+indexP*flagrho)*Vars);
				
			un=Mxn/(*(rhoPhase+indexP*flagrho)*Varn);
			vn=Myn/(*(rhoPhase+indexP*flagrho)*Varn);
		
			
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
			*(F+face_index+5) = 0.0;
			//printf("\n Fw-%e$Fe-%e$Fs-%e$Fn-%e$",*(F+face_index),*(F+face_index+1),*(F+face_index+2),*(F+face_index+3));
			indexP++;
			face_index = face_index + 6;
			area_index = area_index + 18;		
		}/*End of looping over NI,NJ,NK*///Finally we will be having the disperse mass flux in zone->Flowsol->MassfluxDisperse;
	}/*End of looping over the zones */
	return 0;
}

/**
This function computes the volume flux for disperse phase in 3d geometry
*/
int volumefluxMPFV(Base *base)
{
	int nZones,numZones;
	int k,j,i;
	int indexB,indexT,indexE,indexW,indexS,indexN,indexP;
	int K_Increment,J_Increment,I_Increment,face_index,area_index;
	double *domainVolume,*domainMomentumx,*domainMomentumy,*domainMomentumz,*domainRhoDisperse,*domainVartheta;
	double *slabVolume[6],*slabMomentumx[6],*slabMomentumy[6],*slabMomentumz[6],*slabRhoDisperse[6],*slabVartheta[6];	
	double *SurfArea,*F;
	Zone* zone;
	nZones=(base->Number_of_zones);
	
	int NI,NJ,NK,*ptrSize;
	double uw,vw,ww,ue,ve,we,un,vn,wn,us,vs,ws,ut,vt,wt,ub,vb,wb;
	double u_P,u_E,u_W,u_S,u_N,u_T,u_B,v_P,v_E,v_W,v_S,v_N,v_T,v_B,w_P,w_E,w_W,w_S,w_N,w_T,w_B;
	double vol_P,vol_E,vol_W,vol_N,vol_S,vol_T,vol_B;
	
	for(numZones=0;numZones<nZones;numZones++)
	{
		zone=(base->CgnsZone+numZones);
		
		ptrSize=zone->Size;

		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);

		F = zone->PrevFlowsol->MassfluxDisperse;

		SurfArea=zone->ZGeom->SurfaceArea;

		domainMomentumx=MakeDomainArray(zone,62);
		domainMomentumy=MakeDomainArray(zone,63);
		domainMomentumz=MakeDomainArray(zone,64);
		domainVolume=MakeDomainArray(zone,5);
		domainRhoDisperse=MakeDomainArray(zone,66);//RhoDisperse
		domainVartheta=MakeDomainArray(zone,65);

		MakeSlabPtrArray(zone,slabMomentumx,62);
		MakeSlabPtrArray(zone,slabMomentumy,63);
		MakeSlabPtrArray(zone,slabMomentumz,64);
		MakeSlabPtrArray(zone,slabVolume,5);
		MakeSlabPtrArray(zone,slabRhoDisperse,66);//RhoDisperse
		MakeSlabPtrArray(zone,slabVartheta,65);
		
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
			
			*(zone->Flowsol->U_Disperse+indexP)=*(domainMomentumx+indexP)/((*(domainRhoDisperse+indexP)*(*(domainVartheta+indexP))));
			*(zone->Flowsol->V_Disperse+indexP)=*(domainMomentumy+indexP)/((*(domainRhoDisperse+indexP)*(*(domainVartheta+indexP))));
			*(zone->Flowsol->W_Disperse+indexP)=*(domainMomentumz+indexP)/((*(domainRhoDisperse+indexP)*(*(domainVartheta+indexP))));
			indexP++;
		}
		
		indexP=0;
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

			u_P=*(domainMomentumx+indexP)/((*(domainRhoDisperse+indexP)*(*(domainVartheta+indexP))));
			v_P=*(domainMomentumy+indexP)/((*(domainRhoDisperse+indexP)*(*(domainVartheta+indexP))));
			w_P=*(domainMomentumz+indexP)/((*(domainRhoDisperse+indexP)*(*(domainVartheta+indexP))));
			vol_P=*(domainVolume+indexP);

		if(i==(NI-2))
			{
			u_E=*(slabMomentumx[0]+(NJ-1)*k+j)/((*(slabRhoDisperse[0]+(NJ-1)*k+j)*(*(slabVartheta[0]+(NJ-1)*k+j))));
			v_E=*(slabMomentumy[0]+(NJ-1)*k+j)/((*(slabRhoDisperse[0]+(NJ-1)*k+j)*(*(slabVartheta[0]+(NJ-1)*k+j))));
			w_E=*(slabMomentumz[0]+(NJ-1)*k+j)/((*(slabRhoDisperse[0]+(NJ-1)*k+j)*(*(slabVartheta[0]+(NJ-1)*k+j))));
			vol_E=*(slabVolume[0]+(NJ-1)*k+j);
			}
		else
			{
			u_E=*(domainMomentumx+indexE)/((*(domainRhoDisperse+indexE)*(*(domainVartheta+indexE))));
			v_E=*(domainMomentumy+indexE)/((*(domainRhoDisperse+indexE)*(*(domainVartheta+indexE))));
			w_E=*(domainMomentumz+indexE)/((*(domainRhoDisperse+indexE)*(*(domainVartheta+indexE))));
			vol_E=*(domainVolume+indexE);
			}

		if(i==0)
			{
			u_W =*(slabMomentumx[1]+(NJ-1)*k+j)/((*(slabRhoDisperse[1]+(NJ-1)*k+j)*(*(slabVartheta[1]+(NJ-1)*k+j))));
			v_W =*(slabMomentumy[1]+(NJ-1)*k+j)/((*(slabRhoDisperse[1]+(NJ-1)*k+j)*(*(slabVartheta[1]+(NJ-1)*k+j))));
			w_W =*(slabMomentumz[1]+(NJ-1)*k+j)/((*(slabRhoDisperse[1]+(NJ-1)*k+j)*(*(slabVartheta[1]+(NJ-1)*k+j))));
			vol_W =*(slabVolume[1]+(NJ-1)*k+j);
			}
		else
			{
			u_W=*(domainMomentumx+indexW)/((*(domainRhoDisperse+indexW)*(*(domainVartheta+indexW))));
			v_W=*(domainMomentumy+indexW)/((*(domainRhoDisperse+indexW)*(*(domainVartheta+indexW))));
			w_W=*(domainMomentumz+indexW)/((*(domainRhoDisperse+indexW)*(*(domainVartheta+indexW))));
			vol_W=*(domainVolume+indexW);
			}

		if(j==(NJ-2))
			{
			u_N=*(slabMomentumx[2]+(NI-1)*k+i)/((*(slabRhoDisperse[2]+(NI-1)*k+i)*(*(slabVartheta[2]+(NI-1)*k+i))));
			v_N=*(slabMomentumy[2]+(NI-1)*k+i)/((*(slabRhoDisperse[2]+(NI-1)*k+i)*(*(slabVartheta[2]+(NI-1)*k+i))));
			w_N=*(slabMomentumz[2]+(NI-1)*k+i)/((*(slabRhoDisperse[2]+(NI-1)*k+i)*(*(slabVartheta[2]+(NI-1)*k+i))));
			vol_N=*(slabVolume[2]+(NI-1)*k+i);
			}
		else
			{
			u_N=*(domainMomentumx+indexN)/((*(domainRhoDisperse+indexN)*(*(domainVartheta+indexN))));
			v_N=*(domainMomentumy+indexN)/((*(domainRhoDisperse+indexN)*(*(domainVartheta+indexN))));
			w_N=*(domainMomentumz+indexN)/((*(domainRhoDisperse+indexN)*(*(domainVartheta+indexN))));
			vol_N=*(domainVolume+indexN);
			}

		if(j==0)
			{
			u_S=*(slabMomentumx[3]+(NI-1)*k+i)/((*(slabRhoDisperse[3]+(NI-1)*k+i)*(*(slabVartheta[3]+(NI-1)*k+i))));
			v_S=*(slabMomentumy[3]+(NI-1)*k+i)/((*(slabRhoDisperse[3]+(NI-1)*k+i)*(*(slabVartheta[3]+(NI-1)*k+i))));
			w_S=*(slabMomentumz[3]+(NI-1)*k+i)/((*(slabRhoDisperse[3]+(NI-1)*k+i)*(*(slabVartheta[3]+(NI-1)*k+i))));
			vol_S=*(slabVolume[3]+(NI-1)*k+i);
			}
		else
			{
			u_S=*(domainMomentumx+indexS)/((*(domainRhoDisperse+indexS)*(*(domainVartheta+indexS))));
			v_S=*(domainMomentumy+indexS)/((*(domainRhoDisperse+indexS)*(*(domainVartheta+indexS))));
			w_S=*(domainMomentumz+indexS)/((*(domainRhoDisperse+indexS)*(*(domainVartheta+indexS))));
			vol_S=*(domainVolume+indexS);
			}

		if(k==(NK-2))
			{
			u_T=*(slabMomentumx[4]+(NI-1)*j+i)/((*(slabRhoDisperse[4]+(NI-1)*j+i)*(*(slabVartheta[4]+(NI-1)*j+i))));
			v_T=*(slabMomentumy[4]+(NI-1)*j+i)/((*(slabRhoDisperse[4]+(NI-1)*j+i)*(*(slabVartheta[4]+(NI-1)*j+i))));
			w_T=*(slabMomentumz[4]+(NI-1)*j+i)/((*(slabRhoDisperse[4]+(NI-1)*j+i)*(*(slabVartheta[4]+(NI-1)*j+i))));
			vol_T=*(slabVolume[4]+(NI-1)*j+i);
			}
		else
			{
			u_T=*(domainMomentumx+indexT)/((*(domainRhoDisperse+indexT)*(*(domainVartheta+indexT))));
			v_T=*(domainMomentumy+indexT)/((*(domainRhoDisperse+indexT)*(*(domainVartheta+indexT))));
			w_T=*(domainMomentumz+indexT)/((*(domainRhoDisperse+indexT)*(*(domainVartheta+indexT))));
			vol_T=*(domainVolume+indexT);
			}

		if(k==0)
			{
			u_B=*(slabMomentumx[5]+(NI-1)*j+i)/((*(slabRhoDisperse[5]+(NI-1)*j+i)*(*(slabVartheta[5]+(NI-1)*j+i))));
			v_B=*(slabMomentumy[5]+(NI-1)*j+i)/((*(slabRhoDisperse[5]+(NI-1)*j+i)*(*(slabVartheta[5]+(NI-1)*j+i))));
			w_B=*(slabMomentumz[5]+(NI-1)*j+i)/((*(slabRhoDisperse[5]+(NI-1)*j+i)*(*(slabVartheta[5]+(NI-1)*j+i))));
			vol_B=*(slabVolume[5]+(NI-1)*j+i);
			}
		else
			{
			u_B=*(domainMomentumx+indexB)/((*(domainRhoDisperse+indexB)*(*(domainVartheta+indexB))));
			v_B=*(domainMomentumy+indexB)/((*(domainRhoDisperse+indexB)*(*(domainVartheta+indexB))));
			w_B=*(domainMomentumz+indexB)/((*(domainRhoDisperse+indexB)*(*(domainVartheta+indexB))));
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
		
			indexP++;
			face_index = face_index + 6;
			area_index = area_index + 18;	

			
		indexP++;
		}/*End of looping over NI,NJ,NK*/
	}/*End of looping over the zones */
	return 0;
}


