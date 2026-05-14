#include "APKLiftMPFV.h"
/**
This function calculates the Lift force for the case of disperse phase flow. The Lift is calculated and stored in the data structure under the source term. This function is called once 
per time step(Lift force added as in Moukalled et al.)
 */
double LiftMPFV(Base *base,int CounterTstep)
{
	double mu,rho,rho1,v1,v2,dia;
	double *ptrVolume;
	int cntZ,NI,NJ,NK,*ptrSize,numZones,i,j,k;
	Zone *zone;
	int flagDen,flagMu,flagDia;//Flags
	int totalcells,index,nophase,nop;/*,flagrhoDisperse=0;*/

    	int P,E,W,N,S,T,B;
	double UE,UW,UN,US,UT,UB,UP,VE,VW,VN,VS,VT,VB,VP,WE,WW,WN,WS,WT,WB,WP;
	double DUx,DUy,DUz,DVx,DVy,DVz,DWx,DWy,DWz,pi[3],*deltainv;

	double Vc_X,Vc_Y,Vc_Z,Vd_X,Vd_Y,Vd_Z,Uc_X,Uc_Y,Uc_Z,Ud_X,Ud_Y,Ud_Z;
	double Udiff[3],Curl[3],C1,C1a,Min,Px,Py,Pz,Lift_X,Lift_Y,Lift_Z;
	
	double *variable,*variable1,*variable2,*slabVariable[6],*slabVariable1[6],*slabVariable2[6],*vartheta,*slabvartheta[6];
	
	
    	
	numZones= base->Number_of_zones;

	nop = base->nphase;
	
	for(nophase=0;nophase<nop;nophase++)
	{	
		flagMu = base->BaseParameters->flagViscoeffPhase[nophase];
		flagDen = base->BaseParameters->flagRhoPhase[nophase];
	}
	
	for(cntZ=0;cntZ<numZones;cntZ++)
	{
		zone = (base->CgnsZone+cntZ);
		
		variable=MakeDomainArray(zone,41);
		MakeSlabPtrArray(zone,slabVariable,41);
	
		variable1=MakeDomainArray(zone,42);
		MakeSlabPtrArray(zone,slabVariable1,42);

		vartheta=MakeDomainArray(zone,34);
		MakeSlabPtrArray(zone,slabvartheta,34);
		
		if(base->BaseParameters->flag2d==0)
		{
			variable2=MakeDomainArray(zone,43);
			MakeSlabPtrArray(zone,slabVariable2,43);
		}
		deltainv = zone->ZGeom->DeltaInverse;
		ptrSize = zone->Size;
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		totalcells = (NI-1)*(NJ-1)*(NK-1);
		
		flagDia=base->BaseParameters->flagVariableDia;
		ptrVolume= zone->ZGeom->Volume;
		index=0;
		
        	for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{	
			
		
			v1 = *(zone->Flowsol->VarthetaPhase[0]+index);//volfrac 0
			v2 = *(zone->Flowsol->VarthetaPhase[1]+index);//volfrac 1
			
			mu= *(zone->Flowsol->ViscosityofPhase[0]+index*flagMu);
			
			rho= *(zone->Flowsol->RhoofPhase[0]+index*flagDen);//density 0
		    	rho1= *(zone->Flowsol->RhoofPhase[1]+index*flagDen);//density 1
			
			dia = *(zone->Flowsol->ParticleDia+index*flagDia);
			
			
			Vc_X=(*((zone->Flowsol)->Momentumofphasex[0]+index));
			Vc_Y=(*((zone->Flowsol)->Momentumofphasey[0]+index));
			if((zone->PrevFlowsol->Momentumofphasez[0])!=NULL)
			Vc_Z=(*((zone->Flowsol)->Momentumofphasez[0]+index));
			
			Vd_X=(*((zone->Flowsol)->Momentumofphasex[1]+index));
			Vd_Y=(*((zone->Flowsol)->Momentumofphasey[1]+index));
			if((zone->Flowsol->Momentumofphasez[1])!=NULL)
			Vd_Z=(*((zone->Flowsol)->Momentumofphasez[1]+index));
			
			Uc_X=Vc_X/(v1*rho);
			Uc_Y=Vc_Y/(v1*rho);
			if((zone->Flowsol->Momentumofphasez[0])!=NULL)
			Uc_Z=Vc_Z/(v1*rho);
			
			Ud_X=Vd_X/(v2*rho1);
			Ud_Y=Vd_Y/(v2*rho1);
			if((zone->Flowsol->Momentumofphasez[1])!=NULL)
			Ud_Z=Vd_Z/(v2*rho1);
			
            		Udiff[0]=Uc_X-Ud_X;
			Udiff[1]=Uc_Y-Ud_Y;
			Udiff[2]=Uc_Z-Ud_Z;
			//printf("\nk=%d,j=%d,i=%d",k,j,i);
			//printf("\nUdiff[0]=%lf$Udiff[1]=%lf,Udiff[2]=%lf",Udiff[0],Udiff[1],Udiff[2]);
			if(0.2>v2)
			Min=v2;
			else
			Min=0.2;
			
			C1=0.075*(1-(2.78*Min));
			//printf("\nC1=%lf",C1);
			P=(NI-1)*(NJ-1)*k+(NI-1)*j+i; 
			E=(NI-1)*(NJ-1)*k+(NI-1)*j+(i+1);
			W=(NI-1)*(NJ-1)*k+(NI-1)*j+(i-1);
			N=(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i;
			S=(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i;
			T=(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i;
			B=(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i;
           
			if(base->BaseParameters->flag2d==1)
			{
				UP=*(variable+P)/(rho* *(vartheta+P));
				VP=*(variable1+P)/(rho**(vartheta+P));
				
				if(i==(NI-2))
				{
					UE=*(slabVariable[0]+(NJ-1)*k+j)/(rho* *(slabvartheta[0]+(NJ-1)*k+j));
					VE=*(slabVariable1[0]+(NJ-1)*k+j)/(rho* *(slabvartheta[0]+(NJ-1)*k+j));
					//printf("\nk=%d,j=%d,i=%d",k,j,i);
					//printf("\n*(slabVariable[0]+(NJ-1)*k+j)=%e",*(slabVariable[0]+(NJ-1)*k+j));
				}
				else
				{
					UE=*(variable+E)/(rho**(vartheta+E));
					VE=*(variable1+E)/(rho**(vartheta+E));
				}
				if(i==0)
				{
					UW=*(slabVariable[1]+(NJ-1)*k+j)/(rho**(slabvartheta[1]+(NJ-1)*k+j));
					VW=*(slabVariable1[1]+(NJ-1)*k+j)/(rho**(slabvartheta[1]+(NJ-1)*k+j));
				}
				else
				{
					UW=*(variable+W)/(rho**(vartheta+W));
					VW=*(variable1+W)/(rho**(vartheta+W));
				}
				if(j==(NJ-2))
				{
					UN=*(slabVariable[2]+(NI-1)*k+i)/(rho* *(slabvartheta[2]+(NI-1)*k+i));
					VN=*(slabVariable1[2]+(NI-1)*k+i)/(rho* *(slabvartheta[2]+(NI-1)*k+i));
				}
				else
				{
					UN=*(variable+N)/(rho* *(vartheta+N));
					VN=*(variable1+N)/(rho* *(vartheta+N));
				}
				if(j==0)
				{
					US=*(slabVariable[3]+(NI-1)*k+i)/(rho* *(slabvartheta[3]+(NI-1)*k+i));
					VS=*(slabVariable1[3]+(NI-1)*k+i)/(rho* *(slabvartheta[3]+(NI-1)*k+i));
				}
				else
				{
					US=*(variable+S)/(rho* *(vartheta+S));
					VS=*(variable1+S)/(rho* *(vartheta+S));
				}
				if(k==0)
				{
					UB=*(slabVariable[5]+(NI-1)*j+i)/(rho* *(slabvartheta[5]+(NI-1)*j+i));
					VB=*(slabVariable1[5]+(NI-1)*j+i)/(rho* *(slabvartheta[5]+(NI-1)*j+i));
				}
				else
				{
					UB=*(variable+B)/(rho**(vartheta+B));
					VB=*(variable1+B)/(rho**(vartheta+B));
				}
				if(k==(NK-2))
				{
					UT=*(slabVariable[4]+(NI-1)*j+i)/(rho* *(slabvartheta[4]+(NI-1)*j+i));
					VT=*(slabVariable1[4]+(NI-1)*j+i)/(rho* *(slabvartheta[4]+(NI-1)*j+i));
				}
				else
				{
					UT=*(variable+T)/(rho* *(vartheta+T));
					VT=*(variable1+T)/(rho* *(vartheta+T));
				}
				
				//printf("\nUE=%lf$UW=%lf$UN=%lfUS=%lf",UE,UW,UN,US);
				//printf("\nVE=%lf$VW=%lf$VN=%lfVS=%lf",VE,VW,VN,VS);
				pi[0]=UE-UW;
				pi[1]=UN-US;
				DUx=*(deltainv+P*4) * pi[0] + *((deltainv+P*4)+1) * pi[1];
				DUy=*((deltainv+P*4)+2) * pi[0] + *((deltainv+P*4)+3) * pi[1];
				
				pi[0]=VE-VW;
				pi[1]=VN-VS;
				DVx=*(deltainv+P*4) * pi[0] + *((deltainv+P*4)+1) * pi[1];
				DVy=*((deltainv+P*4)+2) * pi[0] + *((deltainv+P*4)+3) * pi[1];
				
				Curl[0]=0;
				Curl[1]=0;
				Curl[2]=DVx-DUy;
				
				Px=Udiff[0]*Curl[0];
				Py=Udiff[1]*Curl[1];
				Pz=Udiff[2]*Curl[2];
			
			Lift_X=C1*rho*v2*((Udiff[1]*Curl[2])-(Udiff[2]*Curl[1]));
			Lift_Y=C1*rho*v2*((Udiff[2]*Curl[0])-(Udiff[0]*Curl[2]));
			Lift_Z=C1*rho*v2*((Udiff[0]*Curl[1])-(Udiff[1]*Curl[0]));//Zero Always as  Curl[0]=0; Curl[1]=0;
			
			*(zone->LiftX+index)=Lift_X;
			*(zone->LiftY+index)=Lift_Y;	
			*(zone->LiftZ+index)=Lift_Z;
			
			//printf("\nLift_X=%e,Lift_Y=%e,Lift_Z=%e",Lift_X,Lift_Y,Lift_Z);
			
			}
			else
			{
				UP=*(variable+P)/(rho**(vartheta+P));
				VP=*(variable1+P)/(rho**(vartheta+P));
				WP=*(variable2+P)/(rho**(vartheta+P));
				if(i==(NI-2))
				{
					UE=*(slabVariable[0]+(NJ-1)*k+j)/(rho**(slabvartheta[0]+(NJ-1)*k+j));
					VE=*(slabVariable1[0]+(NJ-1)*k+j)/(rho**(slabvartheta[0]+(NJ-1)*k+j));
					WE=*(slabVariable2[0]+(NJ-1)*k+j)/(rho**(slabvartheta[0]+(NJ-1)*k+j));
				}
				else
				{
					UE=*(variable+E)/(rho**(vartheta+E));
					VE=*(variable1+E)/(rho**(vartheta+E));
					WE=*(variable2+E)/(rho**(vartheta+E));
				}
				if(i==0)
				{
					UW=*(slabVariable[1]+(NJ-1)*k+j)/(rho**(slabvartheta[1]+(NJ-1)*k+j));
					VW=*(slabVariable1[1]+(NJ-1)*k+j)/(rho**(slabvartheta[1]+(NJ-1)*k+j));
					WW=*(slabVariable2[1]+(NJ-1)*k+j)/(rho**(slabvartheta[1]+(NJ-1)*k+j));
				}
				else
				{
					UW=*(variable+W)/(rho**(vartheta+W));
					VW=*(variable1+W)/(rho**(vartheta+W));
					WW=*(variable2+W)/(rho**(vartheta+W));
				}
				if(j==(NJ-2))
				{
					UN=*(slabVariable[2]+(NI-1)*k+i)/(rho**(slabvartheta[2]+(NI-1)*k+i));
					VN=*(slabVariable1[2]+(NI-1)*k+i)/(rho**(slabvartheta[2]+(NI-1)*k+i));
					WN=*(slabVariable2[2]+(NI-1)*k+i)/(rho**(slabvartheta[2]+(NI-1)*k+i));
				}
				else
				{
					UN=*(variable+N)/(rho**(vartheta+N));
					VN=*(variable1+N)/(rho**(vartheta+N));
					WN=*(variable2+N)/(rho**(vartheta+N));
				}
				if(j==0)
				{
					US=*(slabVariable[3]+(NI-1)*k+i)/(rho**(slabvartheta[3]+(NI-1)*k+i));
					VS=*(slabVariable1[3]+(NI-1)*k+i)/(rho**(slabvartheta[3]+(NI-1)*k+i));
					WS=*(slabVariable2[3]+(NI-1)*k+i)/(rho**(slabvartheta[3]+(NI-1)*k+i));
				}
				else
				{
					US=*(variable+S)/(rho**(vartheta+S));
					VS=*(variable1+S)/(rho**(vartheta+S));
					WS=*(variable2+S)/(rho**(vartheta+S));
				}
				if(k==0)
				{
					UB=*(slabVariable[5]+(NI-1)*j+i)/(rho**(slabvartheta[5]+(NI-1)*j+i));
					VB=*(slabVariable1[5]+(NI-1)*j+i)/(rho**(slabvartheta[5]+(NI-1)*j+i));
					WB=*(slabVariable2[5]+(NI-1)*j+i)/(rho**(slabvartheta[5]+(NI-1)*j+i));
				}
				else
				{
					UB=*(variable+B)/(rho**(vartheta+B));
					VB=*(variable1+B)/(rho**(vartheta+B));
					WB=*(variable2+B)/(rho**(vartheta+B));
				}
				if(k==(NK-2))
				{
					UT=*(slabVariable[4]+(NI-1)*j+i)/(rho**(slabvartheta[4]+(NI-1)*j+i));
					VT=*(slabVariable1[4]+(NI-1)*j+i)/(rho**(slabvartheta[4]+(NI-1)*j+i));
					WT=*(slabVariable2[4]+(NI-1)*j+i)/(rho**(slabvartheta[4]+(NI-1)*j+i));
				}
				else
				{
					UT=*(variable+T)/(rho**(vartheta+T));
					VT=*(variable1+T)/(rho**(vartheta+T));
					WT=*(variable2+T)/(rho**(vartheta+T));
				}
					pi[0]=UE-UW;
					pi[1]=UN-US;
					pi[2]=UT-UB;
			
					DUx=*((deltainv+P*9)+0) * pi[0] + *((deltainv+P*9)+1) * pi[1] + *((deltainv+P*9)+2) * pi[2];
					DUy=*((deltainv+P*9)+3) * pi[0] + *((deltainv+P*9)+4) * pi[1] + *((deltainv+P*9)+5) * pi[2];
					DUz=*((deltainv+P*9)+6) * pi[0] + *((deltainv+P*9)+7) * pi[1] + *((deltainv+P*9)+8) * pi[2];
			
					pi[0]=VE-VW;
					pi[1]=VN-VS;
					pi[2]=VT-VB;
			
					DVx=*((deltainv+P*9)+0) * pi[0] + *((deltainv+P*9)+1) * pi[1] + *((deltainv+P*9)+2) * pi[2];
					DVy=*((deltainv+P*9)+3) * pi[0] + *((deltainv+P*9)+4) * pi[1] + *((deltainv+P*9)+5) * pi[2];
					DVz=*((deltainv+P*9)+6) * pi[0] + *((deltainv+P*9)+7) * pi[1] + *((deltainv+P*9)+8) * pi[2];
			
					pi[0]=WE-WW;
					pi[1]=WN-WS;
					pi[2]=WT-WB;
			
					DWx=*((deltainv+P*9)+0) * pi[0] + *((deltainv+P*9)+1) * pi[1] + *((deltainv+P*9)+2) * pi[2];
					DWy=*((deltainv+P*9)+3) * pi[0] + *((deltainv+P*9)+4) * pi[1] + *((deltainv+P*9)+5) * pi[2];
					DWz=*((deltainv+P*9)+6) * pi[0] + *((deltainv+P*9)+7) * pi[1] + *((deltainv+P*9)+8) * pi[2];
			
					Curl[0]=DWy-DVz;
					Curl[1]=DUz-DWx;
					Curl[2]=DVx-DUy;
			
					Px=Udiff[0]*Curl[0];
					Py=Udiff[1]*Curl[1];
					Pz=Udiff[2]*Curl[2];
			
					Lift_X=C1*rho*v2*((Udiff[1]*Curl[2])-(Udiff[2]*Curl[1]));
					Lift_Y=C1*rho*v2*((Udiff[2]*Curl[0])-(Udiff[0]*Curl[2]));
					Lift_Z=C1*rho*v2*((Udiff[0]*Curl[1])-(Udiff[1]*Curl[0]));
					
					*(zone->LiftX+index)=Lift_X;
					*(zone->LiftY+index)=Lift_Y;	
					*(zone->LiftZ+index)=Lift_Z;
			}
           
		index++;
		}
			
	}
	
	return 0;

}
