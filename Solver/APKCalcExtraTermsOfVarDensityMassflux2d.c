#include "APKCalcExtraTermsOfVarDensityMassflux2d.h"

/**
In variable density case, the pressure Poisson equation has some extra terms which add up to diagonal and to right hand
side of the equation. These extra terms are calculated by this function which are used in UpdatePressure2dVarDensityMassFlux(). It
is called only for 2D case with Two phase flow boiling and combustion problems with variable density.
*/

double CalcExtraTermsOfVarDensityMassflux2d(Base *base,Zone *zone)
{
	double sumMassFlux,*massFlux,prevRho,currentRho,currentP;
	double CONSTANT,volumeP,prevT,prevP;
	int NI,NJ,NK,i,j,k,index[3],*ptrSize,flagrho,indexP,iFace;
	double delta_t,*domainRho,*slabRho[6],*newMassFlux;
	
	delta_t = base->BaseParameters->TimeStep;
	ptrSize = zone->Size ;

	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
	
	flagrho = zone->Flowsol->flagRho;
	
	newMassFlux = (double *)malloc(6*(NI-1)*(NJ-1)*(NK-1)*sizeof(double));
	CalcMassFlux(newMassFlux,zone);
        domainRho=MakeDomainArray(zone,710);/*Density*/
        MakeSlabPtrArray(zone,slabRho,710);//Density


	for(k=0;k<NK-1;k++)
	for(j=0;j<NJ-1;j++)
	for(i=0;i<NI-1;i++)
	{
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		

		currentRho = *(zone->Flowsol->Rho+indexP*flagrho);
		prevRho = *(zone->PrevFlowsol->Rho+indexP*flagrho);
		currentP = *(zone->Flowsol->P+indexP);//Contains pressure of current cell
		prevP = *(zone->PrevFlowsol->P+indexP);//Contains Prevpressure of current cell
		volumeP = *(zone->ZGeom->Volume+indexP);
		if(zone->Flowsol->T!=NULL)
		{
			if(base->flagtypeVariableDensity==1)
			prevT = 273.0;
			else
			prevT = *(zone->PrevFlowsol->T+indexP);
		}
		else
		prevT = 1.0;
		
		massFlux=newMassFlux+indexP*6;
		
		CONSTANT = userEnterconstant(zone,indexP,base);//See userDefinedVarCoefficients.c function
		// A constant linking density to pressure; constant = MolWt./ Unviversal Gas Constant.
		
		sumMassFlux=0.0;
		for(iFace=0;iFace<4;iFace++)
		{
			sumMassFlux=sumMassFlux +(*(massFlux+iFace));
		}
//		sumMassFlux=(*(massFlux+0))+(*(massFlux+1))+(*(massFlux+2))+(*(massFlux+3));
		
		*((zone->grid+0)->GRHS+indexP) = sumMassFlux - (CONSTANT*volumeP/(delta_t*prevT))*prevP;
		
		*(zone->Flowsol->VarDensityCoeff+indexP) = -(CONSTANT*volumeP/(delta_t*prevT));
	}
	
	free(newMassFlux);
	return(0);
}


int CalcMassFlux(double *newMassflux, Zone *zone)
{
	int k,j,i,face;
	int Index,indexB,indexT,indexE,indexW,indexS,indexN,indexP;

	double *domainVolume,*domainUmass,*domainVmass,*domainWmass;
	double *slabVolume[6],*slabUmass[6],*slabVmass[6],*slabWmass[6];
	double *slabRho[6],*domainRho;

	double *SurfArea,*F,*rho;
	int NI,NJ,NK,*ptrSize;
	
	double u_P,u_E,u_W,u_S,u_N,u_T,u_B,v_P,v_E,v_W,v_S,v_N,v_T,v_B,w_P,w_E,w_W,w_S,w_N,w_T,w_B;
	double vol_P,vol_E,vol_W,vol_N,vol_S,vol_T,vol_B;
	double rho_P,rho_E,rho_W,rho_N,rho_S,rho_T,rho_B;
	double rhoue,rhove,rhowe,rhouw,rhovw,rhoww,rhoun,rhovn,rhown,rhous,rhovs,rhows,rhout,rhovt,rhowt,rhoub,rhovb,rhowb;
	
	int K_Increment,J_Increment;
	
	ptrSize=zone->Size;
	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
	K_Increment = (NI-1) *(NJ-1);
	J_Increment = NI-1;
	F = newMassflux;

	SurfArea=zone->ZGeom->SurfaceArea;

	rho=zone->Flowsol->Rho;
	
	if(zone->Flowsol->w!=NULL)
	{
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
	
		indexP=0;
		for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{
//			indexP=(NJ-1)*(NI-1)*k+(NI-1)*j+i;
// 			indexE=(NJ-1)*(NI-1)*k+(NI-1)*j+(i+1);
// 			indexW=(NJ-1)*(NI-1)*k+(NI-1)*j+(i-1);
// 			indexN=(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i;
// 			indexS=(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i;
// 			indexT=(NJ-1)*(NI-1)*(k+1)+(NI-1)*j+i;
// 			indexB=(NJ-1)*(NI-1)*(k-1)+(NI-1)*j+i;
			
			indexE=indexP+1;
			indexW=indexP-1;
			indexN=indexP+J_Increment;
			indexS=indexP-J_Increment;
			indexT=indexP+K_Increment;
			indexB=indexP-K_Increment;
	
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
				indexP++;
		}/*End of looping over NI,NJ,NK*/
	}
	else
	{
		domainUmass=MakeDomainArray(zone,100);
		domainVmass=MakeDomainArray(zone,200);
		domainVolume=MakeDomainArray(zone,5);
		domainRho=MakeDomainArray(zone,710);/*Density*/
	
		MakeSlabPtrArray(zone,slabUmass,100);
		MakeSlabPtrArray(zone,slabVmass,200);
		MakeSlabPtrArray(zone,slabVolume,5);
		MakeSlabPtrArray(zone,slabRho,710);
		indexP=0;
		for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{
// 			indexP=(NJ-1)*(NI-1)*k+(NI-1)*j+i;
// 			indexE=(NJ-1)*(NI-1)*k+(NI-1)*j+(i+1);
// 			indexW=(NJ-1)*(NI-1)*k+(NI-1)*j+(i-1);
// 			indexN=(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i;
// 			indexS=(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i;
// 			indexT=(NJ-1)*(NI-1)*(k+1)+(NI-1)*j+i;
// 			indexB=(NJ-1)*(NI-1)*(k-1)+(NI-1)*j+i;
			
			indexE=indexP+1;
			indexW=indexP-1;
			indexN=indexP+J_Increment;
			indexS=indexP-J_Increment;
			indexT=indexP+K_Increment;
			indexB=indexP-K_Increment;
	
			u_P=*(domainUmass+indexP);
			v_P=*(domainVmass+indexP);
			rho_P=*(domainRho+indexP);
			vol_P=*(domainVolume+indexP);
	
			if(i==(NI-2))
			{
				u_E=*(slabUmass[0]+(NJ-1)*k+j);
				v_E=*(slabVmass[0]+(NJ-1)*k+j);
				rho_E=*(slabRho[0]+(NJ-1)*k+j);
				vol_E=*(slabVolume[0]+(NJ-1)*k+j);
			}
			else
			{
				u_E=*(domainUmass+indexE);
				v_E=*(domainVmass+indexE);
				rho_E=*(domainRho+indexE);
				vol_E=*(domainVolume+indexE);
			}
	
			if(i==0)
			{
				u_W =*(slabUmass[1]+(NJ-1)*k+j);
				v_W =*(slabVmass[1]+(NJ-1)*k+j);
				rho_W =*(slabRho[1]+(NJ-1)*k+j);
				vol_W =*(slabVolume[1]+(NJ-1)*k+j);
			}
			else
			{
				u_W=*(domainUmass+indexW);
				v_W=*(domainVmass+indexW);
				rho_W=*(domainRho+indexW);
				vol_W=*(domainVolume+indexW);
			}
	
			if(j==(NJ-2))
			{
				u_N=*(slabUmass[2]+(NI-1)*k+i);
				v_N=*(slabVmass[2]+(NI-1)*k+i);
				rho_N=*(slabRho[2]+(NI-1)*k+i);
				vol_N=*(slabVolume[2]+(NI-1)*k+i);
			}
			else
			{
				u_N=*(domainUmass+indexN);
				v_N=*(domainVmass+indexN);
				rho_N=*(domainRho+indexN);
				vol_N=*(domainVolume+indexN);
			}
	
			if(j==0)
			{
				u_S=*(slabUmass[3]+(NI-1)*k+i);
				v_S=*(slabVmass[3]+(NI-1)*k+i);
				rho_S=*(slabRho[3]+(NI-1)*k+i);
				vol_S=*(slabVolume[3]+(NI-1)*k+i);
			}
			else
			{
				u_S=*(domainUmass+indexS);
				v_S=*(domainVmass+indexS);
				rho_S=*(domainRho+indexS);
				vol_S=*(domainVolume+indexS);
			}
	
			rhoue=((rho_E*u_E*vol_P+rho_P*u_P*vol_E)/(vol_P+vol_E));
			rhove=((rho_E*v_E*vol_P+rho_P*v_P*vol_E)/(vol_P+vol_E));
	
			rhouw=((rho_W*u_W*vol_P+rho_P*u_P*vol_W)/(vol_P+vol_W));
			rhovw=((rho_W*v_W*vol_P+rho_P*v_P*vol_W)/(vol_P+vol_W));
	
			rhoun=((rho_N*u_N*vol_P+rho_P*u_P*vol_N)/(vol_P+vol_N));
			rhovn=((rho_N*v_N*vol_P+rho_P*v_P*vol_N)/(vol_P+vol_N));
	
			rhous=((rho_S*u_S*vol_P+rho_P*u_P*vol_S)/(vol_P+vol_S));
			rhovs=((rho_S*v_S*vol_P+rho_P*v_P*vol_S)/(vol_P+vol_S));

	
			Index=indexP*18;
			face = 0;  // west face
			*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = (rhouw*(*(SurfArea+Index+0))+rhovw*(*(SurfArea+Index+1)));
				//fprintf(fp1,"west=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
			face = 1;  // east face
			*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = (rhoue*(*(SurfArea+Index+3))+rhove*(*(SurfArea+Index+4)));
				//fprintf(fp1,"east=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
			face = 2;  // south face
			*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = (rhous*(*(SurfArea+Index+6))+rhovs*(*(SurfArea+Index+7)));
				//fprintf(fp1,"south=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
			face = 3;  // north face
			*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face) = (rhoun*(*(SurfArea+Index+9))+rhovn*(*(SurfArea+Index+10)));
				//fprintf(fp1,"north=%lf\t",*(F+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face));
			indexP++;
		}/*End of looping over NI,NJ,NK*/
	}

return 0;
}
