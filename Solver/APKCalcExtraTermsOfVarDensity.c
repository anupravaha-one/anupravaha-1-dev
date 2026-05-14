#include "APKCalcExtraTermsOfVarDensity.h"

/**
In variable density case, the pressure Poisson equation has some extra terms which add up to diagonal and to right hand
side of the equation. These extra terms are calculated by this function which are used in UpdatePressureVarDensity(). It
is called only for 3D case.
The function udotdelRho() computes the last term on the right hand side given in Eqn.(5.40) as given in Project
Report No 1, Dec 2005, Section 5.3.1.2, pp71.
*/

double CalcExtraTermsOfVarDensity(Base *base,Zone *zone)
{
	double sumMassFlux,*massFlux,prevRho,currentRho;
	double CONSTANT,volumeP,prevT,prevP,delta_t,*sumU_pdotdelRho;
	int NI,NJ,NK,i,j,k,iFace,index[3],indexP,*ptrSize,flagrho;

	delta_t = base->BaseParameters->TimeStep;
	ptrSize = zone->Size ;

	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);

	flagrho = zone->Flowsol->flagRho;

	sumU_pdotdelRho = (double *)malloc((NI-1)*(NJ-1)*(NK-1)*sizeof(double));
	
	udotdelRho(zone,sumU_pdotdelRho);

	for(k=0;k<NK-1;k++)
	for(j=0;j<NJ-1;j++)
	for(i=0;i<NI-1;i++)
	{
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		
		index[0]=k,index[1]=j,index[2]=i;
		
		currentRho = *(zone->Flowsol->Rho+indexP*flagrho);
		prevRho = *(zone->PrevFlowsol->Rho+indexP*flagrho);
		prevP = *(zone->PrevFlowsol->P+indexP);
		volumeP = *(zone->ZGeom->Volume+indexP);
		if(zone->Flowsol->T!=NULL)
		{
			if(base->flagtypeVariableDensity==1)
			prevT = 273.0;
			else
			prevT = *(zone->PrevFlowsol->T+indexP);
		}
		else
			prevT = 273.0;
		
		massFlux=zone->Flowsol->Massflux+indexP*6;
		
		CONSTANT = userEnterconstant(zone,indexP,base);//See userDefinedVarCoefficients.c function
		// A constant linking density to pressure; constant = MolWt./ Unviversal Gas Constant.
		
		sumMassFlux=0.0;
		for(iFace=0;iFace<6;iFace++)
		{
			sumMassFlux=sumMassFlux +(*(massFlux+iFace));
		}
		
		if(base->flagSolidification==1)
		{
			currentRho = *(zone->Flowsol->Rho+indexP*flagrho);
			prevRho = *(zone->PrevFlowsol->Rho_prev+indexP*flagrho);
			*((zone->grid+0)->GRHS+indexP) = sumMassFlux + (1.0/prevRho)*(*(sumU_pdotdelRho+indexP)) + (1.0/(prevRho*delta_t)) * volumeP*(currentRho - prevRho);
			*(zone->Flowsol->VarDensityCoeff+indexP) = 0.0;
		}
		else
		{
			*((zone->grid+0)->GRHS+indexP) = sumMassFlux + (1.0/prevRho)*(*(sumU_pdotdelRho+indexP)) - (CONSTANT*volumeP/(delta_t*prevRho*prevT))*prevP;
			
			*(zone->Flowsol->VarDensityCoeff+indexP) = -(CONSTANT*volumeP/(delta_t*prevRho*prevT));
		}
	}
	free(sumU_pdotdelRho);
	return(0);
}

/**
The function udotdelRho() computes the last term on the right hand side given in Eqn.(5.40) as given in Project
Report No 1, Dec 2005, Section 5.3.1.2, pp71.
*/

int udotdelRho(Zone *zone,double *sumU_pdotdelRho)
{
	int k,j,i;
	int Index,indexB,indexT,indexE,indexW,indexS,indexN,indexP;
	double *domainVolume,*domainU,*domainV,*domainW;
	double *slabVolume[6],*slabU[6],*slabV[6],*slabW[6];
	double *slabRho[6],*domainRho;
	double *SurfArea;
	int NI,NJ,NK,*ptrSize;
	double u_P,v_P,w_P;
	double vol_P,vol_E,vol_W,vol_N,vol_S,vol_T,vol_B;
	double rho_P,rho_E,rho_W,rho_N,rho_S,rho_T,rho_B;//Cell centered centered values of Rho
	double rhoe,rhow,rhon,rhos,rhot,rhob;//face centered values of Rho
	double sum = 0.0;
	
	ptrSize=zone->Size;
	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);

	SurfArea=zone->ZGeom->SurfaceArea;

	domainU= MakeDomainArray(zone,1);
	domainV=MakeDomainArray(zone,2);
	domainW=MakeDomainArray(zone,3);
	domainVolume=MakeDomainArray(zone,5);
	domainRho=MakeDomainArray(zone,710);/*Density*/

	MakeSlabPtrArray(zone,slabU,1);
	MakeSlabPtrArray(zone,slabV,2);
	MakeSlabPtrArray(zone,slabW,3);
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

		u_P=*(domainU+indexP);
		v_P=*(domainV+indexP);
		w_P=*(domainW+indexP);
		rho_P=*(domainRho+indexP);
		vol_P=*(domainVolume+indexP);

		if(i==(NI-2))
		{
			rho_E = *(slabRho[0]+(NJ-1)*k+j);
			vol_E = *(slabVolume[0]+(NJ-1)*k+j);
		}
		else
		{
			rho_E = *(domainRho+indexE);
			vol_E = *(domainVolume+indexE);
		}

		if(i==0)
		{
			rho_W =*(slabRho[1]+(NJ-1)*k+j);
			vol_W =*(slabVolume[1]+(NJ-1)*k+j);
		}
		else
		{
			rho_W=*(domainRho+indexW);
			vol_W=*(domainVolume+indexW);
		}

		if(j==(NJ-2))
		{
			rho_N=*(slabRho[2]+(NI-1)*k+i);
			vol_N=*(slabVolume[2]+(NI-1)*k+i);
		}
		else
		{
			rho_N=*(domainRho+indexN);
			vol_N=*(domainVolume+indexN);
		}

		if(j==0)
		{
			rho_S=*(slabRho[3]+(NI-1)*k+i);
			vol_S=*(slabVolume[3]+(NI-1)*k+i);
		}
		else
		{
			rho_S=*(domainRho+indexS);
			vol_S=*(domainVolume+indexS);
		}

		if(k==(NK-2))
		{
			rho_T=*(slabRho[4]+(NI-1)*j+i);
			vol_T=*(slabVolume[4]+(NI-1)*j+i);
		}
		else
		{
			rho_T=*(domainRho+indexT);
			vol_T=*(domainVolume+indexT);
		}

		if(k==0)
		{
			rho_B=*(slabRho[5]+(NI-1)*j+i);
			vol_B=*(slabVolume[5]+(NI-1)*j+i);
		}
		else
		{
			rho_B=*(domainRho+indexB);
			vol_B=*(domainVolume+indexB);
		}

		rhoe =((rho_E*vol_P+rho_P*vol_E)/(vol_P+vol_E));
		rhow =((rho_W*vol_P+rho_P*vol_W)/(vol_P+vol_W));
		rhon=((rho_N*vol_P+rho_P*vol_N)/(vol_P+vol_N));
		rhos=((rho_S*vol_P+rho_P*vol_S)/(vol_P+vol_S));
		rhot=((rho_T*vol_P+rho_P*vol_T)/(vol_P+vol_T));
		rhob=((rho_B*vol_P+rho_P*vol_B)/(vol_P+vol_B));

		Index = (NJ-1)*(NI-1)*k+(NI-1)*j+i;
		sum = 0.0;
		  // west face
			sum+= rhow*(u_P*(*(SurfArea+Index*6*3+3*0+0)) + v_P*(*(SurfArea+Index*6*3+3*0+1)) + w_P*(*(SurfArea+Index*6*3+3*0+2)));

		  // east face
			sum+= rhoe*(u_P*(*(SurfArea+Index*6*3+3*1+0)) + v_P*(*(SurfArea+Index*6*3+3*1+1)) + w_P*(*(SurfArea+Index*6*3+3*1+2)));

		  // south face
			sum+= rhos*(u_P*(*(SurfArea+Index*6*3+3*2+0)) + v_P*(*(SurfArea+Index*6*3+3*2+1)) + w_P*(*(SurfArea+Index*6*3+3*2+2)));

		  // north face
			sum+= rhon*(u_P*(*(SurfArea+Index*6*3+3*3+0)) + v_P*(*(SurfArea+Index*6*3+3*3+1)) + w_P*(*(SurfArea+Index*6*3+3*3+2)));

		  // bottom face
			sum+= rhob*(u_P*(*(SurfArea+Index*6*3+3*4+0)) + v_P*(*(SurfArea+Index*6*3+3*4+1)) + w_P*(*(SurfArea+Index*6*3+3*4+2)));

		  // top face
			sum+= rhot*(u_P*(*(SurfArea+Index*6*3+3*5+0)) + v_P*(*(SurfArea+Index*6*3+3*5+1)) + w_P*(*(SurfArea+Index*6*3+3*5+2)));
		
		*(sumU_pdotdelRho+indexP) = sum;
	}/*End of looping over NI,NJ,NK*/
	return 0;
}

