#include "APKUpdatePressureArtificial.h"
/**
The UpdatePressureArtificial() does the Gauss Seidel Iterations for modified Pressure Poisson equation for Artificial Compressibility
in 3 dimensional geometries. See Eqn(3.70), Section 3.9.0.1, pp 46, in Project Report No 1, Dec 2005.
*/
/*
double UpdatePressureArtificial(Base *base,Zone *zone)
{
	double residual=0.0,residZone=0.0,omega=1.2,*pressureDomain,*prevpressureDomain,sumMassFlux,sumDiff,diffdia;
	double *massFlux,pressureDiff[6];
	int NI,NJ,NK,i,j,k,iFace,index[3],indexP,*ptrSize,flagrho;
	double *slabVolume[6],*slabPressure[6],volume;
	double *beta,*volumeDomain,*diffusionCoeffdomain,rho,betatau=1.0;
	double delta_t;
	
	delta_t = base->BaseParameters->TimeStep;
	pressureDomain = zone->Flowsol->P;
	prevpressureDomain = zone->PrevFlowsol->P;
	volumeDomain = zone->ZGeom->Volume;
	ptrSize = zone->Size ;
	beta=zone->ZGeom->Beta;
	
	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
	
	MakeSlabPtrArray(zone,slabPressure,4);//Called for making slabs of Pressure
	MakeSlabPtrArray(zone,slabVolume,5);//Called for making slabs of Volumes
	diffusionCoeffdomain = MakeDiffusionCoeffDomain(zone,4);
	flagrho = zone->Flowsol->flagRho;
		
	for(k=0;k<NK-1;k++)
	for(j=0;j<NJ-1;j++)
	for(i=0;i<NI-1;i++)
	{
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		
		volume=*(volumeDomain+indexP);
		
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		
		diffdia = (-delta_t/rho) *Diffusion3dWithConstDiffCoeff(index,pressureDomain,volumeDomain,ptrSize,slabPressure,slabVolume,diffusionCoeffdomain,beta,pressureDiff);
		
		massFlux=zone->Flowsol->Massflux+indexP*6;//assign the pointer here
		
		sumDiff=0.0;
		sumMassFlux=0.0;
		for(iFace=0;iFace<6;iFace++)
		{
			sumDiff=sumDiff+pressureDiff[iFace];
			sumMassFlux=sumMassFlux+(*(massFlux+iFace)); 
		}
		diffdia=diffdia + (volume/(rho*betatau));
		//residual= -sumMassFlux+delta_t*sumDiff+(volume*(*(prevpressureDomain+indexP))/betatau)-(volume*(*(pressureDomain+indexP))/betatau);
		residual= -sumMassFlux+(delta_t/rho)*sumDiff+(volume/(rho*betatau))*(*(prevpressureDomain+indexP) - *(pressureDomain+indexP));
		
		*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*residual/diffdia;
		
		residZone=residZone+residual*residual;
	}
	return(residZone);
}
*/
double UpdatePressureArtificial(Base *base,Zone *zone)
{
//	printf("Entered APKUpdatePressure2d.c\n");
	double residZone=0.0,*pressureDomain,*diffdia,*volumeDomain;
	double *LHS,*RHS,resid,rho,diagonal,volume,betatau=1.0;
	int NI,NJ,NK,i,j,k,index[3],indexP,*ptrSize,flagrho,cntDir;
	double delta_t,*CoeffN[18],V[18],Vg,omega=1.0,*slabPressure[6];
	
	FILE *Pomega;
	
	LHS=(double *)malloc(sizeof(double));
	RHS=(double *)malloc(sizeof(double));
	Pomega = fopen("Settings/PressureRelaxation.dat","r");
	fscanf(Pomega,"%lf", &omega);
	fclose(Pomega);
	
	diffdia = (zone->grid)->DiffusionCoeff_P;
	for(i=0;i<18;i++)
	{
		CoeffN[i] = (zone->grid)->DiffusionCoeff_nb[i];
	}
	delta_t = base->BaseParameters->TimeStep;
	pressureDomain = zone->Flowsol->P;
	ptrSize = zone->Size ;
	
	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
	
	for(cntDir=0;cntDir<6;cntDir++)
	{
		slabPressure[cntDir] = zone->Flowsol->SlabPress[cntDir];
	}
//	MakeSlabPtrArray(zone,slabPressure,4);//Called for making slabs of Pressure
	flagrho = zone->Flowsol->flagRho;

//--------------------------------------------- Boundary cells -----------------------------------(1)--------

	k = 0; j = 0; i = 0;
	indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
	index[0]=k,index[1]=j,index[2]=i;
	
	rho = *(zone->Flowsol->Rho+indexP*flagrho);
	Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
	V[0] = *(slabPressure[1]+(NJ-1)*k+j);
	V[1] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
	V[2] = *(slabPressure[3]+(NI-1)*k+i);
	V[3] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
	V[4] = *(slabPressure[5]+(NI-1)*j+i);
	V[5] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
	V[6] = *(slabPressure[1]+(NJ-1)*k+j+1);
	V[7] = 0.0;
	V[8] = *(slabPressure[1]+(NJ-1)*(k+1)+j);
	V[9] = 0.0;
	V[10] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1);
	V[11] = *(slabPressure[3]+(NI-1)*k+i+1);
	V[12] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1);
	V[13] = *(slabPressure[5]+(NI-1)*j+i+1);
	V[14] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i);
	V[15] = *(slabPressure[3]+(NI-1)*(k+1)+i);
	V[16] = *(slabPressure[5]+(NI-1)*(j+1)+i);
	V[17] = 0.0;
	
	volume = *(volumeDomain+indexP);
	
	MatrixMultCorner(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
	
	resid = *((zone->grid)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS)-(volume/(rho*betatau))*Vg;
	
	diagonal = (delta_t/rho)*(*(diffdia+indexP)) + volume/(rho*betatau);

	*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/diagonal;
	
	residZone=residZone+resid*resid;
//--------------------------------------------------------------------------------------------------(2)---------
	k = 0; j = 0;
	
	for(i=1;i<NI-2;i++)
	{
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		V[0] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
		V[1] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
		V[2] = *(slabPressure[3]+(NI-1)*k+i);
		V[3] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
		V[4] = *(slabPressure[5]+(NI-1)*j+i);
		V[5] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
		V[6] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i-1);
		V[7] = *(slabPressure[3]+(NI-1)*k+i-1);
		V[8] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i-1);
		V[9] = *(slabPressure[5]+(NI-1)*j+i-1);
		V[10] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1);
		V[11] = *(slabPressure[3]+(NI-1)*k+i+1);
		V[12] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1);
		V[13] = *(slabPressure[5]+(NI-1)*j+i+1);
		V[14] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i);
		V[15] = *(slabPressure[3]+(NI-1)*(k+1)+i);
		V[16] = *(slabPressure[5]+(NI-1)*(j+1)+i);
		V[17] = 0.0;
		
		volume = *(volumeDomain+indexP);
		
		MatrixMultCorner(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
		
		resid = *((zone->grid)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS)-(volume/(rho*betatau))*Vg;
		
		diagonal = (delta_t/rho)*(*(diffdia+indexP)) + volume/(rho*betatau);
	
		*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/diagonal;
		
		residZone=residZone+resid*resid;
	}
//------------------------------------------------------------------------------------------------------(3)----------
	k = 0; j = NJ-2;
	for(i=1;i<NI-2;i++)
	{
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		V[0] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
		V[1] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
		V[2] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
		V[3] = *(slabPressure[2]+(NI-1)*k+i);
		V[4] = *(slabPressure[5]+(NI-1)*j+i);
		V[5] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
		V[6] = *(slabPressure[2]+(NI-1)*k+i-1);
		V[7] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i-1);
		V[8] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i-1);
		V[9] = *(slabPressure[5]+(NI-1)*j+i-1);
		V[10] = *(slabPressure[2]+(NI-1)*k+i+1);
		V[11] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i+1);
		V[12] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1);
		V[13] = *(slabPressure[5]+(NI-1)*j+i+1);
		V[14] = *(slabPressure[2]+(NI-1)*(k+1)+i);
		V[15] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j-1)+i);
		V[16] = 0.0;
		V[17] = *(slabPressure[5]+(NI-1)*(j-1)+i);
		
		volume = *(volumeDomain+indexP);
		
		MatrixMultCorner(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
		
		resid = *((zone->grid)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS)-(volume/(rho*betatau))*Vg;
		
		diagonal = (delta_t/rho)*(*(diffdia+indexP)) + volume/(rho*betatau);
	
		*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/diagonal;
		
		residZone=residZone+resid*resid;
	}
//----------------------------------------------------------------------------------------------------(4)-----------
	j = NJ-2;
	for(k=1;k<NK-2;k++)
	{
		for(i=1;i<NI-2;i++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			V[0] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			V[1] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			V[2] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			V[3] = *(slabPressure[2]+(NI-1)*k+i);
			V[4] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
			V[5] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
			V[6] = *(slabPressure[2]+(NI-1)*k+i-1);
			V[7] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i-1);
			V[8] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i-1);
			V[9] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i-1);
			V[10] = *(slabPressure[2]+(NI-1)*k+i+1);
			V[11] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i+1);
			V[12] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1);
			V[13] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i+1);
			V[14] = *(slabPressure[2]+(NI-1)*(k+1)+i);
			V[15] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j-1)+i);
			V[16] = *(slabPressure[2]+(NI-1)*(k-1)+i);
			V[17] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+i);
			
			volume = *(volumeDomain+indexP);
			
			MatrixMultCorner(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
			
			resid = *((zone->grid)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS)-(volume/(rho*betatau))*Vg;
			
			diagonal = (delta_t/rho)*(*(diffdia+indexP)) + volume/(rho*betatau);
		
			*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/diagonal;
			
			residZone=residZone+resid*resid;
		}
	}
//------------------------------------------------------------------------------------------------------(5)----------
	k = 0; j = NJ-2; i = NI-2; 
	indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
	index[0]=k,index[1]=j,index[2]=i;
	
	rho = *(zone->Flowsol->Rho+indexP*flagrho);
	Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
	V[0] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
	V[1] = *(slabPressure[0]+(NJ-1)*k+j);
	V[2] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
	V[3] = *(slabPressure[2]+(NI-1)*k+i);
	V[4] = *(slabPressure[5]+(NI-1)*j+i);
	V[5] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
	V[6] = *(slabPressure[2]+(NI-1)*k+i-1);
	V[7] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i-1);
	V[8] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i-1);
	V[9] = *(slabPressure[5]+(NI-1)*j+i-1);
	V[10] = 0.0;
	V[11] = *(slabPressure[0]+(NJ-1)*k+j-1);
	V[12] = *(slabPressure[0]+(NJ-1)*(k+1)+j);
	V[13] = 0.0;
	V[14] = *(slabPressure[2]+(NI-1)*(k+1)+i);
	V[15] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j-1)+i);
	V[16] = 0.0;
	V[17] = *(slabPressure[5]+(NI-1)*(j-1)+i);
	
	volume = *(volumeDomain+indexP);
	
	MatrixMultCorner(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
	
	resid = *((zone->grid)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS)-(volume/(rho*betatau))*Vg;
	
	diagonal = (delta_t/rho)*(*(diffdia+indexP)) + volume/(rho*betatau);

	*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/diagonal;
	
	residZone=residZone+resid*resid;
//-----------------------------------------------------------------------------------------------------(6)----------
	j = NJ-2; i = NI-2;
	for(k=1;k<NK-2;k++)
	{
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		V[0] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
		V[1] = *(slabPressure[0]+(NJ-1)*k+j);
		V[2] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
		V[3] = *(slabPressure[2]+(NI-1)*k+i);
		V[4] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
		V[5] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
		V[6] = *(slabPressure[2]+(NI-1)*k+i-1);
		V[7] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i-1);
		V[8] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i-1);
		V[9] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i-1);
		V[10] = 0.0;
		V[11] = *(slabPressure[0]+(NJ-1)*k+j-1);
		V[12] = *(slabPressure[0]+(NJ-1)*(k+1)+j);
		V[13] = *(slabPressure[0]+(NJ-1)*(k-1)+j);
		V[14] = *(slabPressure[2]+(NI-1)*(k+1)+i);
		V[15] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j-1)+i);
		V[16] = *(slabPressure[2]+(NI-1)*(k-1)+i);
		V[17] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+i);
		
		volume = *(volumeDomain+indexP);
		
		MatrixMultCorner(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
		
		resid = *((zone->grid)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS)-(volume/(rho*betatau))*Vg;
		
		diagonal = (delta_t/rho)*(*(diffdia+indexP)) + volume/(rho*betatau);
	
		*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/diagonal;
		
		residZone=residZone+resid*resid;
	}
//------------------------------------------------------------------------------------------------------(7)-----------
	k = NK-2; j = 0;
	for(i=1;i<NI-2;i++)
	{
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		V[0] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
		V[1] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
		V[2] = *(slabPressure[3]+(NI-1)*k+i);
		V[3] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
		V[4] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
		V[5] = *(slabPressure[4]+(NI-1)*j+i);
		V[6] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i-1);
		V[7] = *(slabPressure[3]+(NI-1)*k+i-1);
		V[8] = *(slabPressure[4]+(NI-1)*j+i-1);
		V[9] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i-1);
		V[10] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1);
		V[11] = *(slabPressure[3]+(NI-1)*k+i+1);
		V[12] = *(slabPressure[4]+(NI-1)*j+i+1);
		V[13] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i+1);
		V[14] = *(slabPressure[4]+(NI-1)*(j+1)+i);
		V[15] = 0.0;
		V[16] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j+1)+i);
		V[17] = *(slabPressure[3]+(NI-1)*(k-1)+i);
		
		volume = *(volumeDomain+indexP);
		
		MatrixMultCorner(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
		
		resid = *((zone->grid)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS)-(volume/(rho*betatau))*Vg;
		
		diagonal = (delta_t/rho)*(*(diffdia+indexP)) + volume/(rho*betatau);
	
		*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/diagonal;
		
		residZone=residZone+resid*resid;
	}
//-----------------------------------------------------------------------------------------------------(8)--------------
	k = 0; j = 0; i = NI-2;
	indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
	index[0]=k,index[1]=j,index[2]=i;
	
	rho = *(zone->Flowsol->Rho+indexP*flagrho);
	Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
	V[0] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
	V[1] = *(slabPressure[0]+(NJ-1)*k+j);
	V[2] = *(slabPressure[3]+(NI-1)*k+i);
	V[3] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
	V[4] = *(slabPressure[5]+(NI-1)*j+i);
	V[5] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
	V[6] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i-1);
	V[7] = *(slabPressure[3]+(NI-1)*k+i-1);
	V[8] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i-1);
	V[9] = *(slabPressure[5]+(NI-1)*j+i-1);
	V[10] = *(slabPressure[0]+(NJ-1)*k+j+1);
	V[11] = 0.0;
	V[12] = *(slabPressure[0]+(NJ-1)*(k+1)+j);
	V[13] = 0.0;
	V[14] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i);
	V[15] = *(slabPressure[3]+(NI-1)*(k+1)+i);
	V[16] = *(slabPressure[5]+(NI-1)*(j+1)+i);
	V[17] = 0.0;
	
	volume = *(volumeDomain+indexP);
	
	MatrixMultCorner(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
	
	resid = *((zone->grid)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS)-(volume/(rho*betatau))*Vg;
	
	diagonal = (delta_t/rho)*(*(diffdia+indexP)) + volume/(rho*betatau);

	*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/diagonal;
	
	residZone=residZone+resid*resid;
//------------------------------------------------------------------------------------------------------(9)---------
	j = 0; i = 0;
	for(k=1;k<NK-2;k++)
	{
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		V[0] = *(slabPressure[1]+(NJ-1)*k+j);
		V[1] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
		V[2] = *(slabPressure[3]+(NI-1)*k+i);
		V[3] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
		V[4] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
		V[5] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
		V[6] = *(slabPressure[1]+(NJ-1)*k+j+1);
		V[7] = 0.0;
		V[8] = *(slabPressure[1]+(NJ-1)*(k+1)+j);
		V[9] = *(slabPressure[1]+(NJ-1)*(k-1)+j);
		V[10] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1);
		V[11] = *(slabPressure[3]+(NI-1)*k+i+1);
		V[12] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1);
		V[13] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i+1);
		V[14] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i);
		V[15] = *(slabPressure[3]+(NI-1)*(k+1)+i);
		V[16] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j+1)+i);
		V[17] = *(slabPressure[3]+(NI-1)*(k-1)+i);
		
		volume = *(volumeDomain+indexP);
		
		MatrixMultCorner(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
		
		resid = *((zone->grid)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS)-(volume/(rho*betatau))*Vg;
		
		diagonal = (delta_t/rho)*(*(diffdia+indexP)) + volume/(rho*betatau);
	
		*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/diagonal;
		
		residZone=residZone+resid*resid;
	}
//--------------------------------------------------------------------------------------------------(10)------------
	j = 0; i = NI-2;
	for(k=1;k<NK-2;k++)
	{
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		V[0] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
		V[1] = *(slabPressure[0]+(NJ-1)*k+j);
		V[2] = *(slabPressure[3]+(NI-1)*k+i);
		V[3] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
		V[4] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
		V[5] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
		V[6] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i-1);
		V[7] = *(slabPressure[3]+(NI-1)*k+i-1);
		V[8] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i-1);
		V[9] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i-1);
		V[10] = *(slabPressure[0]+(NJ-1)*k+j+1);
		V[11] = 0.0;
		V[12] = *(slabPressure[0]+(NJ-1)*(k+1)+j);
		V[13] = *(slabPressure[0]+(NJ-1)*(k-1)+j);
		V[14] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i);
		V[15] = *(slabPressure[3]+(NI-1)*(k+1)+i);
		V[16] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j+1)+i);
		V[17] = *(slabPressure[3]+(NI-1)*(k-1)+i);
		
		volume = *(volumeDomain+indexP);
		
		MatrixMultCorner(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
		
		resid = *((zone->grid)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS)-(volume/(rho*betatau))*Vg;
		
		diagonal = (delta_t/rho)*(*(diffdia+indexP)) + volume/(rho*betatau);
	
		*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/diagonal;
		
		residZone=residZone+resid*resid;
	}
//---------------------------------------------------------------------------------------------------(11)-----------
	j = NJ-2; i = 0;
	for(k=1;k<NK-2;k++)
	{
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		V[0] = *(slabPressure[1]+(NJ-1)*k+j);
		V[1] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
		V[2] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
		V[3] = *(slabPressure[2]+(NI-1)*k+i);
		V[4] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
		V[5] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
		V[6] = 0.0;
		V[7] = *(slabPressure[1]+(NJ-1)*k+j-1);
		V[8] = *(slabPressure[1]+(NJ-1)*(k+1)+j);
		V[9] = *(slabPressure[1]+(NJ-1)*(k-1)+j);
		V[10] = *(slabPressure[2]+(NI-1)*k+i+1);
		V[11] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i+1);
		V[12] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1);
		V[13] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i+1);
		V[14] = *(slabPressure[2]+(NI-1)*(k+1)+i);
		V[15] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j-1)+i);
		V[16] = *(slabPressure[2]+(NI-1)*(k-1)+i);
		V[17] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+i);
		
		volume = *(volumeDomain+indexP);
		
		MatrixMultCorner(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
		
		resid = *((zone->grid)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS)-(volume/(rho*betatau))*Vg;
		
		diagonal = (delta_t/rho)*(*(diffdia+indexP)) + volume/(rho*betatau);
	
		*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/diagonal;
		
		residZone=residZone+resid*resid;
	}
//-----------------------------------------------------------------------------------------------------(12)---------
	j = 0; i = 0; k = NK-2;
	indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
	index[0]=k,index[1]=j,index[2]=i;
	
	rho = *(zone->Flowsol->Rho+indexP*flagrho);
	Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
	V[0] = *(slabPressure[1]+(NJ-1)*k+j);
	V[1] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
	V[2] = *(slabPressure[3]+(NI-1)*k+i);
	V[3] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
	V[4] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
	V[5] = *(slabPressure[4]+(NI-1)*j+i);
	V[6] = *(slabPressure[1]+(NJ-1)*k+j+1);
	V[7] = 0.0;
	V[8] = 0.0;
	V[9] = *(slabPressure[1]+(NJ-1)*(k-1)+j);
	V[10] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1);
	V[11] = *(slabPressure[3]+(NI-1)*k+i+1);
	V[12] = *(slabPressure[4]+(NI-1)*j+i+1);
	V[13] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i+1);
	V[14] = *(slabPressure[4]+(NI-1)*(j+1)+i);
	V[15] = 0.0;
	V[16] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j+1)+i);
	V[17] = *(slabPressure[3]+(NI-1)*(k-1)+i);
	
	volume = *(volumeDomain+indexP);
	
	MatrixMultCorner(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
	
	resid = *((zone->grid)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS)-(volume/(rho*betatau))*Vg;
	
	diagonal = (delta_t/rho)*(*(diffdia+indexP)) + volume/(rho*betatau);

	*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/diagonal;
	
	residZone=residZone+resid*resid;
//-----------------------------------------------------------------------------------------------------(13)---------
	k = 0; i = 0;
	for(j=1;j<NJ-2;j++)
	{
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		V[0] = *(slabPressure[1]+(NJ-1)*k+j);
		V[1] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
		V[2] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
		V[3] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
		V[4] = *(slabPressure[5]+(NI-1)*j+i);
		V[5] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
		V[6] = *(slabPressure[1]+(NJ-1)*k+j+1);
		V[7] = *(slabPressure[1]+(NJ-1)*k+j-1);
		V[8] = *(slabPressure[1]+(NJ-1)*(k+1)+j);
		V[9] = 0.0;
		V[10] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1);
		V[11] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i+1);
		V[12] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1);
		V[13] = *(slabPressure[5]+(NI-1)*j+i+1);
		V[14] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i);
		V[15] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j-1)+i);
		V[16] = *(slabPressure[5]+(NI-1)*(j+1)+i);
		V[17] = *(slabPressure[5]+(NI-1)*(j-1)+i);
		
		volume = *(volumeDomain+indexP);
		
		MatrixMultCorner(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
		
		resid = *((zone->grid)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS)-(volume/(rho*betatau))*Vg;
		
		diagonal = (delta_t/rho)*(*(diffdia+indexP)) + volume/(rho*betatau);
	
		*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/diagonal;
		
		residZone=residZone+resid*resid;
	}
//-----------------------------------------------------------------------------------------------------(14)----------
	k = NK-2; i = 0;
	for(j=1;j<NJ-2;j++)
	{
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		V[0] = *(slabPressure[1]+(NJ-1)*k+j);
		V[1] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
		V[2] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
		V[3] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
		V[4] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
		V[5] = *(slabPressure[4]+(NI-1)*j+i);
		V[6] = *(slabPressure[1]+(NJ-1)*k+j+1);
		V[7] = *(slabPressure[1]+(NJ-1)*k+j-1);
		V[8] = 0.0;
		V[9] = *(slabPressure[1]+(NJ-1)*(k-1)+j);
		V[10] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1);
		V[11] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i+1);
		V[12] = *(slabPressure[4]+(NI-1)*j+i+1);
		V[13] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i+1);
		V[14] = *(slabPressure[4]+(NI-1)*(j+1)+i);
		V[15] = *(slabPressure[4]+(NI-1)*(j-1)+i);
		V[16] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j+1)+i);
		V[17] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+i);
		
		volume = *(volumeDomain+indexP);
		
		MatrixMultCorner(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
		
		resid = *((zone->grid)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS)-(volume/(rho*betatau))*Vg;
		
		diagonal = (delta_t/rho)*(*(diffdia+indexP)) + volume/(rho*betatau);
	
		*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/diagonal;
		
		residZone=residZone+resid*resid;
	}
//---------------------------------------------------------------------------------------------------(15)------------
	k = 0; i = NI-2;
	for(j=1;j<NJ-2;j++)
	{
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		V[0] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
		V[1] = *(slabPressure[0]+(NJ-1)*k+j);
		V[2] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
		V[3] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
		V[4] = *(slabPressure[5]+(NI-1)*j+i);
		V[5] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
		V[6] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i-1);
		V[7] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i-1);
		V[8] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i-1);
		V[9] = *(slabPressure[5]+(NI-1)*j+i-1);
		V[10] = *(slabPressure[0]+(NJ-1)*k+j+1);
		V[11] = *(slabPressure[0]+(NJ-1)*k+j-1);
		V[12] = *(slabPressure[0]+(NJ-1)*(k+1)+j);
		V[13] = 0.0;
		V[14] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i);
		V[15] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j-1)+i);
		V[16] = *(slabPressure[5]+(NI-1)*(j+1)+i);
		V[17] = *(slabPressure[5]+(NI-1)*(j-1)+i);
		
		volume = *(volumeDomain+indexP);
		
		MatrixMultCorner(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
		
		resid = *((zone->grid)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS)-(volume/(rho*betatau))*Vg;
		
		diagonal = (delta_t/rho)*(*(diffdia+indexP)) + volume/(rho*betatau);
	
		*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/diagonal;
		
		residZone=residZone+resid*resid;
	}

//-----------------------------------------------------------------------------------------------------(16)--------
	k = 0; i = 0; j = NJ-2;
	indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
	index[0]=k,index[1]=j,index[2]=i;
	
	rho = *(zone->Flowsol->Rho+indexP*flagrho);
	Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
	V[0] = *(slabPressure[1]+(NJ-1)*k+j);
	V[1] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
	V[2] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
	V[3] = *(slabPressure[2]+(NI-1)*k+i);
	V[4] = *(slabPressure[5]+(NI-1)*j+i);
	V[5] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
	V[6] = 0.0;
	V[7] = *(slabPressure[1]+(NJ-1)*k+j-1);
	V[8] = *(slabPressure[1]+(NJ-1)*(k+1)+j);
	V[9] = 0.0;
	V[10] = *(slabPressure[2]+(NI-1)*k+i+1);
	V[11] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i+1);
	V[12] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1);
	V[13] = *(slabPressure[5]+(NI-1)*j+i+1);
	V[14] = *(slabPressure[2]+(NI-1)*(k+1)+i);
	V[15] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j-1)+i);
	V[16] = 0.0;
	V[17] = *(slabPressure[5]+(NI-1)*(j-1)+i);
	
	volume = *(volumeDomain+indexP);
	
	MatrixMultCorner(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
	
	resid = *((zone->grid)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS)-(volume/(rho*betatau))*Vg;
	
	diagonal = (delta_t/rho)*(*(diffdia+indexP)) + volume/(rho*betatau);

	*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/diagonal;
	
	residZone=residZone+resid*resid;
//------------------------------------------------------------------------------------------------------(17)---------
	j = 0;
	for(k=1;k<NK-2;k++)
	{
		for(i=1;i<NI-2;i++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			V[0] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			V[1] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			V[2] = *(slabPressure[3]+(NI-1)*k+i);
			V[3] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			V[4] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
			V[5] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
			V[6] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i-1);
			V[7] = *(slabPressure[3]+(NI-1)*k+i-1);
			V[8] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i-1);
			V[9] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i-1);
			V[10] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1);
			V[11] = *(slabPressure[3]+(NI-1)*k+i+1);
			V[12] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1);
			V[13] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i+1);
			V[14] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i);
			V[15] = *(slabPressure[3]+(NI-1)*(k+1)+i);
			V[16] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j+1)+i);
			V[17] = *(slabPressure[3]+(NI-1)*(k-1)+i);
			
			volume = *(volumeDomain+indexP);
			
			MatrixMultCorner(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
			
			resid = *((zone->grid)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS)-(volume/(rho*betatau))*Vg;
			
			diagonal = (delta_t/rho)*(*(diffdia+indexP)) + volume/(rho*betatau);
		
			*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/diagonal;
			
			residZone=residZone+resid*resid;
		}
	}
//-------------------------------------------------------------------------------------------------(18)-----------
	k = 0;
	for(j=1;j<NJ-2;j++)
	{
		for(i=1;i<NI-2;i++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			V[0] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			V[1] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			V[2] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			V[3] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			V[4] = *(slabPressure[5]+(NI-1)*j+i);
			V[5] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
			V[6] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i-1);
			V[7] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i-1);
			V[8] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i-1);
			V[9] = *(slabPressure[5]+(NI-1)*j+i-1);
			V[10] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1);
			V[11] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i+1);
			V[12] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1);
			V[13] = *(slabPressure[5]+(NI-1)*j+i+1);
			V[14] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i);
			V[15] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j-1)+i);
			V[16] = *(slabPressure[5]+(NI-1)*(j+1)+i);
			V[17] = *(slabPressure[5]+(NI-1)*(j-1)+i);
			
			volume = *(volumeDomain+indexP);
			
			MatrixMultCorner(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
			
			resid = *((zone->grid)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS)-(volume/(rho*betatau))*Vg;
			
			diagonal = (delta_t/rho)*(*(diffdia+indexP)) + volume/(rho*betatau);
		
			*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/diagonal;
			
			residZone=residZone+resid*resid;
		}
	}
//--------------------------------------------------------------------------------------------------------(19)--------
	i=0;
	for(k=1;k<NK-2;k++)
	{
		for(j=1;j<NJ-2;j++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			V[0] = *(slabPressure[1]+(NJ-1)*k+j);
			V[1] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			V[2] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			V[3] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			V[4] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
			V[5] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
			V[6] = *(slabPressure[1]+(NJ-1)*k+j+1);
			V[7] = *(slabPressure[1]+(NJ-1)*k+j-1);
			V[8] = *(slabPressure[1]+(NJ-1)*(k+1)+j);
			V[9] = *(slabPressure[1]+(NJ-1)*(k-1)+j);
			V[10] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1);
			V[11] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i+1);
			V[12] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1);
			V[13] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i+1);
			V[14] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i);
			V[15] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j-1)+i);
			V[16] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j+1)+i);
			V[17] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+i);
			
			volume = *(volumeDomain+indexP);
			
			MatrixMultCorner(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
			
			resid = *((zone->grid)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS)-(volume/(rho*betatau))*Vg;
			
			diagonal = (delta_t/rho)*(*(diffdia+indexP)) + volume/(rho*betatau);
		
			*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/diagonal;
			
			residZone=residZone+resid*resid;
		}
	}
//-----------------------------------------------------------------------------------------------------(20)---------
	k = NK-2; j = NJ-2; i = 0;
	indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
	index[0]=k,index[1]=j,index[2]=i;
	
	rho = *(zone->Flowsol->Rho+indexP*flagrho);
	Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
	V[0] = *(slabPressure[1]+(NJ-1)*k+j);
	V[1] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
	V[2] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
	V[3] = *(slabPressure[2]+(NI-1)*k+i);
	V[4] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
	V[5] = *(slabPressure[4]+(NI-1)*j+i);
	V[6] = 0.0;
	V[7] = *(slabPressure[1]+(NJ-1)*k+j-1);
	V[8] = 0.0;
	V[9] = *(slabPressure[1]+(NJ-1)*(k-1)+j);
	V[10] = *(slabPressure[2]+(NI-1)*k+i+1);
	V[11] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i+1);
	V[12] = *(slabPressure[4]+(NI-1)*j+i+1);
	V[13] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i+1);
	V[14] = 0.0;
	V[15] = *(slabPressure[4]+(NI-1)*(j-1)+i);
	V[16] = *(slabPressure[2]+(NI-1)*(k-1)+i);
	V[17] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+i);
	
	volume = *(volumeDomain+indexP);
	
	MatrixMultCorner(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
	
	resid = *((zone->grid)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS)-(volume/(rho*betatau))*Vg;
	
	diagonal = (delta_t/rho)*(*(diffdia+indexP)) + volume/(rho*betatau);

	*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/diagonal;
	
	residZone=residZone+resid*resid;
//-----------------------------------------------------------------------------------------------------(21)---------
	k = NK-2; j = NJ-2;
	for(i=1;i<NI-2;i++)
	{
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		V[0] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
		V[1] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
		V[2] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
		V[3] = *(slabPressure[2]+(NI-1)*k+i);
		V[4] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
		V[5] = *(slabPressure[4]+(NI-1)*j+i);
		V[6] = *(slabPressure[2]+(NI-1)*k+i-1);
		V[7] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i-1);
		V[8] = *(slabPressure[4]+(NI-1)*j+i-1);
		V[9] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i-1);
		V[10] = *(slabPressure[2]+(NI-1)*k+i+1);
		V[11] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i+1);
		V[12] = *(slabPressure[4]+(NI-1)*j+i+1);
		V[13] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i+1);
		V[14] = 0.0;
		V[15] = *(slabPressure[4]+(NI-1)*(j-1)+i);
		V[16] = *(slabPressure[2]+(NI-1)*(k-1)+i);
		V[17] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+i);
		
		volume = *(volumeDomain+indexP);
		
		MatrixMultCorner(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
		
		resid = *((zone->grid)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS)-(volume/(rho*betatau))*Vg;
		
		diagonal = (delta_t/rho)*(*(diffdia+indexP)) + volume/(rho*betatau);
	
		*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/diagonal;
		
		residZone=residZone+resid*resid;
	}
//----------------------------------------------------------------------------------------------------(22)----------
	k = NK-2;
	for(j=1;j<NJ-2;j++)
	{
		for(i=1;i<NI-2;i++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			V[0] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			V[1] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
			V[2] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			V[3] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			V[4] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
			V[5] = *(slabPressure[4]+(NI-1)*j+i);
			V[6] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i-1);
			V[7] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i-1);
			V[8] = *(slabPressure[4]+(NI-1)*j+i-1);
			V[9] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i-1);
			V[10] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1);
			V[11] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i+1);
			V[12] = *(slabPressure[4]+(NI-1)*j+i+1);
			V[13] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i+1);
			V[14] = *(slabPressure[4]+(NI-1)*(j+1)+i);
			V[15] = *(slabPressure[4]+(NI-1)*(j-1)+i);
			V[16] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j+1)+i);
			V[17] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+i);
			
			volume = *(volumeDomain+indexP);
			
			MatrixMultCorner(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
			
			resid = *((zone->grid)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS)-(volume/(rho*betatau))*Vg;
			
			diagonal = (delta_t/rho)*(*(diffdia+indexP)) + volume/(rho*betatau);
		
			*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/diagonal;
			
			residZone=residZone+resid*resid;
		}
	}
//----------------------------------------------------------------------------------------------------(23)----------
	k = NK-2; j = 0; i = NI-2;
	indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
	index[0]=k,index[1]=j,index[2]=i;
	
	rho = *(zone->Flowsol->Rho+indexP*flagrho);
	Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
	V[0] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
	V[1] = *(slabPressure[0]+(NJ-1)*k+j);
	V[2] = *(slabPressure[3]+(NI-1)*k+i);
	V[3] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
	V[4] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
	V[5] = *(slabPressure[4]+(NI-1)*j+i);
	V[6] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i-1);
	V[7] = *(slabPressure[3]+(NI-1)*k+i-1);
	V[8] = *(slabPressure[4]+(NI-1)*j+i-1);
	V[9] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i-1);
	V[10] = *(slabPressure[0]+(NJ-1)*k+j+1);
	V[11] = 0.0;
	V[12] = 0.0;
	V[13] = *(slabPressure[0]+(NJ-1)*(k-1)+j);
	V[14] = *(slabPressure[4]+(NI-1)*(j+1)+i);
	V[15] = 0.0;
	V[16] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j+1)+i);
	V[17] = *(slabPressure[3]+(NI-1)*(k-1)+i);
	
	volume = *(volumeDomain+indexP);
	
	MatrixMultCorner(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
	
	resid = *((zone->grid)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS)-(volume/(rho*betatau))*Vg;
	
	diagonal = (delta_t/rho)*(*(diffdia+indexP)) + volume/(rho*betatau);

	*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/diagonal;
	
	residZone=residZone+resid*resid;
//-----------------------------------------------------------------------------------------------------(24)---------
	k = NK-2; i = NI-2;
	for(j=1;j<NJ-2;j++)
	{
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		V[0] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
		V[1] = *(slabPressure[0]+(NJ-1)*k+j);
		V[2] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
		V[3] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
		V[4] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
		V[5] = *(slabPressure[4]+(NI-1)*j+i);
		V[6] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i-1);
		V[7] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i-1);
		V[8] = *(slabPressure[4]+(NI-1)*j+i-1);
		V[9] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i-1);
		V[10] = *(slabPressure[0]+(NJ-1)*k+j+1);
		V[11] = *(slabPressure[0]+(NJ-1)*k+j-1);
		V[12] = 0.0;
		V[13] = *(slabPressure[0]+(NJ-1)*(k-1)+j);
		V[14] = *(slabPressure[4]+(NI-1)*(j+1)+i);
		V[15] = *(slabPressure[4]+(NI-1)*(j-1)+i);
		V[16] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j+1)+i);
		V[17] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+i);
		
		volume = *(volumeDomain+indexP);
		
		MatrixMultCorner(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
		
		resid = *((zone->grid)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS)-(volume/(rho*betatau))*Vg;
		
		diagonal = (delta_t/rho)*(*(diffdia+indexP)) + volume/(rho*betatau);
	
		*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/diagonal;
		
		residZone=residZone+resid*resid;
	}
//-----------------------------------------------------------------------------------------------------(25)---------
	i = NI-2;
	for(k=1;k<NK-2;k++)
	{
		for(j=1;j<NJ-2;j++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
			V[0] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
			V[1] = *(slabPressure[0]+(NJ-1)*k+j);
			V[2] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
			V[3] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
			V[4] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
			V[5] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
			V[6] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i-1);
			V[7] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i-1);
			V[8] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i-1);
			V[9] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i-1);
			V[10] = *(slabPressure[0]+(NJ-1)*k+j+1);
			V[11] = *(slabPressure[0]+(NJ-1)*k+j-1);
			V[12] = *(slabPressure[0]+(NJ-1)*(k+1)+j);
			V[13] = *(slabPressure[0]+(NJ-1)*(k-1)+j);
			V[14] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i);
			V[15] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j-1)+i);
			V[16] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j+1)+i);
			V[17] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+i);
			
			volume = *(volumeDomain+indexP);
			
			MatrixMultCorner(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
			
			resid = *((zone->grid)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS)-(volume/(rho*betatau))*Vg;
			
			diagonal = (delta_t/rho)*(*(diffdia+indexP)) + volume/(rho*betatau);
		
			*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/diagonal;
			
			residZone=residZone+resid*resid;
		}
	}
//---------------------------------------------------------------------------------------------------(26)-----------
	k = NK-2; j = NJ-2; i = NI-2;
	indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
	index[0]=k,index[1]=j,index[2]=i;
	
	rho = *(zone->Flowsol->Rho+indexP*flagrho);
	Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
	V[0] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
	V[1] = *(slabPressure[0]+(NJ-1)*k+j);
	V[2] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
	V[3] = *(slabPressure[2]+(NI-1)*k+i);
	V[4] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
	V[5] = *(slabPressure[4]+(NI-1)*j+i);
	V[6] = *(slabPressure[2]+(NI-1)*k+i-1);
	V[7] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i-1);
	V[8] = *(slabPressure[4]+(NI-1)*j+i-1);
	V[9] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i-1);
	V[10] = 0.0;
	V[11] = *(slabPressure[0]+(NJ-1)*k+j-1);
	V[12] = 0.0;
	V[13] = *(slabPressure[0]+(NJ-1)*(k-1)+j);
	V[14] = 0.0;
	V[15] = *(slabPressure[4]+(NI-1)*(j-1)+i);
	V[16] = *(slabPressure[2]+(NI-1)*(k-1)+i);
	V[17] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+i);
	
	volume = *(volumeDomain+indexP);
	
	MatrixMultCorner(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
	
	resid = *((zone->grid)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS)-(volume/(rho*betatau))*Vg;
	
	diagonal = (delta_t/rho)*(*(diffdia+indexP)) + volume/(rho*betatau);

	*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/diagonal;
	
	residZone=residZone+resid*resid;
//------------------------------------------- End of Boundary cells ----------------------------------------
	
//------------------------------------------- Interior cells -----------------------------------------(27)------
	for(k=1;k<NK-2;k++)
	for(j=1;j<NJ-2;j++)
	for(i=1;i<NI-2;i++)
	{
		indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
	
		index[0]=k,index[1]=j,index[2]=i;
		
		rho = *(zone->Flowsol->Rho+indexP*flagrho);
		Vg = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i);
		V[0] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i-1);
		V[1] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*j+i+1);
		V[2] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i);
		V[3] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i);
		V[4] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i);
		V[5] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i);
		V[6] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i-1);
		V[7] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i-1);
		V[8] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i-1);
		V[9] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i-1);
		V[10] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1);
		V[11] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i+1);
		V[12] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*j+i+1);
		V[13] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*j+i+1);
		V[14] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i);
		V[15] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j-1)+i);
		V[16] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j+1)+i);
		V[17] = *(pressureDomain+(NI-1)*(NJ-1)*(k-1)+(NI-1)*(j-1)+i);
		
		volume = *(volumeDomain+indexP);
		
		MatrixMultCorner(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
		
		resid = *((zone->grid)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS)-(volume/(rho*betatau))*Vg;
		
		diagonal = (delta_t/rho)*(*(diffdia+indexP)) + volume/(rho*betatau);
	
		*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/diagonal;
		
		residZone=residZone+resid*resid;
	}
//------------------------------------------- End of Interior cells ---------------------------------------------------
	residZone=sqrt(residZone/((NI-1)*(NJ-1)*(NK-1)));
	
	free(LHS);
	free(RHS);
	return(residZone);
}
/**
This function CountCells1() counts the no of cells in a given zone
*/
int CountCells1(Zone *zone)
{
	int totalCells;
	int *size;
	size = zone->Size;
	totalCells=(*(size+0)-1)*(*(size+1)-1)*(*(size+2)-1);
	
	return (totalCells);
}
