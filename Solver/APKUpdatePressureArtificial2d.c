#include "APKUpdatePressureArtificial2d.h"
/**
The UpdatePressureArtificial2d() does the Gauss Seidel Iterations for modified Pressure Poisson equation for Artificial Compressibility
in 2 dimensional geometries. See Eqn(3.70), Section 3.9.0.1, pp 46, in Project Report No 1, Dec 2005.
*/
/*
double UpdatePressureArtificial2d(Base *base,Zone *zone)
{
	double residual=0.0,residZone=0.0,*pressureDomain,sumMassFlux,sumDiff,diffdia;
	double *prevpressureDomain,volume;
	double *massFlux,pressureDiff[6];
	int NI,NJ,NK,i,j,k,iFace,index[3],indexP,*ptrSize,flagrho;
	double *slabVolume[6],*slabPressure[6];
	double *beta,*volumeDomain,*diffusionCoeffdomain,rho;
	double delta_t,CoeffN[6],betatau=1.0;
	double *Residue,*CorrectionX,*DiagonalC,*NCoeff;
	
	delta_t = base->BaseParameters->TimeStep;
	pressureDomain = zone->Flowsol->P;
	prevpressureDomain = zone->PrevFlowsol->P;
	volumeDomain = zone->ZGeom->Volume;
	ptrSize = zone->Size;
	beta=zone->ZGeom->Beta;
	
	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
	
	Residue=(double *)malloc((NJ-1)*sizeof(double));
	CorrectionX=(double *)malloc((NJ-1)*sizeof(double));
	DiagonalC=(double *)malloc((NJ-1)*sizeof(double));
	NCoeff=(double *)malloc(6*(NJ-1)*sizeof(double));
	
	MakeSlabPtrArray(zone,slabPressure,4);//Called for making slabs of Pressure
	MakeSlabPtrArray(zone,slabVolume,5);//Called for making slabs of Volumes
	diffusionCoeffdomain = MakeDiffusionCoeffDomain(zone,4);
	flagrho = zone->Flowsol->flagRho;
	//printf("flagrho = %d\n",flagrho);
	for(k=0;k<NK-1;k++)
	for(i=0;i<NI-1;i++)
	{
		for(j=0;j<NJ-1;j++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			index[0]=k,index[1]=j,index[2]=i;
			
			volume=*(volumeDomain+indexP);
			
			rho = *(zone->Flowsol->Rho+indexP*flagrho);
			
			diffdia = (-delta_t/rho)*Diffusion2dWithConstDiffCoeff(index,pressureDomain,volumeDomain,ptrSize,slabPressure,slabVolume,diffusionCoeffdomain,beta,pressureDiff,CoeffN);
			//printf("diffdia = %lf\n",diffdia);
			massFlux=zone->Flowsol->Massflux+indexP*6;//assign the pointer here
			
			sumDiff=0.0;
			sumMassFlux=0.0;
			for(iFace=0;iFace<4;iFace++)
			{
				sumDiff=sumDiff+pressureDiff[iFace];
				sumMassFlux=sumMassFlux +(*(massFlux+iFace)); 
				*(NCoeff+j*6+iFace) = (delta_t/rho)*CoeffN[iFace];
			}
			residual= -sumMassFlux+(delta_t/rho)*sumDiff+(volume/(rho*betatau))*(*(prevpressureDomain+indexP) - *(pressureDomain+indexP));
			//printf("rho= %lf  sumDiff= %lf  residual = %lf\n",rho,sumDiff,residual);
			*(Residue+j)=residual;
			
			*(DiagonalC+j)=diffdia + (volume/(rho*betatau));
		}
		TDMA(NJ-1,Residue,DiagonalC,NCoeff,CorrectionX);
		for(j=0;j<NJ-1;j++)
		{
			indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
			*(pressureDomain+indexP) = *(pressureDomain+indexP) + *(CorrectionX+j);
			residZone=residZone+(*(CorrectionX+j))*(*(CorrectionX+j));
		}
	}
	residZone=sqrt(residZone/((NI-1)*(NJ-1)*(NK-1)));
	
	free(Residue);
	free(CorrectionX);
	free(DiagonalC);
	free(NCoeff);
	return(residZone);
}
*/
double UpdatePressureArtificial2d(Base *base,Zone *zone)
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
	volumeDomain = zone->ZGeom->Volume;
	
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
	V[5] = *(slabPressure[4]+(NI-1)*j+i);
	V[6] = *(slabPressure[1]+(NJ-1)*k+j+1);
	V[7] = 0.0;
	V[8] = 0.0;
	V[9] = 0.0;
	V[10] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1);
	V[11] = *(slabPressure[3]+(NI-1)*k+i+1);
	V[12] = *(slabPressure[4]+(NI-1)*j+i+1);
	V[13] = *(slabPressure[5]+(NI-1)*j+i+1);
	V[14] = *(slabPressure[4]+(NI-1)*(j+1)+i);
	V[15] = 0.0;
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
		V[5] = *(slabPressure[4]+(NI-1)*j+i);
		V[6] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i-1);
		V[7] = *(slabPressure[3]+(NI-1)*k+i-1);
		V[8] = *(slabPressure[4]+(NI-1)*j+i-1);
		V[9] = *(slabPressure[5]+(NI-1)*j+i-1);
		V[10] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1);
		V[11] = *(slabPressure[3]+(NI-1)*k+i+1);
		V[12] = *(slabPressure[4]+(NI-1)*j+i+1);
		V[13] = *(slabPressure[5]+(NI-1)*j+i+1);
		V[14] = *(pressureDomain+(NI-1)*(NJ-1)*(k+1)+(NI-1)*(j+1)+i);
		V[15] = 0.0;
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
		V[5] = *(slabPressure[4]+(NI-1)*j+i);
		V[6] = *(slabPressure[2]+(NI-1)*k+i-1);
		V[7] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i-1);
		V[8] = *(slabPressure[4]+(NI-1)*j+i-1);
		V[9] = *(slabPressure[5]+(NI-1)*j+i-1);
		V[10] = *(slabPressure[2]+(NI-1)*k+i+1);
		V[11] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i+1);
		V[12] = *(slabPressure[4]+(NI-1)*j+i+1);
		V[13] = *(slabPressure[5]+(NI-1)*j+i+1);
		V[14] = 0.0;
		V[15] = *(slabPressure[4]+(NI-1)*(j-1)+i);
		V[16] = 0.0;
		V[17] = *(slabPressure[5]+(NI-1)*(j-1)+i);
		
		volume = *(volumeDomain+indexP);
		
		MatrixMultCorner(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
		
		resid = *((zone->grid)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS)-(volume/(rho*betatau))*Vg;
		
		diagonal = (delta_t/rho)*(*(diffdia+indexP)) + volume/(rho*betatau);
	
		*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/diagonal;
		
		residZone=residZone+resid*resid;
	}
//------------------------------------------------------------------------------------------------------(4)----------
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
	V[5] = *(slabPressure[4]+(NI-1)*j+i);
	V[6] = *(slabPressure[2]+(NI-1)*k+i-1);
	V[7] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i-1);
	V[8] = *(slabPressure[4]+(NI-1)*j+i-1);
	V[9] = *(slabPressure[5]+(NI-1)*j+i-1);
	V[10] = 0.0;
	V[11] = *(slabPressure[0]+(NJ-1)*k+j-1);
	V[12] = 0.0;
	V[13] = 0.0;
	V[14] = 0.0;
	V[15] = *(slabPressure[4]+(NI-1)*(j-1)+i);
	V[16] = 0.0;
	V[17] = *(slabPressure[5]+(NI-1)*(j-1)+i);
	
	volume = *(volumeDomain+indexP);
	
	MatrixMultCorner(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
	
	resid = *((zone->grid)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS)-(volume/(rho*betatau))*Vg;
	
	diagonal = (delta_t/rho)*(*(diffdia+indexP)) + volume/(rho*betatau);

	*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/diagonal;
	
	residZone=residZone+resid*resid;
//-----------------------------------------------------------------------------------------------------(5)--------------
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
	V[5] = *(slabPressure[4]+(NI-1)*j+i);
	V[6] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i-1);
	V[7] = *(slabPressure[3]+(NI-1)*k+i-1);
	V[8] = *(slabPressure[4]+(NI-1)*j+i-1);
	V[9] = *(slabPressure[5]+(NI-1)*j+i-1);
	V[10] = *(slabPressure[0]+(NJ-1)*k+j+1);
	V[11] = 0.0;
	V[12] = 0.0;
	V[13] = 0.0;
	V[14] = *(slabPressure[4]+(NI-1)*(j-1)+i);
	V[15] = 0.0;
	V[16] = *(slabPressure[5]+(NI-1)*(j+1)+i);
	V[17] = 0.0;
	
	volume = *(volumeDomain+indexP);
	
	MatrixMultCorner(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
	
	resid = *((zone->grid)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS)-(volume/(rho*betatau))*Vg;
	
	diagonal = (delta_t/rho)*(*(diffdia+indexP)) + volume/(rho*betatau);

	*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/diagonal;
	
	residZone=residZone+resid*resid;
//------------------------------------------------------------------------------------------------------(6)---------
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
		V[5] = *(slabPressure[4]+(NI-1)*j+i);
		V[6] = *(slabPressure[1]+(NJ-1)*k+j+1);
		V[7] = *(slabPressure[1]+(NJ-1)*k+j-1);
		V[8] = 0.0;
		V[9] = 0.0;
		V[10] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1);
		V[11] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i+1);
		V[12] = *(slabPressure[4]+(NI-1)*j+i+1);
		V[13] = *(slabPressure[5]+(NI-1)*j+i+1);
		V[14] = *(slabPressure[4]+(NI-1)*(j+1)+i);
		V[15] = *(slabPressure[4]+(NI-1)*(j-1)+i);
		V[16] = *(slabPressure[5]+(NI-1)*(j+1)+i);
		V[17] = *(slabPressure[5]+(NI-1)*(j-1)+i);
		
		volume = *(volumeDomain+indexP);
		
		MatrixMultCorner(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
		
		resid = *((zone->grid)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS)-(volume/(rho*betatau))*Vg;
		
		diagonal = (delta_t/rho)*(*(diffdia+indexP)) + volume/(rho*betatau);
	
		*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/diagonal;
		
		residZone=residZone+resid*resid;
	}
//-----------------------------------------------------------------------------------------------------(7)----------
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
		V[5] = *(slabPressure[4]+(NI-1)*j+i);
		V[6] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i-1);
		V[7] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i-1);
		V[8] = *(slabPressure[4]+(NI-1)*j+i-1);
		V[9] = *(slabPressure[5]+(NI-1)*j+i-1);
		V[10] = *(slabPressure[0]+(NJ-1)*k+j+1);
		V[11] = *(slabPressure[0]+(NJ-1)*k+j-1);
		V[12] = 0.0;
		V[13] = 0.0;
		V[14] = *(slabPressure[4]+(NI-1)*(j+1)+i);
		V[15] = *(slabPressure[4]+(NI-1)*(j-1)+i);
		V[16] = *(slabPressure[5]+(NI-1)*(j+1)+i);
		V[17] = *(slabPressure[5]+(NI-1)*(j-1)+i);
		
		volume = *(volumeDomain+indexP);
		
		MatrixMultCorner(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
		
		resid = *((zone->grid)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS)-(volume/(rho*betatau))*Vg;
		
		diagonal = (delta_t/rho)*(*(diffdia+indexP)) + volume/(rho*betatau);
	
		*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/diagonal;
		
		residZone=residZone+resid*resid;
	}
//-----------------------------------------------------------------------------------------------------(8)--------
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
	V[5] = *(slabPressure[4]+(NI-1)*j+i);
	V[6] = 0.0;
	V[7] = *(slabPressure[1]+(NJ-1)*k+j-1);
	V[8] = 0.0;
	V[9] = 0.0;
	V[10] = *(slabPressure[2]+(NI-1)*k+i+1);
	V[11] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i+1);
	V[12] = *(slabPressure[4]+(NI-1)*j+i+1);
	V[13] = *(slabPressure[5]+(NI-1)*j+i+1);
	V[14] = 0.0;
	V[15] = *(slabPressure[4]+(NI-1)*(j-1)+i);
	V[16] = 0.0;
	V[17] = *(slabPressure[5]+(NI-1)*(j-1)+i);
	
	volume = *(volumeDomain+indexP);
	
	MatrixMultCorner(index,ptrSize,Vg,V,diffdia,CoeffN,LHS,RHS);
	
	resid = *((zone->grid)->GRHS+indexP)-(delta_t/rho)*(*RHS)-(delta_t/rho)*(*LHS)-(volume/(rho*betatau))*Vg;
	
	diagonal = (delta_t/rho)*(*(diffdia+indexP)) + volume/(rho*betatau);

	*(pressureDomain+indexP) = *(pressureDomain+indexP) + omega*resid/diagonal;
	
	residZone=residZone+resid*resid;
//------------------------------------------------End Boundary Cells-------------------------------------------------------------
//--------------------------------------------------Interior Cells----------------------------------------------------(9)---------
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
			V[5] = *(slabPressure[4]+(NI-1)*j+i);
			V[6] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i-1);
			V[7] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i-1);
			V[8] = *(slabPressure[4]+(NI-1)*j+i-1);
			V[9] = *(slabPressure[5]+(NI-1)*j+i-1);
			V[10] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j+1)+i+1);
			V[11] = *(pressureDomain+(NI-1)*(NJ-1)*k+(NI-1)*(j-1)+i+1);
			V[12] = *(slabPressure[4]+(NI-1)*j+i+1);
			V[13] = *(slabPressure[5]+(NI-1)*j+i+1);
			V[14] = *(slabPressure[4]+(NI-1)*(j+1)+i);
			V[15] = *(slabPressure[4]+(NI-1)*(j-1)+i);
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
//------------------------------------------- End of Interior cells ---------------------------------------------------
	residZone=sqrt(residZone/((NI-1)*(NJ-1)*(NK-1)));
	
	free(LHS);
	free(RHS);
	return(residZone);
}


