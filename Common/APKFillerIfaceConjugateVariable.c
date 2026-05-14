#include "APKFillerIfaceConjugateVariable.h"

int FillerIfaceConjugateVariable(Base *base,int variableTag,Zone *zone,Zone *zoneNeighbor,Interface *facePtr,int iface)
{
	//printf("i am in FillerIfaceConjugateVariable \n");
	int i,j,k,l,m,n,indexP,indexN,indices1[3],indices2[3],indP[3],indN[3];
	int ni,nj,nk,NI,NJ,NK,*ptrsizeP,*ptrsizeN,ifaceN,slabCellIndexP;
	double *betaP,*betaN,digonal,diffCoefP,diffCoefN;//X;
	double *flowsolP,*flowsolN,sum=0.0,res=0.0,convcri=0.0,relaxation=1.0;
	double FdiffP[6],FdiffN[6],*volP,*volN;
	double *slabVelocityP[6],*slabVolumeP[6],*slabVelocityN[6],*slabVolumeN[6];
	
	flowsolP=MakeDomainArray(zone,variableTag);
	flowsolN=MakeDomainArray(zoneNeighbor,variableTag);
	ifaceN = FindConjuIfaceNeighbor(variableTag,zone,zoneNeighbor);
	
	betaP = zone->ZGeom->Beta;
	betaN=zoneNeighbor->ZGeom->Beta;
	volP= zone->ZGeom->Volume;
	volN= zoneNeighbor->ZGeom->Volume;
	ptrsizeP = zone->Size;
	ptrsizeN=zoneNeighbor->Size;
	
	ni = *(ptrsizeP+0);
	nj = *(ptrsizeP+1);
	nk = *(ptrsizeP+2);
	
	NI = *(ptrsizeN+0);
	NJ = *(ptrsizeN+1);
	NK = *(ptrsizeN+2);

	//printf("\nvariable tag = %d iface = %d\n",variableTag,iface);
	
	MakeSlabPtrArray(zone,slabVelocityP,variableTag);/*Called for making slabs of u */
	MakeSlabPtrArray(zone,slabVolumeP,5);/*Called for making slabs of Volumes*/
	MakeSlabPtrArray(zoneNeighbor,slabVelocityN,variableTag);/*Called for making slabs of u */
	MakeSlabPtrArray(zoneNeighbor,slabVolumeN,5);/*Called for making slabs of Volumes*/
	
	convcri = find_convcri(variableTag,base);
	do
	{
		sum=0.0;
		slabCellIndexP = -1;
		switch(iface)
		{
			case 0://Imax interface
				i=ni-2;
				for (k=0;k<nk-1;k++)
				{ 
					for (j=0;j<nj-1;j++)
					{
						slabCellIndexP++;
						indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
						diffCoefP = CalcDiffusionCoefCommon(base,zone,variableTag,indexP);
						indices1[0]=ni,indices1[1]=j+1,indices1[2]=k+1;
						TransformCellIndicesCommon(indices1,facePtr,indices2);
						l=indices2[0]-1;
						m=indices2[1]-1;
						n=indices2[2]-1;
						indexN = (NI-1)*(NJ-1)*n+(NI-1)*m+l;
						indP[0]=k,indP[1]=j,indP[2]=i;
						indN[0]=n,indN[1]=m,indN[2]=l;
						diffCoefN = CalcDiffusionCoefCommon(base,zoneNeighbor,variableTag,indexN);
						DiffusionCommon(indP,flowsolP,volP,ptrsizeP,slabVelocityP,slabVolumeP,betaP,FdiffP);
						DiffusionCommon(indN,flowsolN,volN,ptrsizeN,slabVelocityN,slabVolumeN,betaN,FdiffN);
						switch(ifaceN)
						{
							case 0://east/* Imax*/
								digonal=(*(betaP+(indexP)*6*3+3*1+0)*diffCoefP) - (*(betaN+(indexN)*6*3+3*1+0)*diffCoefN);
								if(digonal <= 0.0)
									digonal = -1.0*digonal;
								res=FdiffN[1]*(diffCoefN)+FdiffP[1]*(-diffCoefP);
								*(slabVelocityP[0]+slabCellIndexP) = *(slabVelocityP[0]+slabCellIndexP) + (res/(digonal));
								*(slabVelocityN[0]+(NJ-1)*n+m) = *(slabVelocityN[0]+(NJ-1)*n+m) + (res/(digonal));
								/*digonal = (*(volP+(indexP))*diffCoefN) + (*(volN+(indexN))*diffCoefP);
								X = (*(volP+(indexP))*diffCoefN)*(*(flowsolN+(indexN))) + (*(volN+(indexN))*diffCoefP)*(*(flowsolP+(indexP)));
								*(slabVelocityP[0]+slabCellIndexP) = X/digonal;*/
								//*(slabVelocityN[0]+(NJ-1)*n+m) = X/digonal;
							break;
							case 1://west/*Imin*/
								digonal = (*(betaP+(indexP)*6*3+3*1+0)*diffCoefP) - (*(betaN+(indexN)*6*3+3*0+0)*diffCoefN);
								if(digonal <= 0.0)
									digonal = -1.0*digonal;
								res=(FdiffN[0]*(-diffCoefN)+FdiffP[1]*(-diffCoefP));
								*(slabVelocityP[0]+slabCellIndexP) = *(slabVelocityP[0]+slabCellIndexP) + (res/(digonal));
								*(slabVelocityN[1]+(NJ-1)*n+m) = *(slabVelocityN[1]+(NJ-1)*n+m) + (res/(digonal));
								/*digonal = (*(volP+(indexP))*diffCoefN) + (*(volN+(indexN))*diffCoefP);
								X = (*(volP+(indexP))*diffCoefN)*(*(flowsolN+(indexN))) + (*(volN+(indexN))*diffCoefP)*(*(flowsolP+(indexP)));
								*(slabVelocityP[0]+slabCellIndexP) = X/digonal;*/
								//*(slabVelocityN[1]+(NJ-1)*n+m) = X/digonal;
							break;
							case 2://north/*Jmax*/
								digonal=(*(betaP+(indexP)*6*3+3*1+0)*diffCoefP) - (*(betaN+(indexN)*6*3+3*3+1)*diffCoefN);
								if(digonal <= 0.0)
									digonal = -1.0*digonal;
								res=(FdiffN[3]*(diffCoefN)+FdiffP[1]*(-diffCoefP));
								*(slabVelocityP[0]+slabCellIndexP) = *(slabVelocityP[0]+slabCellIndexP) + (res/(digonal));
								*(slabVelocityN[2]+(NI-1)*n+l) = *(slabVelocityN[2]+(NI-1)*n+l) + (res/(digonal));
								/**digonal = (*(volP+(indexP))*diffCoefN) + (*(volN+(indexN))*diffCoefP);
								X = (*(volP+(indexP))*diffCoefN)*(*(flowsolN+(indexN))) + (*(volN+(indexN))*diffCoefP)*(*(flowsolP+(indexP)));
								*(slabVelocityP[0]+slabCellIndexP) = X/digonal;*/
								//*(slabVelocityN[2]+(NI-1)*n+l) = X/digonal;
							break;
							case 3://south	/*Jmin*/
								digonal=(*(betaP+(indexP)*6*3+3*1+0)*diffCoefP) - (*(betaN+(indexN)*6*3+3*2+1)*diffCoefN);
								if(digonal <= 0.0)
									digonal = -1.0*digonal;
								res=(FdiffN[2]*(-diffCoefN)+FdiffP[1]*(-diffCoefP));
								*(slabVelocityP[0]+slabCellIndexP) = *(slabVelocityP[0]+slabCellIndexP) + (res/(digonal));
								*(slabVelocityN[3]+(NI-1)*n+l) = *(slabVelocityN[3]+(NI-1)*n+l) + (res/(digonal));
								/*digonal = (*(volP+(indexP))*diffCoefN) + (*(volN+(indexN))*diffCoefP);
								X = (*(volP+(indexP))*diffCoefN)*(*(flowsolN+(indexN))) + (*(volN+(indexN))*diffCoefP)*(*(flowsolP+(indexP)));
								*(slabVelocityP[0]+slabCellIndexP) = X/digonal;*/
								//*(slabVelocityN[3]+(NI-1)*n+l) = X/digonal;
							break;
							case 4://top/*Kmax*/
								digonal=(*(betaP+(indexP)*6*3+3*1+0)*diffCoefP) - (*(betaN+(indexN)*6*3+3*5+2)*diffCoefN);
								if(digonal <= 0.0)
									digonal = -1.0*digonal;
								res=(FdiffN[5]*(diffCoefN)+FdiffP[1]*(-diffCoefP));
								*(slabVelocityP[0]+slabCellIndexP) = *(slabVelocityP[0]+slabCellIndexP) + (res/(digonal));
								*(slabVelocityN[4]+(NI-1)*m+l) = *(slabVelocityN[4]+(NI-1)*m+l) + (res/(digonal));
								/*digonal = (*(volP+(indexP))*diffCoefN) + (*(volN+(indexN))*diffCoefP);
								X = (*(volP+(indexP))*diffCoefN)*(*(flowsolN+(indexN))) + (*(volN+(indexN))*diffCoefP)*(*(flowsolP+(indexP)));
								*(slabVelocityP[0]+slabCellIndexP) = X/digonal;*/
								//*(slabVelocityN[4]+(NI-1)*m+l) = X/digonal;
							break;
							case 5://bottom/*Kmin*/
								digonal=(*(betaP+(indexP)*6*3+3*1+0)*diffCoefP) - (*(betaN+(indexN)*6*3+3*4+2)*diffCoefN);
								if(digonal <= 0.0)
									digonal = -1.0*digonal;
								res=(FdiffN[4]*(-diffCoefN)+FdiffP[1]*(-diffCoefP));
								*(slabVelocityP[0]+slabCellIndexP) = *(slabVelocityP[0]+slabCellIndexP) + (res/(digonal));
								*(slabVelocityN[5]+(NI-1)*m+l) = *(slabVelocityN[5]+(NI-1)*m+l) + (res/(digonal));
								/*digonal = (*(volP+(indexP))*diffCoefN) + (*(volN+(indexN))*diffCoefP);
								X = (*(volP+(indexP))*diffCoefN)*(*(flowsolN+(indexN))) + (*(volN+(indexN))*diffCoefP)*(*(flowsolP+(indexP)));
								*(slabVelocityP[0]+slabCellIndexP) = X/digonal;*/
								//*(slabVelocityN[5]+(NI-1)*m+l) = X/digonal;
							break;
						}
						sum=sum+res*res;
					}
				}
				sum=sqrt(sum/((nk-1)*(nj-1)));//These two lines should be written outside the k-loop <<< Narsimha.
			break;
	
			case 1://Imin interface
				i=0;
				for (k=0;k<nk-1;k++) 
				{
					for (j=0;j<nj-1;j++)
					{
						slabCellIndexP++;
						indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
						diffCoefP = CalcDiffusionCoefCommon(base,zone,variableTag,indexP);
						indices1[0]=0,indices1[1]=j+1,indices1[2]=k+1;
						TransformCellIndicesCommon(indices1,facePtr,indices2);
						l=indices2[0]-1;
						m=indices2[1]-1;
						n=indices2[2]-1;
						indexN = (NI-1)*(NJ-1)*n+(NI-1)*m+l;
						indP[0]=k,indP[1]=j,indP[2]=i;
						indN[0]=n,indN[1]=m,indN[2]=l;
						diffCoefN = CalcDiffusionCoefCommon(base,zoneNeighbor,variableTag,indexN);
						DiffusionCommon(indP,flowsolP,volP,ptrsizeP,slabVelocityP,slabVolumeP,betaP,FdiffP);
						DiffusionCommon(indN,flowsolN,volN,ptrsizeN,slabVelocityN,slabVolumeN,betaN,FdiffN);
						switch(ifaceN)
						{
							case 0://east/* Imax*/
								digonal=-(*(betaP+(indexP)*6*3+3*0+0)*diffCoefP) + (*(betaN+(indexN)*6*3+3*1+0)*diffCoefN);
								if(digonal <= 0.0)
									digonal = -1.0*digonal;
								res=(FdiffN[1]*(-diffCoefN)+FdiffP[0]*(-diffCoefP));
								*(slabVelocityP[1]+slabCellIndexP) = *(slabVelocityP[1]+slabCellIndexP) + (res/(digonal));
								*(slabVelocityN[0]+(NJ-1)*n+m) = *(slabVelocityN[0]+(NJ-1)*n+m) + (res/(digonal));
								/*digonal = (*(volP+(indexP))*diffCoefN) + (*(volN+(indexN))*diffCoefP);
								X = (*(volP+(indexP))*diffCoefN)*(*(flowsolN+(indexN))) + (*(volN+(indexN))*diffCoefP)*(*(flowsolP+(indexP)));
								*(slabVelocityP[1]+slabCellIndexP) = X/digonal;*/
								//*(slabVelocityN[0]+(NJ-1)*n+m) = X/digonal;
							break;
							case 1://west/*Imin*/
								digonal=-(*(betaP+(indexP)*6*3+3*0+0)*diffCoefP) + (*(betaN+(indexN)*6*3+3*0+0)*diffCoefN);
								if(digonal <= 0.0)
									digonal = -1.0*digonal;
								res=(FdiffN[0]*(diffCoefN)+FdiffP[0]*(-diffCoefP));
								*(slabVelocityP[1]+slabCellIndexP) = *(slabVelocityP[1]+slabCellIndexP) + (res/(digonal));
								*(slabVelocityN[1]+(NJ-1)*n+m) = *(slabVelocityN[1]+(NJ-1)*n+m) + (res/(digonal));
								/*digonal = (*(volP+(indexP))*diffCoefN) + (*(volN+(indexN))*diffCoefP);
								X = (*(volP+(indexP))*diffCoefN)*(*(flowsolN+(indexN))) + (*(volN+(indexN))*diffCoefP)*(*(flowsolP+(indexP)));
								*(slabVelocityP[1]+slabCellIndexP) = X/digonal;*/
								//*(slabVelocityN[1]+(NJ-1)*n+m) = X/digonal;
							break;
							case 2://north/*Jmax*/
								digonal=-(*(betaP+(indexP)*6*3+3*0+0)*diffCoefP) + (*(betaN+(indexN)*6*3+3*3+1)*diffCoefN);
								if(digonal <= 0.0)
									digonal = -1.0*digonal;
								res=(FdiffN[3]*(-diffCoefN)+FdiffP[0]*(-diffCoefP));
								*(slabVelocityP[1]+slabCellIndexP) = *(slabVelocityP[1]+slabCellIndexP) + (res/(digonal));
								*(slabVelocityN[2]+(NI-1)*n+l) = *(slabVelocityN[2]+(NI-1)*n+l) + (res/(digonal));
								/*digonal = (*(volP+(indexP))*diffCoefN) + (*(volN+(indexN))*diffCoefP);
								X = (*(volP+(indexP))*diffCoefN)*(*(flowsolN+(indexN))) + (*(volN+(indexN))*diffCoefP)*(*(flowsolP+(indexP)));
								*(slabVelocityP[1]+slabCellIndexP) = X/digonal;*/
								//*(slabVelocityN[2]+(NI-1)*n+l) = X/digonal;
							break;
							case 3://south	/*Jmin*/
								digonal=-(*(betaP+(indexP)*6*3+3*0+0)*diffCoefP) + (*(betaN+(indexN)*6*3+3*2+1)*diffCoefN);
								if(digonal <= 0.0)
									digonal = -1.0*digonal;
								res=(FdiffN[2]*(diffCoefN)+FdiffP[0]*(-diffCoefP));
								*(slabVelocityP[1]+slabCellIndexP) = *(slabVelocityP[1]+slabCellIndexP) + (res/(digonal));
								*(slabVelocityN[3]+(NI-1)*n+l) = *(slabVelocityN[3]+(NI-1)*n+l) + (res/(digonal));
								/*digonal = (*(volP+(indexP))*diffCoefN) + (*(volN+(indexN))*diffCoefP);
								X = (*(volP+(indexP))*diffCoefN)*(*(flowsolN+(indexN))) + (*(volN+(indexN))*diffCoefP)*(*(flowsolP+(indexP)));
								*(slabVelocityP[1]+slabCellIndexP) = X/digonal;*/
								//*(slabVelocityN[3]+(NI-1)*n+l) = X/digonal;
							break;
							case 4://top/*Kmax*/
								digonal=-(*(betaP+(indexP)*6*3+3*0+0)*diffCoefP) + (*(betaN+(indexN)*6*3+3*5+2)*diffCoefN);
								if(digonal <= 0.0)
									digonal = -1.0*digonal;
								res=(FdiffN[5]*(-diffCoefN)+FdiffP[0]*(-diffCoefP));
								*(slabVelocityP[1]+slabCellIndexP) = *(slabVelocityP[1]+slabCellIndexP) + (res/(digonal));
								*(slabVelocityN[4]+(NI-1)*m+l) = *(slabVelocityN[4]+(NI-1)*m+l) + (res/(digonal));
								/*digonal = (*(volP+(indexP))*diffCoefN) + (*(volN+(indexN))*diffCoefP);
								X = (*(volP+(indexP))*diffCoefN)*(*(flowsolN+(indexN))) + (*(volN+(indexN))*diffCoefP)*(*(flowsolP+(indexP)));
								*(slabVelocityP[1]+slabCellIndexP) = X/digonal;*/
								//*(slabVelocityN[4]+(NI-1)*m+l) = X/digonal;
							break;
							case 5://bottom/*Kmin*/
								digonal=-(*(betaP+(indexP)*6*3+3*0+0)*diffCoefP) + (*(betaN+(indexN)*6*3+3*4+2)*diffCoefN);
								if(digonal <= 0.0)
									digonal = -1.0*digonal;
								res=(FdiffN[4]*(diffCoefN)+FdiffP[0]*(-diffCoefP));
								*(slabVelocityP[1]+slabCellIndexP) = *(slabVelocityP[1]+slabCellIndexP) + (res/(digonal));
								*(slabVelocityN[5]+(NI-1)*m+l) = *(slabVelocityN[5]+(NI-1)*m+l) + (res/(digonal));
								/*digonal = (*(volP+(indexP))*diffCoefN) + (*(volN+(indexN))*diffCoefP);
								X = (*(volP+(indexP))*diffCoefN)*(*(flowsolN+(indexN))) + (*(volN+(indexN))*diffCoefP)*(*(flowsolP+(indexP)));
								*(slabVelocityP[1]+slabCellIndexP) = X/digonal;*/
								//*(slabVelocityN[5]+(NI-1)*m+l) = X/digonal;
							break;
						}
					sum=sum+res*res;
					}
				}
				sum=sqrt(sum/((nk-1)*(nj-1)));
			break;
			
			case 2://Jmax interface
				j=nj-2;
				for (k=0;k<nk-1;k++) 
				{
					for (i=0;i<ni-1;i++)
					{
						slabCellIndexP++;
						indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
						diffCoefP = CalcDiffusionCoefCommon(base,zone,variableTag,indexP);
						indices1[0]=i+1,indices1[1]=nj,indices1[2]=k+1;
						TransformCellIndicesCommon(indices1,facePtr,indices2);
						l=indices2[0]-1;
						m=indices2[1]-1;
						n=indices2[2]-1;
						indexN = (NI-1)*(NJ-1)*n+(NI-1)*m+l;
						indP[0]=k,indP[1]=j,indP[2]=i;
						indN[0]=n,indN[1]=m,indN[2]=l;
						diffCoefN = CalcDiffusionCoefCommon(base,zoneNeighbor,variableTag,indexN);
						DiffusionCommon(indP,flowsolP,volP,ptrsizeP,slabVelocityP,slabVolumeP,betaP,FdiffP);
						DiffusionCommon(indN,flowsolN,volN,ptrsizeN,slabVelocityN,slabVolumeN,betaN,FdiffN);
						switch(ifaceN)
						{
							case 0://east/* Imax*/
								digonal=(*(betaP+(indexP)*6*3+3*3+1)*diffCoefP) - (*(betaN+(indexN)*6*3+3*1+0)*diffCoefN);
								if(digonal <= 0.0)
									digonal = -1.0*digonal;
								res=(FdiffN[1]*(diffCoefN)+FdiffP[3]*(-diffCoefP));
								*(slabVelocityP[2]+slabCellIndexP) = *(slabVelocityP[2]+slabCellIndexP) + (res/(digonal));
								*(slabVelocityN[0]+(NJ-1)*n+m) = *(slabVelocityN[0]+(NJ-1)*n+m) + (res/(digonal));
								/*digonal = (*(volP+(indexP))*diffCoefN) + (*(volN+(indexN))*diffCoefP);
								X = (*(volP+(indexP))*diffCoefN)*(*(flowsolN+(indexN))) + (*(volN+(indexN))*diffCoefP)*(*(flowsolP+(indexP)));
								*(slabVelocityP[2]+slabCellIndexP) = X/digonal;*/
								//*(slabVelocityN[0]+(NJ-1)*n+m) = X/digonal;
							break;
							case 1://west/*Imin*/
								digonal=(*(betaP+(indexP)*6*3+3*3+1)*diffCoefP) - (*(betaN+(indexN)*6*3+3*0+0)*diffCoefN);
								if(digonal <= 0.0)
									digonal = -1.0*digonal;
								res=(FdiffN[0]*(-diffCoefN)+FdiffP[3]*(-diffCoefP));
								*(slabVelocityP[2]+slabCellIndexP) = *(slabVelocityP[2]+slabCellIndexP) + (res/(digonal));
								*(slabVelocityN[1]+(NJ-1)*n+m) = *(slabVelocityN[1]+(NJ-1)*n+m) + (res/(digonal));
								/*digonal = (*(volP+(indexP))*diffCoefN) + (*(volN+(indexN))*diffCoefP);
								X = (*(volP+(indexP))*diffCoefN)*(*(flowsolN+(indexN))) + (*(volN+(indexN))*diffCoefP)*(*(flowsolP+(indexP)));
								*(slabVelocityP[2]+slabCellIndexP) = X/digonal;*/
								//*(slabVelocityN[1]+(NJ-1)*n+m) = X/digonal;
							break;
							case 2://north/*Jmax*/
								digonal=(*(betaP+(indexP)*6*3+3*3+1)*diffCoefP) - (*(betaN+(indexN)*6*3+3*3+1)*diffCoefN);
								if(digonal <= 0.0)
									digonal = -1.0*digonal;
								res=(FdiffN[3]*(diffCoefN)+FdiffP[3]*(-diffCoefP));
								*(slabVelocityP[2]+slabCellIndexP) = *(slabVelocityP[2]+slabCellIndexP) + (res/(digonal));
								*(slabVelocityN[2]+(NI-1)*n+l) = *(slabVelocityN[2]+(NI-1)*n+l) + (res/(digonal));
								/*digonal = (*(volP+(indexP))*diffCoefN) + (*(volN+(indexN))*diffCoefP);
								X = (*(volP+(indexP))*diffCoefN)*(*(flowsolN+(indexN))) + (*(volN+(indexN))*diffCoefP)*(*(flowsolP+(indexP)));
								*(slabVelocityP[2]+slabCellIndexP) = X/digonal;*/
								//*(slabVelocityN[2]+(NI-1)*n+l) = X/digonal;
							break;
							case 3://south	/*Jmin*/
								digonal=(*(betaP+(indexP)*6*3+3*3+1)*diffCoefP) - (*(betaN+(indexN)*6*3+3*2+1)*diffCoefN);
								if(digonal <= 0.0)
									digonal = -1.0*digonal;
								res = FdiffN[2]*(-diffCoefN)+FdiffP[3]*(-diffCoefP);
								*(slabVelocityP[2]+slabCellIndexP) = *(slabVelocityP[2]+slabCellIndexP) + (res/(digonal));
								*(slabVelocityN[3]+(NI-1)*n+l) = *(slabVelocityN[3]+(NI-1)*n+l) + (res/(digonal));
								/*digonal = (*(volP+(indexP))*diffCoefN) + (*(volN+(indexN))*diffCoefP);
								X = (*(volP+(indexP))*diffCoefN)*(*(flowsolN+(indexN))) + (*(volN+(indexN))*diffCoefP)*(*(flowsolP+(indexP)));
								*(slabVelocityP[2]+slabCellIndexP) = X/digonal;*/
								//*(slabVelocityN[3]+(NI-1)*n+l) = X/digonal;
		
							break;
							case 4://top/*Kmax*/
								digonal=(*(betaP+(indexP)*6*3+3*3+1)*diffCoefP) - (*(betaN+(indexN)*6*3+3*5+2)*diffCoefN);
								if(digonal <= 0.0)
									digonal = -1.0*digonal;
								res=(FdiffN[5]*(diffCoefN)+FdiffP[3]*(-diffCoefP));
								*(slabVelocityP[2]+slabCellIndexP) = *(slabVelocityP[2]+slabCellIndexP) + (res/(digonal));
								*(slabVelocityN[4]+(NI-1)*m+l) = *(slabVelocityN[4]+(NI-1)*m+l) + (res/(digonal));
								/*digonal = (*(volP+(indexP))*diffCoefN) + (*(volN+(indexN))*diffCoefP);
								X = (*(volP+(indexP))*diffCoefN)*(*(flowsolN+(indexN))) + (*(volN+(indexN))*diffCoefP)*(*(flowsolP+(indexP)));
								*(slabVelocityP[2]+slabCellIndexP) = X/digonal;*/
								//*(slabVelocityN[4]+(NI-1)*m+l) = X/digonal;
							break;
							case 5://bottom/*Kmin*/
								digonal=(*(betaP+(indexP)*6*3+3*3+1)*diffCoefP) - (*(betaN+(indexN)*6*3+3*4+2)*diffCoefN);
								if(digonal <= 0.0)
									digonal = -1.0*digonal;
								res=(FdiffN[4]*(-diffCoefN)+FdiffP[3]*(-diffCoefP));
								*(slabVelocityP[2]+slabCellIndexP) = *(slabVelocityP[2]+slabCellIndexP) + (res/(digonal));
								*(slabVelocityN[5]+(NI-1)*m+l) = *(slabVelocityN[5]+(NI-1)*m+l) + (res/(digonal));
								/*digonal = (*(volP+(indexP))*diffCoefN) + (*(volN+(indexN))*diffCoefP);
								X = (*(volP+(indexP))*diffCoefN)*(*(flowsolN+(indexN))) + (*(volN+(indexN))*diffCoefP)*(*(flowsolP+(indexP)));
								*(slabVelocityP[2]+slabCellIndexP) = X/digonal;*/
								//*(slabVelocityN[5]+(NI-1)*m+l) = X/digonal;
							break;
						}
						sum=sum+res*res;
					}
				}
				sum = sqrt(sum/((nk-1)*(ni-1)));
			break;
			
			case 3://Jmin interface
				j=0;
				for (k=0;k<nk-1;k++)
				{
					for (i=0;i<ni-1;i++)
					{
						slabCellIndexP++;
						indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
						diffCoefP = CalcDiffusionCoefCommon(base,zone,variableTag,indexP);
						indices1[0]=i+1,indices1[1]=0,indices1[2]=k+1;
						TransformCellIndicesCommon(indices1,facePtr,indices2);
						l=indices2[0]-1;
						m=indices2[1]-1;
						n=indices2[2]-1;
						indexN = (NI-1)*(NJ-1)*n+(NI-1)*m+l;
						indP[0]=k,indP[1]=j,indP[2]=i;
						indN[0]=n,indN[1]=m,indN[2]=l;
						diffCoefN = CalcDiffusionCoefCommon(base,zoneNeighbor,variableTag,indexN);
						DiffusionCommon(indP,flowsolP,volP,ptrsizeP,slabVelocityP,slabVolumeP,betaP,FdiffP);
						DiffusionCommon(indN,flowsolN,volN,ptrsizeN,slabVelocityN,slabVolumeN,betaN,FdiffN);
						switch(ifaceN)
						{
							case 0://east/* Imax*/
								digonal=-(*(betaP+(indexP)*6*3+3*2+1)*diffCoefP) + (*(betaN+(indexN)*6*3+3*1+0)*diffCoefN);
								if(digonal <= 0.0)
									digonal = -1.0*digonal;
								res=(FdiffP[1]*(-diffCoefP)+FdiffN[2]*(-diffCoefN));
								*(slabVelocityP[3]+slabCellIndexP) = *(slabVelocityP[3]+slabCellIndexP) + (res/(digonal));
								*(slabVelocityN[0]+(NJ-1)*n+m) = *(slabVelocityN[0]+(NJ-1)*n+m) + (res/(digonal));
								/*digonal = (*(volP+(indexP))*diffCoefN) + (*(volN+(indexN))*diffCoefP);
								X = (*(volP+(indexP))*diffCoefN)*(*(flowsolN+(indexN))) + (*(volN+(indexN))*diffCoefP)*(*(flowsolP+(indexP)));
								*(slabVelocityP[3]+slabCellIndexP) = X/digonal;*/
								//*(slabVelocityN[0]+(NJ-1)*n+m) = X/digonal;
							break;
							case 1://west/*Imin*/
								digonal=-(*(betaP+(indexP)*6*3+3*2+1)*diffCoefP) + (*(betaN+(indexN)*6*3+3*0+0)*diffCoefN);
								if(digonal <= 0.0)
									digonal = -1.0*digonal;
								res=(FdiffN[0]*(diffCoefN)+FdiffP[2]*(-diffCoefP));
								*(slabVelocityP[3]+slabCellIndexP) = *(slabVelocityP[3]+slabCellIndexP) + (res/(digonal));
								*(slabVelocityN[1]+(NJ-1)*n+m) = *(slabVelocityN[1]+(NJ-1)*n+m) + (res/(digonal));
								/*digonal = (*(volP+(indexP))*diffCoefN) + (*(volN+(indexN))*diffCoefP);
								X = (*(volP+(indexP))*diffCoefN)*(*(flowsolN+(indexN))) + (*(volN+(indexN))*diffCoefP)*(*(flowsolP+(indexP)));
								*(slabVelocityP[3]+slabCellIndexP) = X/digonal;*/
								//*(slabVelocityN[1]+(NJ-1)*n+m) = X/digonal;
							break;
							case 2://north/*Jmax*/
								digonal=-(*(betaP+(indexP)*6*3+3*2+1)*diffCoefP) + (*(betaN+(indexN)*6*3+3*3+1)*diffCoefN);
								if(digonal <= 0.0)
									digonal = -1.0*digonal;
								res = (FdiffN[3]*(-diffCoefN)+FdiffP[2]*(-diffCoefP));
								*(slabVelocityP[3]+slabCellIndexP) = *(slabVelocityP[3]+slabCellIndexP) + (res/(digonal));
								*(slabVelocityN[2]+(NI-1)*n+l) = *(slabVelocityN[2]+(NI-1)*n+l) + (res/(digonal));
								/*digonal = (*(volP+(indexP))*diffCoefN) + (*(volN+(indexN))*diffCoefP);
								X = (*(volP+(indexP))*diffCoefN)*(*(flowsolN+(indexN))) + (*(volN+(indexN))*diffCoefP)*(*(flowsolP+(indexP)));
								*(slabVelocityP[3]+slabCellIndexP) = X/digonal;*/
								//*(slabVelocityN[2]+(NI-1)*n+l) = X/digonal;
							break;
							case 3://south	/*Jmin*/
								digonal=-(*(betaP+(indexP)*6*3+3*2+1)*diffCoefP) + (*(betaN+(indexN)*6*3+3*2+1)*diffCoefN);
								if(digonal <= 0.0)
									digonal = -1.0*digonal;
								res=(FdiffN[2]*(diffCoefN)+FdiffP[2]*(-diffCoefP));
								*(slabVelocityP[3]+slabCellIndexP) = *(slabVelocityP[3]+slabCellIndexP) + (res/(digonal));
								*(slabVelocityN[3]+(NI-1)*n+l) = *(slabVelocityN[3]+(NI-1)*n+l) + (res/(digonal));
								/*digonal = (*(volP+(indexP))*diffCoefN) + (*(volN+(indexN))*diffCoefP);
								X = (*(volP+(indexP))*diffCoefN)*(*(flowsolN+(indexN))) + (*(volN+(indexN))*diffCoefP)*(*(flowsolP+(indexP)));
								*(slabVelocityP[3]+slabCellIndexP) = X/digonal;*/
								//*(slabVelocityN[3]+(NI-1)*n+l) = X/digonal;
							break;
							case 4://top/*Kmax*/
								digonal=-(*(betaP+(indexP)*6*3+3*2+1)*diffCoefP) + (*(betaN+(indexN)*6*3+3*5+2)*diffCoefN);
								if(digonal <= 0.0)
									digonal = -1.0*digonal;
								res=(FdiffN[5]*(-diffCoefN)+FdiffP[2]*(-diffCoefP));
								*(slabVelocityP[3]+slabCellIndexP) = *(slabVelocityP[3]+slabCellIndexP) + (res/(digonal));
								*(slabVelocityN[4]+(NI-1)*m+l) = *(slabVelocityN[4]+(NI-1)*m+l) + (res/(digonal));
								/*digonal = (*(volP+(indexP))*diffCoefN) + (*(volN+(indexN))*diffCoefP);
								X = (*(volP+(indexP))*diffCoefN)*(*(flowsolN+(indexN))) + (*(volN+(indexN))*diffCoefP)*(*(flowsolP+(indexP)));
								*(slabVelocityP[3]+slabCellIndexP) = X/digonal;*/
								//*(slabVelocityN[4]+(NI-1)*m+l) = X/digonal;
							break;
							case 5://bottom/*Kmin*/
								digonal=-(*(betaP+(indexP)*6*3+3*2+1)*diffCoefP) + (*(betaN+(indexN)*6*3+3*4+2)*diffCoefN);
								if(digonal <= 0.0)
									digonal = -1.0*digonal;
								res=(FdiffN[4]*(diffCoefN)+FdiffP[2]*(-diffCoefP));
								*(slabVelocityP[3]+slabCellIndexP) = *(slabVelocityP[3]+slabCellIndexP) + (res/(digonal));
								*(slabVelocityN[5]+(NI-1)*m+l) = *(slabVelocityN[5]+(NI-1)*m+l) + (res/(digonal));
								/*digonal = (*(volP+(indexP))*diffCoefN) + (*(volN+(indexN))*diffCoefP);
								X = (*(volP+(indexP))*diffCoefN)*(*(flowsolN+(indexN))) + (*(volN+(indexN))*diffCoefP)*(*(flowsolP+(indexP)));
								*(slabVelocityP[3]+slabCellIndexP) = X/digonal;*/
								//*(slabVelocityN[5]+(NI-1)*m+l) = X/digonal;
							break;	
						}
					sum=sum+res*res;
					}
				}
				sum=sqrt(sum/((nk-1)*(ni-1)));
			break;
	
			case 4://Kmax interface
				k=nk-2;
				for (j=0;j<nj-1;j++) 
				{
					for (i=0;i<ni-1;i++)
					{
						slabCellIndexP++;
						indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
						diffCoefP = CalcDiffusionCoefCommon(base,zone,variableTag,indexP);
						indices1[0]=i+1,indices1[1]=j+1,indices1[2]=nk;
						TransformCellIndicesCommon(indices1,facePtr,indices2);
						l=indices2[0]-1;
						m=indices2[1]-1;
						n=indices2[2]-1;
						indexN = (NI-1)*(NJ-1)*n+(NI-1)*m+l;
						indP[0]=k,indP[1]=j,indP[2]=i;
						indN[0]=n,indN[1]=m,indN[2]=l;
						diffCoefN = CalcDiffusionCoefCommon(base,zoneNeighbor,variableTag,indexN);
						DiffusionCommon(indP,flowsolP,volP,ptrsizeP,slabVelocityP,slabVolumeP,betaP,FdiffP);
						DiffusionCommon(indN,flowsolN,volN,ptrsizeN,slabVelocityN,slabVolumeN,betaN,FdiffN);
						switch(ifaceN)
						{
							case 0://east/* Imax*/
								digonal=(*(betaP+(indexP)*6*3+3*5+2)*diffCoefP) - (*(betaN+(indexN)*6*3+3*1+0)*diffCoefN);
								if(digonal <= 0.0)
									digonal = -1.0*digonal;
								res=(FdiffN[1]*(diffCoefN)+FdiffP[5]*(-diffCoefP));
								*(slabVelocityP[4]+slabCellIndexP) = *(slabVelocityP[4]+slabCellIndexP) + (res/(digonal));
								*(slabVelocityN[0]+(NJ-1)*n+m) = *(slabVelocityN[0]+(NJ-1)*n+m) + (res/(digonal));
								/*digonal = (*(volP+(indexP))*diffCoefN) + (*(volN+(indexN))*diffCoefP);
								X = (*(volP+(indexP))*diffCoefN)*(*(flowsolN+(indexN))) + (*(volN+(indexN))*diffCoefP)*(*(flowsolP+(indexP)));
								*(slabVelocityP[4]+slabCellIndexP) = X/digonal;*/
								//*(slabVelocityN[0]+(NJ-1)*n+m) = X/digonal;
							break;
							case 1://west/*Imin*/
								digonal = (*(betaP+(indexP)*6*3+3*5+2)*diffCoefP) - (*(betaN+(indexN)*6*3+3*0+0)*diffCoefN);
								if(digonal <= 0.0)
									digonal = -1.0*digonal;
								res=(FdiffP[0]*(-diffCoefP)+FdiffN[5]*(-diffCoefN));
								*(slabVelocityP[4]+slabCellIndexP) = *(slabVelocityP[4]+slabCellIndexP) + (res/(digonal));
								*(slabVelocityN[1]+(NJ-1)*n+m) = *(slabVelocityN[1]+(NJ-1)*n+m) + (res/(digonal));
								/*digonal = (*(volP+(indexP))*diffCoefN) + (*(volN+(indexN))*diffCoefP);
								X = (*(volP+(indexP))*diffCoefN)*(*(flowsolN+(indexN))) + (*(volN+(indexN))*diffCoefP)*(*(flowsolP+(indexP)));
								*(slabVelocityP[4]+slabCellIndexP) = X/digonal;*/
								//*(slabVelocityN[1]+(NJ-1)*n+m) = X/digonal;
							break;
							case 2://north/*Jmax*/
								digonal=(*(betaP+(indexP)*6*3+3*5+2)*diffCoefP) - (*(betaN+(indexN)*6*3+3*3+1)*diffCoefN);
								if(digonal <= 0.0)
									digonal = -1.0*digonal;
								res=(FdiffN[3]*(diffCoefN)+FdiffP[5]*(-diffCoefP));
								*(slabVelocityP[4]+slabCellIndexP) = *(slabVelocityP[4]+slabCellIndexP) + (res/(digonal));
								*(slabVelocityN[2]+(NI-1)*n+l) = *(slabVelocityN[2]+(NI-1)*n+l) + (res/(digonal));
								/*digonal = (*(volP+(indexP))*diffCoefN) + (*(volN+(indexN))*diffCoefP);
								X = (*(volP+(indexP))*diffCoefN)*(*(flowsolN+(indexN))) + (*(volN+(indexN))*diffCoefP)*(*(flowsolP+(indexP)));
								*(slabVelocityP[4]+slabCellIndexP) = X/digonal;*/
								//*(slabVelocityN[2]+(NI-1)*n+l) = X/digonal;
							break;
							case 3://south	/*Jmin*/
								digonal=(*(betaP+(indexP)*6*3+3*5+2)*diffCoefP) - (*(betaN+(indexN)*6*3+3*2+1)*diffCoefN);
								if(digonal <= 0.0)
									digonal = -1.0*digonal;
								res=(FdiffN[2]*(-diffCoefN)+FdiffP[5]*(-diffCoefP));
								*(slabVelocityP[4]+slabCellIndexP) = *(slabVelocityP[4]+slabCellIndexP) + (res/(digonal));
								*(slabVelocityN[3]+(NI-1)*n+l) = *(slabVelocityN[3]+(NI-1)*n+l) + (res/(digonal));
								/*digonal = (*(volP+(indexP))*diffCoefN) + (*(volN+(indexN))*diffCoefP);
								X = (*(volP+(indexP))*diffCoefN)*(*(flowsolN+(indexN))) + (*(volN+(indexN))*diffCoefP)*(*(flowsolP+(indexP)));
								*(slabVelocityP[4]+slabCellIndexP) = X/digonal;*/
								//*(slabVelocityN[3]+(NI-1)*n+l) = X/digonal;
							break;
							case 4://top/*Kmax*/
								digonal=(*(betaP+(indexP)*6*3+3*5+2)*diffCoefP) - (*(betaN+(indexN)*6*3+3*5+2)*diffCoefN);
								if(digonal <= 0.0)
									digonal = -1.0*digonal;
								res=(FdiffN[5]*(diffCoefN)+FdiffP[5]*(-diffCoefP));
								*(slabVelocityP[4]+slabCellIndexP) = *(slabVelocityP[4]+slabCellIndexP) + (res/(digonal));
								*(slabVelocityN[4]+(NI-1)*m+l) = *(slabVelocityN[4]+(NI-1)*m+l) + (res/(digonal));
								/*digonal = (*(volP+(indexP))*diffCoefN) + (*(volN+(indexN))*diffCoefP);
								X = (*(volP+(indexP))*diffCoefN)*(*(flowsolN+(indexN))) + (*(volN+(indexN))*diffCoefP)*(*(flowsolP+(indexP)));
								*(slabVelocityP[4]+slabCellIndexP) = X/digonal;*/
								//*(slabVelocityN[4]+(NI-1)*m+l) = X/digonal;
							break;
							case 5://bottom/*Kmin*/
								digonal=(*(betaP+(indexP)*6*3+3*5+2)*diffCoefP) - (*(betaN+(indexN)*6*3+3*4+2)*diffCoefN);
								if(digonal <= 0.0)
									digonal = -1.0*digonal;
								res=(FdiffN[4]*(-diffCoefN)+FdiffP[5]*(-diffCoefP));
								*(slabVelocityP[4]+slabCellIndexP) = *(slabVelocityP[4]+slabCellIndexP) + (res/(digonal));
								*(slabVelocityN[5]+(NI-1)*m+l) = *(slabVelocityN[5]+(NI-1)*m+l) + (res/(digonal));
								/*digonal = (*(volP+(indexP))*diffCoefN) + (*(volN+(indexN))*diffCoefP);
								X = (*(volP+(indexP))*diffCoefN)*(*(flowsolN+(indexN))) + (*(volN+(indexN))*diffCoefP)*(*(flowsolP+(indexP)));
								*(slabVelocityP[4]+slabCellIndexP) = X/digonal;*/
								//*(slabVelocityN[5]+(NI-1)*m+l) = X/digonal;
							break;
						}
						sum=sum+res*res;
					}
				}
				sum=sqrt(sum/((ni-1)*(nj-1)));
			break;
	
			case 5://Kmin interface
				k=0;
				for (j=0;j<nj-1;j++) 
				{
					for (i=0;i<ni-1;i++)
					{
						slabCellIndexP++;
						indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
						diffCoefP = CalcDiffusionCoefCommon(base,zone,variableTag,indexP);
						indices1[0]=i+1,indices1[1]=j+1,indices1[2]=0;
						TransformCellIndicesCommon(indices1,facePtr,indices2);
						l=indices2[0]-1;
						m=indices2[1]-1;
						n=indices2[2]-1;
						indexN = (NI-1)*(NJ-1)*n+(NI-1)*m+l;
						indP[0]=k,indP[1]=j,indP[2]=i;
						indN[0]=n,indN[1]=m,indN[2]=l;
						diffCoefN = CalcDiffusionCoefCommon(base,zoneNeighbor,variableTag,indexN);
						DiffusionCommon(indP,flowsolP,volP,ptrsizeP,slabVelocityP,slabVolumeP,betaP,FdiffP);
						DiffusionCommon(indN,flowsolN,volN,ptrsizeN,slabVelocityN,slabVolumeN,betaN,FdiffN);
						switch(ifaceN)
						{
							case 0://east/* Imax*/
								digonal=-(*(betaP+(indexP)*6*3+3*4+2)*diffCoefP) + (*(betaN+(indexN)*6*3+3*1+0)*diffCoefN);
								if(digonal <= 0.0)
									digonal = -1.0*digonal;
								res=(FdiffN[1]*(-diffCoefN)+FdiffP[4]*(-diffCoefP));
								*(slabVelocityP[5]+slabCellIndexP) = *(slabVelocityP[5]+slabCellIndexP) + (res/(digonal));
								*(slabVelocityN[0]+(NJ-1)*n+m) = *(slabVelocityN[0]+(NJ-1)*n+m) + (res/(digonal));
								/*digonal = (*(volP+(indexP))*diffCoefN) + (*(volN+(indexN))*diffCoefP);
								X = (*(volP+(indexP))*diffCoefN)*(*(flowsolN+(indexN))) + (*(volN+(indexN))*diffCoefP)*(*(flowsolP+(indexP)));
								*(slabVelocityP[5]+slabCellIndexP) = X/digonal;*/
								//*(slabVelocityN[0]+(NJ-1)*n+m) = X/digonal;
							break;
							case 1://west/*Imin*/
								digonal=-(*(betaP+(indexP)*6*3+3*4+2)*diffCoefP) + (*(betaN+(indexN)*6*3+3*0+0)*diffCoefN);
								if(digonal <= 0.0)
									digonal = -1.0*digonal;
								res=(FdiffN[0]*(diffCoefN)+FdiffP[4]*(-diffCoefP));
								*(slabVelocityP[5]+slabCellIndexP) = *(slabVelocityP[5]+slabCellIndexP) + (res/(digonal));
								*(slabVelocityN[1]+(NJ-1)*n+m) = *(slabVelocityN[1]+(NJ-1)*n+m) + (res/(digonal));
								/*digonal = (*(volP+(indexP))*diffCoefN) + (*(volN+(indexN))*diffCoefP);
								X = (*(volP+(indexP))*diffCoefN)*(*(flowsolN+(indexN))) + (*(volN+(indexN))*diffCoefP)*(*(flowsolP+(indexP)));
								*(slabVelocityP[5]+slabCellIndexP) = X/digonal;*/
								//*(slabVelocityN[1]+(NJ-1)*n+m) = X/digonal;
							break;
							case 2://north/*Jmax*/
								digonal=-(*(betaP+(indexP)*6*3+3*4+2)*diffCoefP) + (*(betaN+(indexN)*6*3+3*3+1)*diffCoefN);
								if(digonal <= 0.0)
									digonal = -1.0*digonal;
								res=(FdiffN[3]*(-diffCoefN)+FdiffP[4]*(-diffCoefP));
								*(slabVelocityP[5]+slabCellIndexP) = *(slabVelocityP[5]+slabCellIndexP) + (res/(digonal));
								*(slabVelocityN[2]+(NI-1)*n+l) = *(slabVelocityN[2]+(NI-1)*n+l) + (res/(digonal));
								/*digonal = (*(volP+(indexP))*diffCoefN) + (*(volN+(indexN))*diffCoefP);
								X = (*(volP+(indexP))*diffCoefN)*(*(flowsolN+(indexN))) + (*(volN+(indexN))*diffCoefP)*(*(flowsolP+(indexP)));
								*(slabVelocityP[5]+slabCellIndexP) = X/digonal;*/
								//*(slabVelocityN[2]+(NI-1)*n+l) = X/digonal;
							break;
							case 3://south	/*Jmin*/
								digonal=-(*(betaP+(indexP)*6*3+3*4+2)*diffCoefP) + (*(betaN+(indexN)*6*3+3*2+1)*diffCoefN);
								if(digonal <= 0.0)
									digonal = -1.0*digonal;
								res=(FdiffN[2]*(diffCoefN)+FdiffP[4]*(-diffCoefP));
								*(slabVelocityP[5]+slabCellIndexP) = *(slabVelocityP[5]+slabCellIndexP) + (res/(digonal));
								*(slabVelocityN[3]+(NI-1)*n+l) = *(slabVelocityN[3]+(NI-1)*n+l) + (res/(digonal));
								/*digonal = (*(volP+(indexP))*diffCoefN) + (*(volN+(indexN))*diffCoefP);
								X = (*(volP+(indexP))*diffCoefN)*(*(flowsolN+(indexN))) + (*(volN+(indexN))*diffCoefP)*(*(flowsolP+(indexP)));
								*(slabVelocityP[5]+slabCellIndexP) = X/digonal;*/
								//*(slabVelocityN[3]+(NI-1)*n+l) = X/digonal;
							break;
							case 4://top/*Kmax*/
								digonal=-(*(betaP+(indexP)*6*3+3*4+2)*diffCoefP) + (*(betaN+(indexN)*6*3+3*5+2)*diffCoefN);
								if(digonal <= 0.0)
									digonal = -1.0*digonal;
								res=(FdiffN[5]*(-diffCoefN)+FdiffP[4]*(-diffCoefP));
								*(slabVelocityP[5]+slabCellIndexP) = *(slabVelocityP[5]+slabCellIndexP) + (res/(digonal));
								*(slabVelocityN[4]+(NI-1)*m+l) = *(slabVelocityN[4]+(NI-1)*m+l) + (res/(digonal));
								/*digonal = (*(volP+(indexP))*diffCoefN) + (*(volN+(indexN))*diffCoefP);
								X = (*(volP+(indexP))*diffCoefN)*(*(flowsolN+(indexN))) + (*(volN+(indexN))*diffCoefP)*(*(flowsolP+(indexP)));
								*(slabVelocityP[5]+slabCellIndexP) = X/digonal;*/
								//*(slabVelocityN[4]+(NI-1)*m+l) = X/digonal;
							break;
							case 5://bottom/*Kmin*/
								digonal=-(*(betaP+(indexP)*6*3+3*4+2)*diffCoefP) + (*(betaN+(indexN)*6*3+3*4+2)*diffCoefN);
								if(digonal <= 0.0)
									digonal = -1.0*digonal;
								res=(FdiffN[4]*(diffCoefN)+FdiffP[4]*(-diffCoefP));
								*(slabVelocityP[5]+slabCellIndexP) = *(slabVelocityP[5]+slabCellIndexP) + (res/(digonal));
								*(slabVelocityN[5]+(NI-1)*m+l) = *(slabVelocityN[5]+(NI-1)*m+l) + (res/(digonal));
								/*digonal = (*(volP+(indexP))*diffCoefN) + (*(volN+(indexN))*diffCoefP);
								X = (*(volP+(indexP))*diffCoefN)*(*(flowsolN+(indexN))) + (*(volN+(indexN))*diffCoefP)*(*(flowsolP+(indexP)));
								*(slabVelocityP[5]+slabCellIndexP) = X/digonal;*/
								//*(slabVelocityN[5]+(NI-1)*m+l) = X/digonal;
							break;
		
						}
						sum=sum+res*res;
					}
				}
				sum=sqrt(sum/((ni-1)*(nj-1)));
			break;
		}
		
		if(sum > 1.0e20)
		{
			printf("variableTag=%d Zone:%s Interface Direction=%d residual (%e) has diverged beyond 1.0e20 ...\n***********\nSuggestion:\n***********\n1) Improve the quality of GRID (reduce SKEWNESS) \nOR\n2) Use a lower UNDER RELAXATION parameter for the variable[%d] Present relaxation is: %lf (variableTag=4 <==> Pressure)\n", variableTag,zone->zonename,facePtr->Direction,sum,variableTag,relaxation);
			exit(24);
		}
	}while(sum > convcri);

	facePtr->ConjugateRms = sum;
	return (0);
}
