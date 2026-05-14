#include "APKNonOrthogonalNeumannIC.h"
/**This function is used to apply Neumann condition at the SOLID-FLUID zone interface having NON-ORTHOGONAL/ORTHOGONAL
GRID. This function is called for the variables 'Pressure', 'k' in standard k-epsilon and RNG k-epsilon models and 'epsilon' in Lam-Bremhorst model.

In general this function can be called for any variable which needs Neumann (non-homogeneous/homogeneous) Condition 
to be applied at the SOLID-FLUID interface having  NON-ORTHOGONAL/ORTHOGONAL GRID.

Note that Inward flux is positive and outward flux is negative */
double NonOrthogonalNeumannIC(Base *base, Zone *z,Interface *facePtr,int vartag, double inwardFlux)
{
	int ni,nj,nk;//,size_ind,nobc,numbc;
	int i,j,k,*ptrsize,indexP,index[3],slabIndex=0;
	double p,q,r,Area,qo,res,diffCoef,convcri=0.0;
	double *vols,*beta,*surf,*flowsol,*dummy;
	double *slabVelocity[6],*slabVolume[6],Fdiff[6];
	double digonal,sum=0.0,Diff;
	
	ptrsize=z->Size;
	ni=*(ptrsize+0);
	nj=*(ptrsize+1);
	nk=*(ptrsize+2); 
	
	qo=inwardFlux;
	vols= z->ZGeom->Volume;
	beta= z->ZGeom->Beta;
	surf= z->ZGeom->SurfaceArea;
	flowsol = MakeDomainArray(z,vartag);
	dummy   = FindSlabPtrIface(facePtr,vartag);
	MakeSlabPtrArray(z,slabVelocity,vartag);//Called for making slabs of 'Variable' corresponding to vartag
	MakeSlabPtrArray(z,slabVolume,5);//Called for making slabs of Volumes
	
	convcri = find_convcri(vartag,base);
	switch(facePtr->Direction)
	{
		case 1://east
			i=ni-2;
			do
			{
				slabIndex=0;
				sum=0.0;
				for (k=0;k<nk-1;k++) 
				for (j=0;j<nj-1;j++)
				{
					indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
					diffCoef = CalcDiffusionCoefCommon(base,z,vartag,indexP);
					p=*(surf+(indexP)*6*3+3*1+0);
					q=*(surf+(indexP)*6*3+3*1+1);
					r=*(surf+(indexP)*6*3+3*1+2);
					Area=sqrt(p*p+q*q+r*r);
					index[0]=k,index[1]=j,index[2]=ni-2;
					Diff=diffCoef*DiffusionCommon(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff);
					digonal=*(beta+(indexP)*6*3+3*1+0);
					if(digonal <= 0.0)
					{
						digonal=-digonal;
						res=(qo*Area)+diffCoef*Fdiff[1];//Diff;
					}
					else
						res=(qo*Area)-diffCoef*Fdiff[1];//Diff;
					
					*(dummy+slabIndex)=*(dummy+slabIndex)+res/(diffCoef*digonal);
					sum=sum+res*res;
					//printf("\nsum inside the loop=%12.8e  and case %d",res,direction);
					slabIndex++;
				}
				sum=sqrt(sum/((nk-1)*(nj-1)));
			}while(sum > convcri);
			break;
			
		case 2://west
			i=0;
			do
			{
				slabIndex=0;
				sum=0.0;
				for (k=0;k<nk-1;k++) 
				for (j=0;j<nj-1;j++)
				{
					indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
					diffCoef = CalcDiffusionCoefCommon(base,z,vartag,indexP);
					p=*(surf+(indexP)*6*3+3*0+0);
					q=*(surf+(indexP)*6*3+3*0+1);
					r=*(surf+(indexP)*6*3+3*0+2);
					Area=sqrt(p*p+q*q+r*r);
					index[0]=k,index[1]=j,index[2]=0;
					Diff=diffCoef*DiffusionCommon(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff);
					digonal=(*(beta+(indexP)*6*3+3*0+0));
					if(digonal <= 0.0)
					{
						digonal=-digonal;
						res=qo*Area-diffCoef*Fdiff[0];
					}
					else
						res=qo*Area+diffCoef*Fdiff[0];
						
					//res=qo*Area+diffCoef*Fdiff[0];
					//*(dummy+((nj-1)*k+j))=*(flowsol+indexP)-(qo*Area)/digonal;
					*(dummy+slabIndex)=*(dummy+slabIndex)+res/(diffCoef*digonal);
					sum=sum+res*res;
					//printf("\nsum inside the loop=%lf  and case %d",digonal,direction);
					slabIndex++;
				}
				sum=sqrt(sum/((nk-1)*(nj-1)));
			}while(sum > convcri);
			break;
		
		case 3://north
			j = nj-2;
			do
			{
				slabIndex=0;
				sum=0.0;
				for (k=0;k<nk-1;k++) 
				for (i=0;i<ni-1;i++)
				{
					indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
					diffCoef = CalcDiffusionCoefCommon(base,z,vartag,indexP);
					p = *(surf+(indexP)*6*3+3*3+0);
					q = *(surf+(indexP)*6*3+3*3+1);
					r = *(surf+(indexP)*6*3+3*3+2);
					Area=sqrt(p*p+q*q+r*r);
					index[0]=k,index[1]=nj-2,index[2]=i;
					Diff=diffCoef*DiffusionCommon(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff);
					digonal=*(beta+(indexP)*6*3+3*3+1);
					if(digonal <= 0.0)
					{
						digonal=-digonal;
						res=(qo*Area)+diffCoef*Fdiff[3];//Diff;
					}
					else
						res=(qo*Area)-diffCoef*Fdiff[3];//Diff;
					//res=-(qo*Area)-diffCoef*Fdiff[3];//Diff;
					*(dummy+slabIndex)=*(dummy+slabIndex)+res/(diffCoef*digonal);
					sum=sum+res*res;
					//printf("\nsum inside the loop=%lf  and case %d",digonal,direction);
					slabIndex++;
				}
				sum=sqrt(sum/((nk-1)*(ni-1)));
			}while(sum > convcri);
			break;
		
		case 4://south	
			j=0;
			do
			{
				slabIndex=0;
				sum=0.0;
				for (k=0;k<nk-1;k++) 
				for (i=0;i<ni-1;i++)
				{
					indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
					diffCoef = CalcDiffusionCoefCommon(base,z,vartag,indexP);
					p=*(surf+(indexP)*6*3+3*2+0);
					q=*(surf+(indexP)*6*3+3*2+1);
					r=*(surf+(indexP)*6*3+3*2+2);
					Area=sqrt(p*p+q*q+r*r);
					index[0]=k,index[1]=0,index[2]=i;
					Diff=diffCoef*DiffusionCommon(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff);
					digonal=(*(beta+(indexP)*6*3+3*2+1));
					if(digonal <= 0.0)
					{
						digonal=-digonal;
						res=qo*Area-diffCoef*Fdiff[2];
					}
					else
						res=qo*Area+diffCoef*Fdiff[2];
					
					//*(dummy+((ni-1)*k+i))=*(flowsol+indexP)-(qo*Area)/digonal;
					*(dummy+slabIndex)=*(dummy+slabIndex)+res/(diffCoef*digonal);
					sum=sum+res*res;
					//printf("\nsum inside the loop=%lf  and case %d",Fdiff[2],direction);
					slabIndex++;
				}
				sum=sqrt(sum/((nk-1)*(ni-1)));
			}while(sum > convcri);
			break;
		
		case 5://top
			k=nk-2;
			do
			{
				slabIndex=0;
				sum=0.0;
				for (j=0;j<nj-1;j++) 
				for (i=0;i<ni-1;i++)
				{
					indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
					diffCoef = CalcDiffusionCoefCommon(base,z,vartag,indexP);
					p=*(surf+(indexP)*6*3+3*5+0);
					q=*(surf+(indexP)*6*3+3*5+1);
					r=*(surf+(indexP)*6*3+3*5+2);
					Area=sqrt(p*p+q*q+r*r);
					index[0]=nk-2,index[1]=j,index[2]=i;
					Diff=diffCoef*DiffusionCommon(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff);
					digonal=*(beta+(indexP)*6*3+3*5+2);
					if(digonal <= 0.0)
					{
						digonal=-digonal;
						res=(qo*Area)+diffCoef*Fdiff[5];//Diff;
					}
					else
						res=(qo*Area)-diffCoef*Fdiff[5];//Diff;
						
					//res=(qo*Area)-diffCoef*Fdiff[5];
					*(dummy+slabIndex)=*(dummy+slabIndex)+res/(diffCoef*digonal);
					sum=sum+res*res;
					//printf("\nsum inside the loop=%lf  and case %d",digonal,direction);
					slabIndex++;
				}
				sum=sqrt(sum/((ni-1)*(nj-1)));
			}while(sum > convcri);
			break;
		
		case 6://bottom
			k=0;
			do
			{
				slabIndex=0;
				sum=0.0;
				for (j=0;j<nj-1;j++) 
				for (i=0;i<ni-1;i++)
				{
					indexP = (ni-1)*(nj-1)*k+(ni-1)*j+i;
					diffCoef = CalcDiffusionCoefCommon(base,z,vartag,indexP);
					p=*(surf+(indexP)*6*3+3*4+0);
					q=*(surf+(indexP)*6*3+3*4+1);
					r=*(surf+(indexP)*6*3+3*4+2);
					Area=sqrt(p*p+q*q+r*r);
					index[0]=0,index[1]=j,index[2]=i;
					Diff=diffCoef*DiffusionCommon(index,flowsol,vols,ptrsize,slabVelocity,slabVolume,beta,Fdiff);
					digonal=(*(beta+(indexP)*6*3+3*4+2));
					if(digonal <= 0.0)
					{
						digonal=-digonal;
						res=qo*Area-diffCoef*Fdiff[4];
					}
					else
						res=qo*Area+diffCoef*Fdiff[4];
						
					//res=qo*Area+diffCoef*Fdiff[4];
					*(dummy+slabIndex)=*(dummy+slabIndex)+res/(diffCoef*digonal);
					sum=sum+res*res;
					//printf("\nsum inside the loop=%lf  and case %d",digonal,direction);
					slabIndex++;
				}
				sum=sqrt(sum/((ni-1)*(nj-1)));
			}while(sum > convcri);
			break;
	}
	
	facePtr->ConjugateRms = sum;
	return(sum);
}

