#include "APKIntExtremeCell3D.h"

/**
This function computes the Deferred correction term for 3 dimensional geometries for those cells which have more than one boundary or 
interface among its six faces. This function returns both the diagonal coefficient and the Deferred correction term
(which is calculated using the previous iteration values). */ 

double IntExtremeCell3D(Base *base, Zone *zone, double *Intface, double *NewInt, double *OldInt,  int index[5], double alpha, int face)
	{
	 int ni,nj,nk,nt,np,i,j,k,l,m,*ptrSize,cellindex;
	 double *Ds, volU, volD, IntU, IntD, IntC, IntCU, Diag[6], *volume;
	  
	 i=index[0],j=index[1],k=index[2],l=index[3],m=index[4];
	 ptrSize=zone->Size;
	 ni=*(ptrSize+0);
	 nj=*(ptrSize+1);
	 nk=*(ptrSize+2);
	 nt=base->BaseParameters->ntheta;
	 np=base->BaseParameters->nphi;
	 volume=zone->ZGeom->Volume;
	 cellindex=(ni-1)*(nj-1)*k+(ni-1)*j+i;
	 
	 Ds=zone->ZGeom->DmF;
	 
	 	 switch(face)
		 	{
			 case 1: //East Face
			 if(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m)>0.0)
			 	{
				 volU=*(volume+(ni-1)*(nj-1)*k+(ni-1)*j+i);
				 volD=*(volume+(ni-1)*(nj-1)*k+(ni-1)*j+i+1);
				 IntU=*(OldInt+((ni-1)*(nj-1)*k+(ni-1)*j+i)*nt*np+np*l+m);
				 IntD=*(OldInt+((ni-1)*(nj-1)*k+(ni-1)*j+i+1)*nt*np+np*l+m);
				 IntCU=*(NewInt+((ni-1)*(nj-1)*k+(ni-1)*j+i)*nt*np+np*l+m);
				 IntC=(volU*IntD+volD*IntU)/(volU+volD);
				 
				 *(Intface+cellindex*6*nt*np+nt*np*face+np*l+m)=((IntC-IntU)*alpha+IntCU)*(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m));
				 Diag[face]=(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m));
				}
			else
				{
				 volU=*(volume+(ni-1)*(nj-1)*k+(ni-1)*j+i+1);
				 volD=*(volume+(ni-1)*(nj-1)*k+(ni-1)*j+i);
				 IntU=*(OldInt+((ni-1)*(nj-1)*k+(ni-1)*j+i+1)*nt*np+np*l+m);
				 IntD=*(OldInt+((ni-1)*(nj-1)*k+(ni-1)*j+i)*nt*np+np*l+m);
				 IntCU=*(NewInt+((ni-1)*(nj-1)*k+(ni-1)*j+i+1)*nt*np+np*l+m);
				 IntC=(volU*IntD+volD*IntU)/(volU+volD);
				 
				 *(Intface+cellindex*6*nt*np+nt*np*face+np*l+m)=((IntC-IntU)*alpha+IntCU)*(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m));
				 Diag[face]=0.0;
				}
			break;
			
			case 0: //West face
			 if(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m)>0.0)
			 	{
				 volU=*(volume+(ni-1)*(nj-1)*k+(ni-1)*j+i);
				 volD=*(volume+(ni-1)*(nj-1)*k+(ni-1)*j+i-1);
				 IntU=*(OldInt+((ni-1)*(nj-1)*k+(ni-1)*j+i)*nt*np+np*l+m);
				 IntD=*(OldInt+((ni-1)*(nj-1)*k+(ni-1)*j+i-1)*nt*np+np*l+m);
				 IntCU=*(NewInt+((ni-1)*(nj-1)*k+(ni-1)*j+i)*nt*np+np*l+m);
				 IntC=(volU*IntD+volD*IntU)/(volU+volD);
				 
				 *(Intface+cellindex*6*nt*np+nt*np*face+np*l+m)=((IntC-IntU)*alpha+IntCU)*(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m));
				 Diag[face]=(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m));
				}
			else
				{
				 volU=*(volume+(ni-1)*(nj-1)*k+(ni-1)*j+i-1);
				 volD=*(volume+(ni-1)*(nj-1)*k+(ni-1)*j+i);
				 IntU=*(OldInt+((ni-1)*(nj-1)*k+(ni-1)*j+i-1)*nt*np+np*l+m);
				 IntD=*(OldInt+((ni-1)*(nj-1)*k+(ni-1)*j+i)*nt*np+np*l+m);
				 IntCU=*(NewInt+((ni-1)*(nj-1)*k+(ni-1)*j+i-1)*nt*np+np*l+m);
				 IntC=(volU*IntD+volD*IntU)/(volU+volD);
				 
				 *(Intface+cellindex*6*nt*np+nt*np*face+np*l+m)=((IntC-IntU)*alpha+IntCU)*(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m));
				 Diag[face]=0.0;
				}
			break;
			
			case 3: //North face
			 if(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m)>0.0)
			 	{
				 volU=*(volume+(ni-1)*(nj-1)*k+(ni-1)*j+i);
				 volD=*(volume+(ni-1)*(nj-1)*k+(ni-1)*(j+1)+i);
				 IntU=*(OldInt+((ni-1)*(nj-1)*k+(ni-1)*j+i)*nt*np+np*l+m);
				 IntD=*(OldInt+((ni-1)*(nj-1)*k+(ni-1)*(j+1)+i)*nt*np+np*l+m);
				 IntCU=*(NewInt+((ni-1)*(nj-1)*k+(ni-1)*j+i)*nt*np+np*l+m);
				 IntC=(volU*IntD+volD*IntU)/(volU+volD);
				 
				 *(Intface+cellindex*6*nt*np+nt*np*face+np*l+m)=((IntC-IntU)*alpha+IntCU)*(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m));
				 Diag[face]=(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m));
				}
			else
				{
				 volU=*(volume+(ni-1)*(nj-1)*k+(ni-1)*(j+1)+i);
				 volD=*(volume+(ni-1)*(nj-1)*k+(ni-1)*j+i);
				 IntU=*(OldInt+((ni-1)*(nj-1)*k+(ni-1)*(j+1)+i)*nt*np+np*l+m);
				 IntD=*(OldInt+((ni-1)*(nj-1)*k+(ni-1)*j+i)*nt*np+np*l+m);
				 IntCU=*(NewInt+((ni-1)*(nj-1)*k+(ni-1)*(j+1)+i)*nt*np+np*l+m);
				 IntC=(volU*IntD+volD*IntU)/(volU+volD);
				 
				 *(Intface+cellindex*6*nt*np+nt*np*face+np*l+m)=((IntC-IntU)*alpha+IntCU)*(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m));
				 Diag[face]=0.0;
				}
			break;
			
			case 2: //South face
			 if(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m)>0.0)
			 	{
				 volU=*(volume+(ni-1)*(nj-1)*k+(ni-1)*j+i);
				 volD=*(volume+(ni-1)*(nj-1)*k+(ni-1)*(j-1)+i);
				 IntU=*(OldInt+((ni-1)*(nj-1)*k+(ni-1)*j+i)*nt*np+np*l+m);
				 IntD=*(OldInt+((ni-1)*(nj-1)*k+(ni-1)*(j-1)+i)*nt*np+np*l+m);
				 IntCU=*(NewInt+((ni-1)*(nj-1)*k+(ni-1)*j+i)*nt*np+np*l+m);
				 IntC=(volU*IntD+volD*IntU)/(volU+volD);
				 
				 *(Intface+cellindex*6*nt*np+nt*np*face+np*l+m)=((IntC-IntU)*alpha+IntCU)*(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m));
				 Diag[face]=(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m));
				}
			else
				{
				 volU=*(volume+(ni-1)*(nj-1)*k+(ni-1)*(j-1)+i);
				 volD=*(volume+(ni-1)*(nj-1)*k+(ni-1)*j+i);
				 IntU=*(OldInt+((ni-1)*(nj-1)*k+(ni-1)*(j-1)+i)*nt*np+np*l+m);
				 IntD=*(OldInt+((ni-1)*(nj-1)*k+(ni-1)*j+i)*nt*np+np*l+m);
				 IntCU=*(NewInt+((ni-1)*(nj-1)*k+(ni-1)*(j-1)+i)*nt*np+np*l+m);
				 IntC=(volU*IntD+volD*IntU)/(volU+volD);
				 
				 *(Intface+cellindex*6*nt*np+nt*np*face+np*l+m)=((IntC-IntU)*alpha+IntCU)*(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m));
				 Diag[face]=0.0;
				}
			break;
			
			case 5: //top face
			 if(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m)>0.0)
			 	{
				 volU=*(volume+(ni-1)*(nj-1)*k+(ni-1)*j+i);
				 volD=*(volume+(ni-1)*(nj-1)*(k+1)+(ni-1)*j+i);
				 IntU=*(OldInt+((ni-1)*(nj-1)*k+(ni-1)*j+i)*nt*np+np*l+m);
				 IntD=*(OldInt+((ni-1)*(nj-1)*(k+1)+(ni-1)*j+i)*nt*np+np*l+m);
				 IntCU=*(NewInt+((ni-1)*(nj-1)*k+(ni-1)*j+i)*nt*np+np*l+m);
				 IntC=(volU*IntD+volD*IntU)/(volU+volD);
				 
				 *(Intface+cellindex*6*nt*np+nt*np*face+np*l+m)=((IntC-IntU)*alpha+IntCU)*(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m));
				 Diag[face]=(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m));
				}
			else
				{
				 volU=*(volume+(ni-1)*(nj-1)*(k+1)+(ni-1)*j+i);
				 volD=*(volume+(ni-1)*(nj-1)*k+(ni-1)*j+i);
				 IntU=*(OldInt+((ni-1)*(nj-1)*(k+1)+(ni-1)*j+i)*nt*np+np*l+m);
				 IntD=*(OldInt+((ni-1)*(nj-1)*k+(ni-1)*j+i)*nt*np+np*l+m);
				 IntCU=*(NewInt+((ni-1)*(nj-1)*(k+1)+(ni-1)*j+i)*nt*np+np*l+m);
				 IntC=(volU*IntD+volD*IntU)/(volU+volD);
				 
				 *(Intface+cellindex*6*nt*np+nt*np*face+np*l+m)=((IntC-IntU)*alpha+IntCU)*(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m));
				 Diag[face]=0.0;
				}
			break;
			
			case 4: //Bottom face
			 if(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m)>0.0)
			 	{
				 volU=*(volume+(ni-1)*(nj-1)*k+(ni-1)*j+i);
				 volD=*(volume+(ni-1)*(nj-1)*(k-1)+(ni-1)*j+i);
				 IntU=*(OldInt+((ni-1)*(nj-1)*k+(ni-1)*j+i)*nt*np+np*l+m);
				 IntD=*(OldInt+((ni-1)*(nj-1)*(k-1)+(ni-1)*j+i)*nt*np+np*l+m);
				 IntCU=*(NewInt+((ni-1)*(nj-1)*k+(ni-1)*j+i)*nt*np+np*l+m);
				 IntC=(volU*IntD+volD*IntU)/(volU+volD);
				 
				 *(Intface+cellindex*6*nt*np+nt*np*face+np*l+m)=((IntC-IntU)*alpha+IntCU)*(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m));
				 Diag[face]=(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m));
				}
			else
				{
				 volU=*(volume+(ni-1)*(nj-1)*(k-1)+(ni-1)*j+i);
				 volD=*(volume+(ni-1)*(nj-1)*k+(ni-1)*j+i);
				 IntU=*(OldInt+((ni-1)*(nj-1)*(k-1)+(ni-1)*j+i)*nt*np+np*l+m);
				 IntD=*(OldInt+((ni-1)*(nj-1)*k+(ni-1)*j+i)*nt*np+np*l+m);
				 IntCU=*(NewInt+((ni-1)*(nj-1)*(k-1)+(ni-1)*j+i)*nt*np+np*l+m);
				 IntC=(volU*IntD+volD*IntU)/(volU+volD);
				 
				 *(Intface+cellindex*6*nt*np+nt*np*face+np*l+m)=((IntC-IntU)*alpha+IntCU)*(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m));
				 Diag[face]=0.0;
				}
			break;
			}
	return(Diag[face]);
	}
