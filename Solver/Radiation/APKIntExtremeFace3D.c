#include "APKIntExtremeFace3D.h"

/**
This function computes the diagonal coefficient and deffered correction term for the exterior face of the zone.
*/ 

double IntExtremeFace3D(Base *base, Zone *zone, double *Intface, double *NewInt, double *OldInt, double *slabInt[6], double *slabvolume[6], int index[5], double alpha, int face)
	{
	 int i,j,k,l,m,ni,nj,nk,nt,np, *ptrSize,cellindex;
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
			 volD=*(slabvolume[0]+(nj-1)*k+j);
			 IntU=*(OldInt+((ni-1)*(nj-1)*k+(ni-1)*j+i)*nt*np+np*l+m);
			 IntD=*(slabInt[0]+((nj-1)*k+j)*nt*np+np*l+m);
			 IntCU=*(NewInt+((ni-1)*(nj-1)*k+(ni-1)*j+i)*nt*np+np*l+m);
			 IntC=(volU*IntD+volD*IntU)/(volU+volD);
				 
			 *(Intface+cellindex*6*nt*np+nt*np*face+np*l+m)=((IntC-IntU)*alpha+IntCU)*(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m));
			 Diag[face]=(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m));
			}
		else
			{
			 volU=*(slabvolume[0]+(nj-1)*k+j);
			 volD=*(volume+(ni-1)*(nj-1)*k+(ni-1)*j+i);
			 IntU=*(slabInt[0]+((nj-1)*k+j)*nt*np+np*l+m);
			 IntD=*(OldInt+((ni-1)*(nj-1)*k+(ni-1)*j+i)*nt*np+np*l+m);
			 IntCU=*(slabInt[0]+((nj-1)*k+j)*nt*np+np*l+m);
			 IntC=(volU*IntD+volD*IntU)/(volU+volD);
			 
			 *(Intface+cellindex*6*nt*np+nt*np*face+np*l+m)=((IntC-IntU)*alpha+IntCU)*(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m));
			 Diag[face]=0.0;
			}
		break;
			
		case 0: //West face
		 if(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m)>0.0)
		 	{
			 volU=*(volume+(ni-1)*(nj-1)*k+(ni-1)*j+i);
			 volD=*(slabvolume[1]+(nj-1)*k+j);
			 IntU=*(OldInt+((ni-1)*(nj-1)*k+(ni-1)*j+i)*nt*np+np*l+m);
			 IntD=*(slabInt[1]+((nj-1)*k+j)*nt*np+np*l+m);
			 IntCU=*(NewInt+((ni-1)*(nj-1)*k+(ni-1)*j+i)*nt*np+np*l+m);
			 IntC=(volU*IntD+volD*IntU)/(volU+volD);
			 
			 *(Intface+cellindex*6*nt*np+nt*np*face+np*l+m)=((IntC-IntU)*alpha+IntCU)*(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m));
			 Diag[face]=(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m));
			}
		else
			{
			 volU=*(slabvolume[1]+(nj-1)*k+j);
			 volD=*(volume+(ni-1)*(nj-1)*k+(ni-1)*j+i);
			 IntU=*(slabInt[1]+((nj-1)*k+j)*nt*np+np*l+m);
			 IntD=*(OldInt+((ni-1)*(nj-1)*k+(ni-1)*j+i)*nt*np+np*l+m);
			 IntCU=*(slabInt[1]+((nj-1)*k+j)*nt*np+np*l+m);
			 IntC=(volU*IntD+volD*IntU)/(volU+volD);
				 
			 *(Intface+cellindex*6*nt*np+nt*np*face+np*l+m)=((IntC-IntU)*alpha+IntCU)*(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m));
			 Diag[face]=0.0;
			}
		break;
		
		case 3: //North face
		 if(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m)>0.0)
		 	{
			 volU=*(volume+(ni-1)*(nj-1)*k+(ni-1)*j+i);
			 volD=*(slabvolume[2]+(ni-1)*k+i);
			 IntU=*(OldInt+((ni-1)*(nj-1)*k+(ni-1)*j+i)*nt*np+np*l+m);
			 IntD=*(slabInt[2]+((ni-1)*k+i)*nt*np+np*l+m);
			 IntCU=*(NewInt+((ni-1)*(nj-1)*k+(ni-1)*j+i)*nt*np+np*l+m);
			 IntC=(volU*IntD+volD*IntU)/(volU+volD);
				 
			 *(Intface+cellindex*6*nt*np+nt*np*face+np*l+m)=((IntC-IntU)*alpha+IntCU)*(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m));
			 Diag[face]=(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m));
			}
		else
			{
			 volU=*(slabvolume[2]+(ni-1)*k+i);
			 volD=*(volume+(ni-1)*(nj-1)*k+(ni-1)*j+i);
			 IntU=*(slabInt[2]+((ni-1)*k+i)*nt*np+np*l+m);
			 IntD=*(OldInt+((ni-1)*(nj-1)*k+(ni-1)*j+i)*nt*np+np*l+m);
			 IntCU=*(slabInt[2]+((ni-1)*k+i)*nt*np+np*l+m);
			 IntC=(volU*IntD+volD*IntU)/(volU+volD);
				 
			 *(Intface+cellindex*6*nt*np+nt*np*face+np*l+m)=((IntC-IntU)*alpha+IntCU)*(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m));
			 Diag[face]=0.0;
			}
		break;
		
		case 2: //South face
		 if(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m)>0.0)
		 	{
			 volU=*(volume+(ni-1)*(nj-1)*k+(ni-1)*j+i);
			 volD=*(slabvolume[3]+(ni-1)*k+i);
			 IntU=*(OldInt+((ni-1)*(nj-1)*k+(ni-1)*j+i)*nt*np+np*l+m);
			 IntD=*(slabInt[3]+((ni-1)*k+i)*nt*np+np*l+m);
			 IntCU=*(NewInt+((ni-1)*(nj-1)*k+(ni-1)*j+i)*nt*np+np*l+m);
			 IntC=(volU*IntD+volD*IntU)/(volU+volD);
			 
			 *(Intface+cellindex*6*nt*np+nt*np*face+np*l+m)=((IntC-IntU)*alpha+IntCU)*(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m));
			 Diag[face]=(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m));
			}
		else
			{
			 volU=*(slabvolume[3]+(ni-1)*k+i);
			 volD=*(volume+(ni-1)*(nj-1)*k+(ni-1)*j+i);
			 IntU=*(slabInt[3]+((ni-1)*k+i)*nt*np+np*l+m);
			 IntD=*(OldInt+((ni-1)*(nj-1)*k+(ni-1)*j+i)*nt*np+np*l+m);
			 IntCU=*(slabInt[3]+((ni-1)*k+i)*nt*np+np*l+m);
			 IntC=(volU*IntD+volD*IntU)/(volU+volD);
				 
			 *(Intface+cellindex*6*nt*np+nt*np*face+np*l+m)=((IntC-IntU)*alpha+IntCU)*(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m));
			 Diag[face]=0.0;
			}
		break;
			
		case 5: //top face
		 if(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m)>0.0)
		 	{
			 volU=*(volume+(ni-1)*(nj-1)*k+(ni-1)*j+i);
			 volD=*(slabvolume[4]+(ni-1)*j+i);
			 IntU=*(OldInt+((ni-1)*(nj-1)*k+(ni-1)*j+i)*nt*np+np*l+m);
			 IntD=*(slabInt[4]+((ni-1)*j+i)*nt*np+np*l+m);
			 IntCU=*(NewInt+((ni-1)*(nj-1)*k+(ni-1)*j+i)*nt*np+np*l+m);
			 IntC=(volU*IntD+volD*IntU)/(volU+volD);
				 
			 *(Intface+cellindex*6*nt*np+nt*np*face+np*l+m)=((IntC-IntU)*alpha+IntCU)*(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m));
			 Diag[face]=(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m));
			}
		else
			{
			 volU=*(slabvolume[4]+(ni-1)*j+i);
			 volD=*(volume+(ni-1)*(nj-1)*k+(ni-1)*j+i);
			 IntU=*(slabInt[4]+((ni-1)*j+i)*nt*np+np*l+m);
			 IntD=*(OldInt+((ni-1)*(nj-1)*k+(ni-1)*j+i)*nt*np+np*l+m);
			 IntCU=*(slabInt[4]+((ni-1)*j+i)*nt*np+np*l+m);
			 IntC=(volU*IntD+volD*IntU)/(volU+volD);
			 
			 *(Intface+cellindex*6*nt*np+nt*np*face+np*l+m)=((IntC-IntU)*alpha+IntCU)*(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m));
			 Diag[face]=0.0;
			}
		break;
			
		case 4: //Bottom face
		 if(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m)>0.0)
		 	{
			 volU=*(volume+(ni-1)*(nj-1)*k+(ni-1)*j+i);
			 volD=*(slabvolume[5]+(ni-1)*j+i);
			 IntU=*(OldInt+((ni-1)*(nj-1)*k+(ni-1)*j+i)*nt*np+np*l+m);
			 IntD=*(slabInt[5]+((ni-1)*j+i)*nt*np+np*l+m);
			 IntCU=*(NewInt+((ni-1)*(nj-1)*k+(ni-1)*j+i)*nt*np+np*l+m);
			 IntC=(volU*IntD+volD*IntU)/(volU+volD);
			 
			 *(Intface+cellindex*6*nt*np+nt*np*face+np*l+m)=((IntC-IntU)*alpha+IntCU)*(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m));
			 Diag[face]=(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m));
			}
		else
			{
			 volU=*(slabvolume[5]+(ni-1)*j+i);
			 volD=*(volume+(ni-1)*(nj-1)*k+(ni-1)*j+i);
			 IntU=*(slabInt[5]+((ni-1)*j+i)*nt*np+np*l+m);
			 IntD=*(OldInt+((ni-1)*(nj-1)*k+(ni-1)*j+i)*nt*np+np*l+m);
			 IntCU=*(slabInt[5]+((ni-1)*j+i)*nt*np+np*l+m);
			 IntC=(volU*IntD+volD*IntU)/(volU+volD);
			 
			 *(Intface+cellindex*6*nt*np+nt*np*face+np*l+m)=((IntC-IntU)*alpha+IntCU)*(*(Ds+cellindex*6*nt*np+nt*np*face+np*l+m));
			 Diag[face]=0.0;
			}
		break;
		}
		
	return (Diag[face]);
	}
