#include "APKDisperisonflux2d.h"
/**
This function is used to calculate dispersion flux and faces of cells and return diagonal term  for 2-D 
Similar to this for 3D , function Disperisonflux3d  is used.

*/
double Disperisonflux2d(int index[3],double *domain,double *volumes,int *size,double *slabVariable[6],double *slabVolume[6],double *diffusionCoeffdomain,double *slabDiffusionCoeff[6],double *SurfArea,double DiffFlux[6],double *GradX,double *GradY,double *slabGradX[6],double *slabGradY[6],double *DeltaI,double *slabDeltaInverse[6])
 {
	//2 D case
	int NK,NJ,NI,k,j,i,iFace,comp;
	double u_p,u_e,u_w,u_s,u_n,vol_p,vol_e,vol_w,vol_s,vol_n;
	double Area[6][3];
	double diff_p[9],diff_e[9],diff_w[9], diff_s[9],diff_n[9];//These are Cell center values
	double diff_E[9],diff_W[9],diff_N[9],diff_S[9];// face centre value
	double GradX_p,GradX_e,GradX_w,GradX_s,GradX_n;
	double GradX_E,GradX_W,GradX_N,GradX_S;
	double GradY_p,GradY_e,GradY_w,GradY_s,GradY_n;
	double GradY_E,GradY_W,GradY_N,GradY_S;
	double D_w[4],D_e[4],D_n[4],D_s[4];
	double Coeff1[6],Coeff2[6];
	double DiagCoeff[4],diagonalW,diagonalE,diagonalS,diagonalN,diagonal;//
	
	
	int c,IND,IND_e,IND_w,IND_n,IND_s;
	
	NI=*(size+0);
	NJ=*(size+1);
	NK=*(size+2);
	
	k=index[0],j=index[1],i=index[2];
	
	IND=(NJ-1)*(NI-1)*k+(NI-1)*j+i;
	
	u_p=*(domain+IND);
	vol_p=*(volumes+IND);
	
	GradX_p=*(GradX+IND);
	GradY_p=*(GradY+IND);
	
	
	for(c=0;c<9;c++)
	diff_p[c]=*(diffusionCoeffdomain+9*(IND)+c);
		
	if(i!=0)
	{        
	       IND_w=(NJ-1)*(NI-1)*k+(NI-1)*j+(i-1);
	         
		u_w=*(domain+IND_w);
		vol_w=*(volumes+IND_w);
		
		GradX_w=*(GradX+IND_w);
		GradY_w=*(GradY+IND_w);
		
		for(c=0;c<4;c++)
		D_w[c]=*(DeltaI+4*IND_w+c);
		
		for(c=0;c<9;c++)
		diff_w[c]=*(diffusionCoeffdomain+9*IND_w+c);		
	}
	else
	{         IND_w=(NJ-1)*k+j;
	           
		u_w =*(slabVariable[1]+IND_w);
		vol_w =*(slabVolume[1]+IND_w);
		
		GradX_w=*(slabGradX[1]+IND_w);
		GradY_w=*(slabGradY[1]+IND_w);
		
		for(c=0;c<4;c++)
		D_w[c]=*(slabDeltaInverse[1]+4*IND_w+c);
		
					
		for(c=0;c<9;c++)
	       diff_w[c]=*(slabDiffusionCoeff[1]+9*IND_w+c);
	}
	if(i!=(NI-2))
	{
		IND_e= (NJ-1)*(NI-1)*k+(NI-1)*j+(i+1);
		u_e=*(domain+IND_e);
		vol_e=*(volumes+IND_e);
		
		GradX_e=*(GradX+IND_e);
		GradY_e=*(GradY+IND_e);
		
		for(c=0;c<4;c++)
		D_e[c]=*(DeltaI+4*IND_e+c);
		
		for(c=0;c<9;c++)
		diff_e[c]=*(diffusionCoeffdomain+9*(IND_e)+c);
	}
	else
	{        IND_e=(NJ-1)*k+j;
		u_e= *(slabVariable[0]+IND_e);
		vol_e=*(slabVolume[0]+IND_e);
	        
		GradX_e=*(slabGradX[0]+IND_e);
		GradY_e=*(slabGradY[0]+IND_e);
		
		for(c=0;c<4;c++)
		D_e[c]=*(slabDeltaInverse[0]+4*IND_e+c);
		 
		 for(c=0;c<9;c++)
	        diff_e[c]=*(slabDiffusionCoeff[0]+9*(IND_e)+c); 
	}   
	if(j!=0)
	{ 	IND_s=(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i;
		u_s=*(domain+IND_s);
		vol_s=*(volumes+IND_s);
		
		GradX_s=*(GradX+IND_s);
		GradY_s=*(GradY+IND_s);
		
		for(c=0;c<4;c++)
		D_s[c]=*(DeltaI+4*IND_s+c);
		
		for(c=0;c<9;c++)
		diff_s[c]=*(diffusionCoeffdomain+9*(IND_s)+c);
	}
	else
	{	IND_s=(NI-1)*k+i;
		u_s=*(slabVariable[3]+IND_s);
		vol_s=*(slabVolume[3]+IND_s);
	
		GradX_s=*(slabGradX[3]+IND_s);
		GradY_s=*(slabGradY[3]+IND_s);
		
		for(c=0;c<4;c++)
		D_s[c]=*(slabDeltaInverse[3]+4*IND_s+c);
		
		for(c=0;c<9;c++)
		{diff_s[c]=*(slabDiffusionCoeff[3]+9*IND_s+c);
		 
		}
		
		
	}
	if(j!=(NJ-2))
	{	IND_n=(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i;
		u_n=*(domain+IND_n);
		vol_n=*(volumes+IND_n);
		
		GradX_n=*(GradX+IND_n);
		GradY_n=*(GradY+IND_n);
		
		for(c=0;c<4;c++)
		D_n[c]=*(DeltaI+4*IND_n+c);
		
		for(c=0;c<9;c++)
		diff_n[c]=*(diffusionCoeffdomain+9*(IND_n)+c);
	}
	else
	{	IND_n=(NI-1)*k+i;
		u_n=*(slabVariable[2]+IND_n);
		vol_n=*(slabVolume[2]+IND_n);
		
		GradX_n=*(slabGradX[2]+IND_n);
		GradY_n=*(slabGradY[2]+IND_n);
		
		for(c=0;c<4;c++)
		D_n[c]=*(slabDeltaInverse[2]+4*IND_n+c);
		
		for(c=0;c<9;c++)
		diff_n[c]=*(slabDiffusionCoeff[2]+9*(IND_n)+c);
	}
	
	
	
	DiagCoeff[0]=vol_p/(vol_p+vol_w);
	DiagCoeff[1]=vol_p/(vol_p+vol_e);
	DiagCoeff[2]=vol_p/(vol_p+vol_s);
	DiagCoeff[3]=vol_p/(vol_p+vol_n);
	
	
	
	for(c=0;c<9;c++)
	{
	diff_E[c] = (diff_e[c]*vol_p+diff_p[c]*vol_e)/(vol_p+vol_e);
	diff_W[c] = (diff_w[c]*vol_p+diff_p[c]*vol_w)/(vol_p+vol_w);
	diff_N[c] = (diff_n[c]*vol_p+diff_p[c]*vol_n)/(vol_p+vol_n);
	diff_S[c] = (diff_s[c]*vol_p+diff_p[c]*vol_s)/(vol_p+vol_s);
//	diff_T[c] = 0.0; //(diff_t[c]*vol_p + diff_p[c]*vol_t)/(vol_p+vol_t);
//	diff_B[c] = 0.0; //(diff_b[c]*vol_p + diff_p[c]*vol_b)/(vol_p+vol_b);
	}
	     
	      
	     
	 
	      
	
	
	GradX_E = (GradX_e*vol_p+GradX_p*vol_e)/(vol_p+vol_e);
	GradX_W = (GradX_w*vol_p+GradX_p*vol_w)/(vol_p+vol_w);
	GradX_N = (GradX_n*vol_p+GradX_p*vol_n)/(vol_p+vol_n);
	GradX_S = (GradX_s*vol_p+GradX_p*vol_s)/(vol_p+vol_s);
	 
	GradY_E = (GradY_e*vol_p+GradY_p*vol_e)/(vol_p+vol_e);
	GradY_W = (GradY_w*vol_p+GradY_p*vol_w)/(vol_p+vol_w);
	GradY_N = (GradY_n*vol_p+GradY_p*vol_n)/(vol_p+vol_n);
	GradY_S = (GradY_s*vol_p+GradY_p*vol_s)/(vol_p+vol_s);
		
	
	        for(iFace=0;iFace<6;iFace++)
		for(comp=0;comp<3;comp++)
		Area[iFace][comp]=*(SurfArea+IND*6*3+3*iFace+comp);
		
		Coeff1[0]= diff_W[0]*Area[0][0]+diff_W[3]*Area[0][1]+diff_W[6]* Area[0][2];//west face 
		Coeff1[1]= diff_E[0]*Area[1][0]+diff_E[3]*Area[1][1]+diff_E[6]* Area[1][2];//east face 
		Coeff1[2]= diff_S[0]*Area[2][0]+diff_S[3]*Area[2][1]+diff_S[6]* Area[2][2];//south face
		Coeff1[3]= diff_N[0]*Area[3][0]+diff_N[3]*Area[3][1]+diff_N[6]* Area[3][2];//north face
		//Coeff1[4]= diff_B[0]*Area[4][0]+diff_B[3]*Area[4][1]+diff_B[6]* Area[4][2];
		//Coeff1[5]= diff_T[0]*Area[5][0]+diff_T[3]*Area[4][1]+diff_T[6]* Area[5][2];
		
		
		Coeff2[0]= diff_W[1]*Area[0][0]+diff_W[4]*Area[0][1]+diff_W[7]* Area[0][2];//west face 
		Coeff2[1]= diff_E[1]*Area[1][0]+diff_E[4]*Area[1][1]+diff_E[7]* Area[1][2];//east face 
		Coeff2[2]= diff_S[1]*Area[2][0]+diff_S[4]*Area[2][1]+diff_S[7]* Area[2][2];//south face
		Coeff2[3]= diff_N[1]*Area[3][0]+diff_N[4]*Area[3][1]+diff_N[7]* Area[3][2];//north face
		//Coeff2[4]= diff_B[1]*Area[4][0]+diff_B[4]*Area[4][1]+diff_B[7]* Area[4][2];
		//Coeff2[5]= diff_T[1]*Area[5][0]+diff_T[4]*Area[4][1]+diff_T[7]* Area[5][2];
		
	
		//for 3d problem these are required for 2D not, however including doenot affetc any thing
		//Coeff3[0]= diff_W[2]*Area[0][0]+diff_W[5]*Area[0][1]+diff_W[8]* Area[0][2];//west face 
		//Coeff3[1]= diff_E[2]*Area[1][0]+diff_E[5]*Area[1][1]+diff_E[8]* Area[1][2];//east face 
		//Coeff3[2]= diff_S[2]*Area[2][0]+diff_S[5]*Area[2][1]+diff_S[8]* Area[2][2];//south face
		//Coeff3[3]= diff_N[2]*Area[3][0]+diff_N[5]*Area[3][1]+diff_N[8]* Area[3][2];//north face
		//Coeff3[4]= diff_B[2]*Area[4][0]+diff_B[5]*Area[4][5]+diff_B[8]* Area[4][2];
		//Coeff3[5]= diff_T[2]*Area[5][0]+diff_T[5]*Area[4][1]+diff_T[8]* Area[5][2];
		
	 DiffFlux[0]=Coeff1[0]*GradX_W+Coeff2[0]*GradY_W;//west face 
			 
     	 DiffFlux[1]=Coeff1[1]*GradX_E+Coeff2[1]*GradY_E;//east face 
      
     	 DiffFlux[2]=Coeff1[2]*GradX_S+Coeff2[2]*GradY_S;//south face
	
     	 DiffFlux[3]=Coeff1[3]*GradX_N+Coeff2[3]*GradY_N;//north face
	 
	 DiffFlux[4]=0.0;
	 
	 DiffFlux[5]=0.0;
	// pi[0]=UE-UW;
	//	pi[1]=UN-US;
	//D[0]=*((zone->ZGeom->DeltaInverse+P*4)+0) * pi[0] + *((zone->ZGeom->DeltaInverse+P*4)+1) * pi[1];
	//D[1]=*((zone->ZGeom->DeltaInverse+P*4)+2) * pi[0] + *((zone->ZGeom->DeltaInverse+P*4)+3) * pi[1];
	 
	diagonalW=  DiagCoeff[0]*(D_w[0]*Coeff1[0]+D_w[2]*Coeff2[0]);
	diagonalE= -DiagCoeff[1]*(D_e[0]*Coeff1[1]+D_e[2]*Coeff2[1]);
	diagonalS=  DiagCoeff[2]*(D_s[1]*Coeff1[2]+D_s[3]*Coeff2[2]);
	diagonalN= -DiagCoeff[3]*(D_n[1]*Coeff1[3]+D_n[3]*Coeff2[3]); 
	
	diagonal=diagonalW+diagonalE+diagonalS+diagonalN;
	//printf("\n %lf",diagonal);
//printf("\nDiff -%12.8lf  %12.8lf %12.8lf %12.8lf  %12.8lf\n",DiffFlux[0],DiffFlux[1],DiffFlux[2],DiffFlux[3],diagonal);
	//exit(1);
	return (diagonal);
}


double Disperisonflux3d(int index[3],double *domain,double *volumes,int *size,double *slabVariable[6],double *slabVolume[6],double *diffusionCoeffdomain,double *slabDiffusionCoeff[6],double *SurfArea,double DiffFlux[6],double *GradX,double *GradY,double *GradZ,double *slabGradX[6],double *slabGradY[6],double *slabGradZ[6],double *DeltaI,double *slabDeltaInverse[6])
 {
	
	int NK,NJ,NI,k,j,i,iFace,comp;
	double u_p,u_e,u_w,u_s,u_n,u_t,u_b,vol_p,vol_e,vol_w,vol_s,vol_n,vol_t,vol_b;
	
	double Area[6][3];
	
	
	double diff_p[9],diff_e[9],diff_w[9], diff_s[9],diff_n[9],diff_t[9],diff_b[9];//These are Cell center values
	double diff_E[9],diff_W[9],diff_N[9],diff_S[9],diff_T[9],diff_B[9];// face centre value
	double GradX_p,GradX_e,GradX_w,GradX_s,GradX_n,GradX_t,GradX_b;
	double GradX_E,GradX_W,GradX_N,GradX_S,GradX_T,GradX_B;
	double GradY_p,GradY_e,GradY_w,GradY_s,GradY_n,GradY_t,GradY_b;
	double GradY_E,GradY_W,GradY_N,GradY_S,GradY_T,GradY_B;
	double GradZ_p,GradZ_e,GradZ_w,GradZ_s,GradZ_n,GradZ_t,GradZ_b;
	double GradZ_E,GradZ_W,GradZ_N,GradZ_S,GradZ_T,GradZ_B;
	
	double D_w[9],D_e[9],D_n[9],D_s[9],D_t[9],D_b[9];
	double Coeff1[6],Coeff2[6],Coeff3[6];
	double DiagCoeff[6],diagonalW,diagonalE,diagonalS,diagonalN,diagonalT,diagonalB,diagonal;
	
	
	int c,IND,IND_e,IND_w,IND_n,IND_s,IND_b,IND_t;
	
	NI=*(size+0);
	NJ=*(size+1);
	NK=*(size+2);
	
	k=index[0],j=index[1],i=index[2];
	
	IND=(NJ-1)*(NI-1)*k+(NI-1)*j+i;
	
	u_p=*(domain+IND);
	vol_p=*(volumes+IND);
	
	GradX_p=*(GradX+IND);
	GradY_p=*(GradY+IND);
	GradZ_p=*(GradZ+IND);
	
	for(c=0;c<9;c++)
	diff_p[c]=*(diffusionCoeffdomain+9*(IND)+c);
		
	if(i!=0)
	{        
	       IND_w=(NJ-1)*(NI-1)*k+(NI-1)*j+(i-1);
	         
		u_w=*(domain+IND_w);
		vol_w=*(volumes+IND_w);
		
		GradX_w=*(GradX+IND_w);
		GradY_w=*(GradY+IND_w);
		GradZ_w=*(GradZ+IND_w);
		
		
		for(c=0;c<9;c++)
		D_w[c]=*(DeltaI+9*IND_w+c);
		
		for(c=0;c<9;c++)
		diff_w[c]=*(diffusionCoeffdomain+9*IND_w+c);		
	}
	else
	{         IND_w=(NJ-1)*k+j;
	           
		u_w =*(slabVariable[1]+IND_w);
		vol_w =*(slabVolume[1]+IND_w);
		
		GradX_w=*(slabGradX[1]+IND_w);
		GradY_w=*(slabGradY[1]+IND_w);
		GradZ_w=*(slabGradZ[1]+IND_w);
		
		
		for(c=0;c<9;c++)
		D_w[c]=*(slabDeltaInverse[1]+9*IND_w+c);
		
					
		for(c=0;c<9;c++)
	       diff_w[c]=*(slabDiffusionCoeff[1]+9*IND_w+c);
	}
	if(i!=(NI-2))
	{
		IND_e= (NJ-1)*(NI-1)*k+(NI-1)*j+(i+1);
		u_e=*(domain+IND_e);
		vol_e=*(volumes+IND_e);
		
		GradX_e=*(GradX+IND_e);
		GradY_e=*(GradY+IND_e);
		GradZ_e=*(GradZ+IND_e);
		
		for(c=0;c<9;c++)
		D_e[c]=*(DeltaI+9*IND_e+c);
		
		for(c=0;c<9;c++)
		diff_e[c]=*(diffusionCoeffdomain+9*(IND_e)+c);
	}
	else
	{        IND_e=(NJ-1)*k+j;
		u_e= *(slabVariable[0]+IND_e);
		vol_e=*(slabVolume[0]+IND_e);
	        
		GradX_e=*(slabGradX[0]+IND_e);
		GradY_e=*(slabGradY[0]+IND_e);
		GradZ_e=*(slabGradZ[0]+IND_e);
		
		
		for(c=0;c<9;c++)
		D_e[c]=*(slabDeltaInverse[0]+9*IND_e+c);
		 
		 for(c=0;c<9;c++)
	        diff_e[c]=*(slabDiffusionCoeff[0]+9*(IND_e)+c); 
	}   
	if(j!=0)
	{ 	IND_s=(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i;
		u_s=*(domain+IND_s);
		vol_s=*(volumes+IND_s);
		
		GradX_s=*(GradX+IND_s);
		GradY_s=*(GradY+IND_s);
		GradZ_s=*(GradZ+IND_s);
		
		
		for(c=0;c<9;c++)
		D_s[c]=*(DeltaI+9*IND_s+c);
		
		for(c=0;c<9;c++)
		diff_s[c]=*(diffusionCoeffdomain+9*(IND_s)+c);
	}
	else
	{	IND_s=(NI-1)*k+i;
		u_s=*(slabVariable[3]+IND_s);
		vol_s=*(slabVolume[3]+IND_s);
	
		GradX_s=*(slabGradX[3]+IND_s);
		GradY_s=*(slabGradY[3]+IND_s);
		GradZ_s=*(slabGradZ[3]+IND_s);
		
		
		for(c=0;c<9;c++)
		D_s[c]=*(slabDeltaInverse[3]+9*IND_s+c);
		
		for(c=0;c<9;c++)
		{diff_s[c]=*(slabDiffusionCoeff[3]+9*IND_s+c);
		 
		}
		
		
	}
	if(j!=(NJ-2))
	{	IND_n=(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i;
		u_n=*(domain+IND_n);
		vol_n=*(volumes+IND_n);
		
		GradX_n=*(GradX+IND_n);
		GradY_n=*(GradY+IND_n);
		GradZ_n=*(GradZ+IND_n);
		
		
		for(c=0;c<9;c++)
		D_n[c]=*(DeltaI+9*IND_n+c);
		
		for(c=0;c<9;c++)
		diff_n[c]=*(diffusionCoeffdomain+9*(IND_n)+c);
	}
	else
	{	IND_n=(NI-1)*k+i;
		u_n=*(slabVariable[2]+IND_n);
		vol_n=*(slabVolume[2]+IND_n);
		
		GradX_n=*(slabGradX[2]+IND_n);
		GradY_n=*(slabGradY[2]+IND_n);
		GradZ_n=*(slabGradZ[2]+IND_n);
		
		for(c=0;c<9;c++)
		D_n[c]=*(slabDeltaInverse[2]+9*IND_n+c);
		
		for(c=0;c<9;c++)
		diff_n[c]=*(slabDiffusionCoeff[2]+9*(IND_n)+c);
	}
	if(k!=0)
	{        IND_b=(NJ-1)*(NI-1)*(k-1)+(NI-1)*j+i;
	        
		u_b=*(domain+IND_b);
		vol_b=*(volumes+IND_b);
		
		GradX_b=*(GradX+IND_b);
		GradY_b=*(GradY+IND_b);
		GradZ_b=*(GradZ+IND_b);
		
		for(c=0;c<9;c++)
		D_b[c]=*(DeltaI+9*IND_b+c);
		
		for(c=0;c<9;c++)
		diff_b[c]=*(diffusionCoeffdomain+9*IND_b+c);
	}
	else
	{	IND_b=(NI-1)*j+i;
		u_b=*(slabVariable[5]+IND_b);
		vol_b=*(slabVolume[5]+IND_b);
		
		GradX_b=*(slabGradX[5]+IND_b);
		GradY_b=*(slabGradY[5]+IND_b);
		GradZ_b=*(slabGradZ[5]+IND_b);
		
		for(c=0;c<9;c++)
		D_b[c]=*(slabDeltaInverse[5]+9*IND_b+c);
		
		for(c=0;c<9;c++)
		diff_b[c]=*(slabDiffusionCoeff[5]+9*IND_b+c);
	}
	if(k!=(NK-2))
	{       
	        IND_t=(NJ-1)*(NI-1)*(k+1)+(NI-1)*j+i;
		
		u_t=*(domain+IND_t);
		vol_t=*(volumes+IND_t);
		
		GradX_t=*(GradX+IND_t);
		GradY_t=*(GradY+IND_t);
		GradZ_t=*(GradZ+IND_t);
		
		for(c=0;c<9;c++)
		diff_t[c]=*(diffusionCoeffdomain+9*IND_t+c);
	}
	else
	{       IND_t=(NI-1)*j+i;
		
		u_t=*(slabVariable[4]+IND_t);
		vol_t=*(slabVolume[4]+IND_t);
		
		GradX_t=*(slabGradX[4]+IND_t);
		GradY_t=*(slabGradY[4]+IND_t);
		GradZ_t=*(slabGradZ[4]+IND_t);
		
		
		for(c=0;c<9;c++)
		D_t[c]=*(slabDeltaInverse[4]+9*IND_t+c);
		
		for(c=0;c<9;c++)
		diff_t[c]=*(slabDiffusionCoeff[4]+9*IND_t+c);
	}
		
	
	
	
	DiagCoeff[0]=vol_p/(vol_p+vol_w);              //diagonal coeficent 
	DiagCoeff[1]=vol_p/(vol_p+vol_e);
	DiagCoeff[2]=vol_p/(vol_p+vol_s);
	DiagCoeff[3]=vol_p/(vol_p+vol_n);
	DiagCoeff[4]=vol_p/(vol_p+vol_b);
	DiagCoeff[5]=vol_p/(vol_p+vol_t);
	
	for(c=0;c<9;c++)
	{
	diff_E[c] = (diff_e[c]*vol_p + diff_p[c]*vol_e)/(vol_p+vol_e);
	diff_W[c] = (diff_w[c]*vol_p + diff_p[c]*vol_w)/(vol_p+vol_w);
	diff_N[c] = (diff_n[c]*vol_p + diff_p[c]*vol_n)/(vol_p+vol_n);
	diff_S[c] = (diff_s[c]*vol_p + diff_p[c]*vol_s)/(vol_p+vol_s);
	diff_B[c] = (diff_b[c]*vol_p + diff_p[c]*vol_b)/(vol_p+vol_b);
	diff_T[c] = (diff_t[c]*vol_p + diff_p[c]*vol_t)/(vol_p+vol_t);
	}
	
	// printf("\n %lf %lf %lf   %lf %lf %lf   %lf %lf %lf ",diff_E[0],diff_E[1],diff_E[2],diff_E[3],diff_E[4],diff_E[5],diff_E[6],diff_E[7],diff_E[8]);
	
	
	
	     	
	GradX_E = (GradX_e*vol_p+GradX_p*vol_e)/(vol_p+vol_e);
	GradX_W = (GradX_w*vol_p+GradX_p*vol_w)/(vol_p+vol_w);
	GradX_N = (GradX_n*vol_p+GradX_p*vol_n)/(vol_p+vol_n);
	GradX_S = (GradX_s*vol_p+GradX_p*vol_s)/(vol_p+vol_s);
	GradX_T=  (GradX_t*vol_p+GradX_p*vol_t)/(vol_p+vol_t);
	GradX_B=  (GradX_b*vol_p+GradX_p*vol_b)/(vol_p+vol_b);
		
	GradY_E = (GradY_e*vol_p+GradY_p*vol_e)/(vol_p+vol_e);
	GradY_W = (GradY_w*vol_p+GradY_p*vol_w)/(vol_p+vol_w);
	GradY_N = (GradY_n*vol_p+GradY_p*vol_n)/(vol_p+vol_n);
	GradY_S = (GradY_s*vol_p+GradY_p*vol_s)/(vol_p+vol_s);
	GradY_T=  (GradY_t*vol_p+GradY_p*vol_t)/(vol_p+vol_t);
	GradY_B=  (GradY_b*vol_p+GradY_p*vol_b)/(vol_p+vol_b);
	
	GradZ_E = (GradZ_e*vol_p+GradZ_p*vol_e)/(vol_p+vol_e);
	GradZ_W = (GradZ_w*vol_p+GradZ_p*vol_w)/(vol_p+vol_w);
	GradZ_N = (GradZ_n*vol_p+GradZ_p*vol_n)/(vol_p+vol_n);
	GradZ_S = (GradZ_s*vol_p+GradZ_p*vol_s)/(vol_p+vol_s);
	GradZ_T=  (GradZ_t*vol_p+GradZ_p*vol_t)/(vol_p+vol_t);
	GradZ_B=  (GradZ_b*vol_p+GradZ_p*vol_b)/(vol_p+vol_b);
	
	        
		for(iFace=0;iFace<6;iFace++)
		for(comp=0;comp<3;comp++)
		Area[iFace][comp]=*(SurfArea+IND*6*3+3*iFace+comp);
		
		
		
		Coeff1[0]= diff_W[0]*Area[0][0]+diff_W[3]*Area[0][1]+diff_W[6]* Area[0][2];//west face 
		Coeff1[1]= diff_E[0]*Area[1][0]+diff_E[3]*Area[1][1]+diff_E[6]* Area[1][2];//east face 
		Coeff1[2]= diff_S[0]*Area[2][0]+diff_S[3]*Area[2][1]+diff_S[6]* Area[2][2];//south face
		Coeff1[3]= diff_N[0]*Area[3][0]+diff_N[3]*Area[3][1]+diff_N[6]* Area[3][2];//north face
		Coeff1[4]= diff_B[0]*Area[4][0]+diff_B[3]*Area[4][1]+diff_B[6]* Area[4][2];//bottom face
		Coeff1[5]= diff_T[0]*Area[5][0]+diff_T[3]*Area[5][1]+diff_T[6]* Area[5][2]; //top face
		
		Coeff2[0]= diff_W[1]*Area[0][0]+diff_W[4]*Area[0][1]+diff_W[7]* Area[0][2];//west face 
		Coeff2[1]= diff_E[1]*Area[1][0]+diff_E[4]*Area[1][1]+diff_E[7]* Area[1][2];//east face 
		Coeff2[2]= diff_S[1]*Area[2][0]+diff_S[4]*Area[2][1]+diff_S[7]* Area[2][2];//south face
		Coeff2[3]= diff_N[1]*Area[3][0]+diff_N[4]*Area[3][1]+diff_N[7]* Area[3][2];//north face
		Coeff2[4]= diff_B[1]*Area[4][0]+diff_B[4]*Area[4][1]+diff_B[7]* Area[4][2]; //bottom face
		Coeff2[5]= diff_T[1]*Area[5][0]+diff_T[4]*Area[5][1]+diff_T[7]* Area[5][2]; //top face
		
			
		Coeff3[0]= diff_W[2]*Area[0][0]+diff_W[5]*Area[0][1]+diff_W[8]* Area[0][2];//west face 
		Coeff3[1]= diff_E[2]*Area[1][0]+diff_E[5]*Area[1][1]+diff_E[8]* Area[1][2];//east face 
		Coeff3[2]= diff_S[2]*Area[2][0]+diff_S[5]*Area[2][1]+diff_S[8]* Area[2][2];//south face
		Coeff3[3]= diff_N[2]*Area[3][0]+diff_N[5]*Area[3][1]+diff_N[8]* Area[3][2];//north face
		Coeff3[4]= diff_B[2]*Area[4][0]+diff_B[5]*Area[4][1]+diff_B[8]* Area[4][2];//bottom face
		Coeff3[5]= diff_T[2]*Area[5][0]+diff_T[5]*Area[5][1]+diff_T[8]* Area[5][2];//top face
		
	 	 
	 
	 
	 DiffFlux[0]=Coeff1[0]*GradX_W+Coeff2[0]*GradY_W+Coeff3[0]*GradZ_W;//west face 
			 
     	 DiffFlux[1]=Coeff1[1]*GradX_E+Coeff2[1]*GradY_E+Coeff3[1]*GradZ_E;//east face 
      
     	 DiffFlux[2]=Coeff1[2]*GradX_S+Coeff2[2]*GradY_S+Coeff3[2]*GradZ_S;//south face
	
     	 DiffFlux[3]=Coeff1[3]*GradX_N+Coeff2[3]*GradY_N+Coeff3[3]*GradZ_N;//north face
	 
	 DiffFlux[4]=Coeff1[4]*GradX_B+Coeff2[4]*GradY_B+Coeff3[4]*GradZ_B;// bottom face
	 
	 DiffFlux[5]=Coeff1[5]*GradX_T+Coeff2[5]*GradY_T+Coeff3[5]*GradZ_T;// top face
	
	 
	 
	diagonalW=  DiagCoeff[0]*(D_w[0]*Coeff1[0]+D_w[3]*Coeff2[0]+D_w[6]*Coeff3[0]);
	diagonalE= -DiagCoeff[1]*(D_e[0]*Coeff1[1]+D_e[3]*Coeff2[1]+D_e[6]*Coeff3[1]);
	
	diagonalS=  DiagCoeff[2]*(D_s[1]*Coeff1[2]+D_s[4]*Coeff2[2]+D_s[7]*Coeff3[2]);
	diagonalN= -DiagCoeff[3]*(D_n[1]*Coeff1[3]+D_n[4]*Coeff2[3]+D_n[7]*Coeff3[3]);
	
	diagonalB=  DiagCoeff[4]*(D_b[2]*Coeff1[3]+D_b[5]*Coeff2[3]+D_b[8]*Coeff3[4]);
	diagonalT= -DiagCoeff[5]*(D_t[2]*Coeff1[3]+D_t[5]*Coeff2[3]+D_b[8]*Coeff3[5]); 
	
	
	diagonal=diagonalW+diagonalE+diagonalS+diagonalN+diagonalT+diagonalB;
	  
	

	return (diagonal);
}














