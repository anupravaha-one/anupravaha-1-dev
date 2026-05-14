#include "APKCalcCoefficientsCorner.h"

/**
The function is used to calculate the coefficients for the 3D pressure poisson equation. This is
called for only for corner, edge and slab cells.
*/

double CalcCoefficientsCorner(int index[3],double vol_p,double u[18],int *size,double *beta,double *DiaN[18])
{
	int NK,NJ,NI,k,j,i,comp,indexP,iFace;
//	double vol_e,vol_w,vol_s,vol_n,vol_t,vol_b,vol_ws,vol_wn,vol_wt,vol_wb,vol_es,vol_en,vol_et,vol_eb,vol_ts,vol_tn,vol_bs,vol_bn;
	double betaLocal[6][3],diag;
//	double coeffUp[12],coeffWp[4],coeffEp[4],coeffNp[4],coeffSp[4],coeffTp[4],coeffBp[4];
//	double coeffPWS,coeffPWN,coeffPWT,coeffPWB,coeffPEN,coeffPES,coeffPET,coeffPEB,coeffPTN,coeffPBN,coeffPTS,coeffPBS;	
	
	NI=*(size+0);
	NJ=*(size+1);
	NK=*(size+2);
	
	k=index[0],j=index[1],i=index[2];
	
/*	vol_w = u[0];
	vol_e = u[1];
	vol_s = u[2];
	vol_n = u[3];
	vol_b = u[4];
	vol_t = u[5];
	vol_wn = u[6];
	vol_ws = u[7];
	vol_wt = u[8];
	vol_wb = u[9];
	vol_en = u[10];
	vol_es = u[11];
	vol_et = u[12];
	vol_eb = u[13];
	vol_tn = u[14];
	vol_ts = u[15];
	vol_bn = u[16];
	vol_bs = u[17];*/
	
	for(iFace=0;iFace<6;iFace++)
		for(comp=0;comp<3;comp++)
			betaLocal[iFace][comp]=*(beta+((NI-1)*(NJ-1)*k+(NI-1)*j+i)*6*3+3*iFace+comp);
			
	indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
	
	diag=(betaLocal[0][0]-betaLocal[1][0]+betaLocal[2][1]-betaLocal[3][1]+betaLocal[4][2]-betaLocal[5][2]
			+(betaLocal[0][2]+betaLocal[1][2]+betaLocal[3][2]+betaLocal[2][2])*((u[5]/(u[5]+vol_p))-(u[4]/(u[4]+vol_p)))
			+(betaLocal[0][1]+betaLocal[1][1]+betaLocal[4][1]+betaLocal[5][1])*((u[3]/(u[3]+vol_p))-(u[2]/(u[2]+vol_p)))
			+(betaLocal[3][0]+betaLocal[2][0]+betaLocal[4][0]+betaLocal[5][0])*((u[1]/(u[1]+vol_p))-(u[0]/(u[0]+vol_p))));
		 
	*(DiaN[0]+indexP)=-1.0*(betaLocal[0][0]+(betaLocal[3][0]+betaLocal[2][0]+betaLocal[4][0]+betaLocal[5][0])*(vol_p/(u[0]+vol_p)));
		     
	*(DiaN[1]+indexP)=(betaLocal[1][0]+(betaLocal[3][0]+betaLocal[2][0]+betaLocal[4][0]+betaLocal[5][0])*(vol_p/(u[1]+vol_p)));
		     
	*(DiaN[2]+indexP)=-1.0*(betaLocal[2][1]+(betaLocal[0][1]+betaLocal[1][1]+betaLocal[4][1]+betaLocal[5][1])*(vol_p/(u[2]+vol_p)));
		     
	*(DiaN[3]+indexP)=(betaLocal[3][1]+(betaLocal[0][1]+betaLocal[1][1]+betaLocal[4][1]+betaLocal[5][1])*(vol_p/(u[3]+vol_p)));
		      
	*(DiaN[4]+indexP)=-1.0*(betaLocal[4][2]+(betaLocal[0][2]+betaLocal[1][2]+betaLocal[3][2]+betaLocal[2][2])*(vol_p/(u[4]+vol_p)));
	
	*(DiaN[5]+indexP)=(betaLocal[5][2]+(betaLocal[0][2]+betaLocal[1][2]+betaLocal[3][2]+betaLocal[2][2])*(vol_p/(u[5]+vol_p)));
	
	*(DiaN[6]+indexP)=0.0;
	
	*(DiaN[7]+indexP)=0.0;
	
	*(DiaN[8]+indexP)=0.0;
	
	*(DiaN[9]+indexP)=0.0;
	
	*(DiaN[10]+indexP)=0.0;
	
	*(DiaN[11]+indexP)=0.0;
	
	*(DiaN[12]+indexP)=0.0;
	
	*(DiaN[13]+indexP)=0.0;
	
	*(DiaN[14]+indexP)=0.0;
	
	*(DiaN[15]+indexP)=0.0;
	
	*(DiaN[16]+indexP)=0.0;
	
	*(DiaN[17]+indexP)=0.0;
	
//	printf("%lf %lf %lf %lf %lf %lf\n",*(DiaN[0]+indexP),*(DiaN[1]+indexP),*(DiaN[2]+indexP),*(DiaN[3]+indexP),*(DiaN[4]+indexP),*(DiaN[5]+indexP));
	
	return (diag);
}
