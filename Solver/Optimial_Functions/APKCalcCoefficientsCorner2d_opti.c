#include "APKCalcCoefficientsCorner2d_opti.h"

/**
The function is used to calculate the coefficients for the 2D pressure poisson equation.  Refer function CalcCoefficients_opti for complete documentation.
This function is replica of that function, apart from the fact that it is used only incase of 2d problems. */


double CalcCoefficientsCorner2d_opti(int index[3],double vol_p,double u[18],int *size,double *beta,double *DiaN[18],double *D[36],double *DN[24],double *Diag[6])
{
	int NK,NJ,NI,k,j,i,comp,indexP,iFace;
//	double vol_e,vol_w,vol_s,vol_n,vol_t,vol_b,vol_ws,vol_wn,vol_wt,vol_wb,vol_es,vol_en,vol_et,vol_eb,vol_ts,vol_tn,vol_bs,vol_bn;
	double betaLocal[6][3],diag;
	double C1,C2,C3;
	
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
	
	diag=(betaLocal[0][0]-betaLocal[1][0]-betaLocal[3][1]+betaLocal[2][1]+(betaLocal[0][1]+betaLocal[1][1])*((u[3]/(u[3]+vol_p))-(u[2]/(u[2]+vol_p)))+(betaLocal[0][2]+betaLocal[1][2]+betaLocal[3][2]+betaLocal[2][2])*((u[5]/(u[5]+vol_p))-(u[4]/(u[4]+vol_p)))+(betaLocal[3][0]+betaLocal[2][0])*((u[1]/(u[1]+vol_p))-(u[0]/(u[0]+vol_p))));
	
	C1 = ((u[3]/(u[3]+vol_p))-(u[2]/(u[2]+vol_p)));
	C2 = ((u[5]/(u[5]+vol_p))-(u[4]/(u[4]+vol_p)));
	C3 = ((u[1]/(u[1]+vol_p))-(u[0]/(u[0]+vol_p)));
	
	diag=(betaLocal[0][0]-betaLocal[1][0]-betaLocal[3][1]+betaLocal[2][1]
	+(betaLocal[0][1]+betaLocal[1][1])*((u[3]/(u[3]+vol_p))-(u[2]/(u[2]+vol_p)))
	+(betaLocal[0][2]+betaLocal[1][2]+betaLocal[3][2]+betaLocal[2][2])*((u[5]/(u[5]+vol_p))-(u[4]/(u[4]+vol_p)))
	+(betaLocal[3][0]+betaLocal[2][0])*((u[1]/(u[1]+vol_p))-(u[0]/(u[0]+vol_p))));
		 
	*(Diag[0]+indexP) =  betaLocal[0][0] + C1 * betaLocal[0][1] + C2 *betaLocal[0][2];
	*(Diag[1]+indexP) = -betaLocal[1][0] + C1 * betaLocal[1][1] + C2 *betaLocal[1][2];
	*(Diag[2]+indexP) =  betaLocal[2][1] + C3 * betaLocal[2][0] + C2 *betaLocal[2][2];
	*(Diag[3]+indexP) = -betaLocal[3][1] + C3 * betaLocal[3][0] + C2 *betaLocal[3][2];
	*(Diag[4]+indexP) = 0.0;
	*(Diag[5]+indexP) = 0.0;
	
	*(DiaN[0]+indexP)=-1.0*(betaLocal[0][0]+(betaLocal[3][0]+betaLocal[2][0])*(vol_p/(u[0]+vol_p)));
	
	*(D[0]+indexP) = -1.0*betaLocal[0][0];
	*(D[1]+indexP) = 0.0;
	*(D[2]+indexP) = -1.0*betaLocal[2][0]*(vol_p/(u[0]+vol_p));
	*(D[3]+indexP) = -1.0*betaLocal[3][0]*(vol_p/(u[0]+vol_p));
	
		     
	*(DiaN[1]+indexP)=(betaLocal[1][0]+(betaLocal[3][0]+betaLocal[2][0])*(vol_p/(u[1]+vol_p)));
	
	*(D[6]+indexP) = 0.0;
	*(D[7]+indexP) = betaLocal[1][0];
	*(D[8]+indexP) = betaLocal[2][0]*(vol_p/(u[1]+vol_p));
	*(D[9]+indexP) = betaLocal[3][0]*(vol_p/(u[1]+vol_p));
	
		     
	*(DiaN[2]+indexP)=-1.0*(betaLocal[2][1]+(betaLocal[0][1]+betaLocal[1][1])*(vol_p/(u[2]+vol_p)));
	

	*(D[12]+indexP) = -1.0*betaLocal[0][1]*(vol_p/(u[2]+vol_p));
	*(D[13]+indexP) = -1.0*betaLocal[1][1]*(vol_p/(u[2]+vol_p));
	*(D[14]+indexP) = -1.0*betaLocal[2][1];
	*(D[15]+indexP) = 0.0;
// 	if(i==73 && j==23)
// 	{
// 		printf("SUM = %e\tactual=%e\n,",*(D[12]+indexP)+*(D[13]+indexP)+*(D[14]+indexP),*(DiaN[2]+indexP));
// 		exit(12);
// 	}
	 
	*(DiaN[3]+indexP)=(betaLocal[3][1]+(betaLocal[0][1]+betaLocal[1][1])*(vol_p/(u[3]+vol_p)));
	
	
	*(D[18]+indexP) = betaLocal[0][1]*(vol_p/(u[3]+vol_p));
	*(D[19]+indexP) = betaLocal[1][1]*(vol_p/(u[3]+vol_p));
	*(D[20]+indexP) = 0.0;
	*(D[21]+indexP) = betaLocal[3][1];

// 	*(DiaN[4]+indexP)=-1.0*(betaLocal[0][2]+betaLocal[1][2]+betaLocal[3][2]+betaLocal[2][2])*(vol_p/(u[4]+vol_p));//-1.0*diffCoeff*(betaLocal[4][2]);
// 	
// 	*(DiaN[5]+indexP)=(betaLocal[0][2]+betaLocal[1][2]+betaLocal[3][2]+betaLocal[2][2])*(vol_p/(u[5]+vol_p));//-1.0*diffCoeff*(-betaLocal[5][2]); 
	
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
	
	return (diag);
}


