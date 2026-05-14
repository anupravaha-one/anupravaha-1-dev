#include "APKCalcCoefficientsCorner_opti.h"

/**
The function is used to calculate the coefficients for the 3D pressure poisson equation. This is called for only for corner, edge and slab cells.
Refer function CalcCoefficients_opti for complete documentation.
This function is replica of that function, apart from the fact that it is used only incase of slabs.*/


double CalcCoefficientsCorner_opti(int index[3],double vol_p,double u[18],int *size,double *beta,double *DiaN[18],double *D[36],double *DN[24],double *Diag[6])
{
	int NK,NJ,NI,k,j,i,comp,indexP,iFace;

	double betaLocal[6][3],diag;
	double C1,C2,C3;
	
	NI=*(size+0);
	NJ=*(size+1);
	NK=*(size+2);
	
	k=index[0],j=index[1],i=index[2];
	
	
	for(iFace=0;iFace<6;iFace++)
		for(comp=0;comp<3;comp++)
			betaLocal[iFace][comp]=*(beta+((NI-1)*(NJ-1)*k+(NI-1)*j+i)*6*3+3*iFace+comp);
			
	indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
	
	diag=(betaLocal[0][0]-betaLocal[1][0]+betaLocal[2][1]-betaLocal[3][1]+betaLocal[4][2]-betaLocal[5][2]
			+(betaLocal[0][2]+betaLocal[1][2]+betaLocal[3][2]+betaLocal[2][2])*((u[5]/(u[5]+vol_p))-(u[4]/(u[4]+vol_p)))
			+(betaLocal[0][1]+betaLocal[1][1]+betaLocal[4][1]+betaLocal[5][1])*((u[3]/(u[3]+vol_p))-(u[2]/(u[2]+vol_p)))
			+(betaLocal[3][0]+betaLocal[2][0]+betaLocal[4][0]+betaLocal[5][0])*((u[1]/(u[1]+vol_p))-(u[0]/(u[0]+vol_p))));
	C1=((u[5]/(u[5]+vol_p))-(u[4]/(u[4]+vol_p)));
	C2=((u[3]/(u[3]+vol_p))-(u[2]/(u[2]+vol_p)));
	C3=((u[1]/(u[1]+vol_p))-(u[0]/(u[0]+vol_p)));
	
	diag=(betaLocal[0][0]-betaLocal[1][0]+betaLocal[2][1]-betaLocal[3][1]+betaLocal[4][2]-betaLocal[5][2]
			+(betaLocal[0][2]+betaLocal[1][2]+betaLocal[3][2]+betaLocal[2][2])*((u[5]/(u[5]+vol_p))-(u[4]/(u[4]+vol_p)))
			+(betaLocal[0][1]+betaLocal[1][1]+betaLocal[4][1]+betaLocal[5][1])*((u[3]/(u[3]+vol_p))-(u[2]/(u[2]+vol_p)))
			+(betaLocal[3][0]+betaLocal[2][0]+betaLocal[4][0]+betaLocal[5][0])*((u[1]/(u[1]+vol_p))-(u[0]/(u[0]+vol_p))));
	
	*(Diag[0]+indexP) =  betaLocal[0][0] + C1 * betaLocal[0][2] + C2 *betaLocal[0][1];
	*(Diag[1]+indexP) = -betaLocal[1][0] + C1 * betaLocal[1][2] + C2 *betaLocal[1][1];
	*(Diag[2]+indexP) =  betaLocal[2][1] + C1 * betaLocal[2][2] + C3 *betaLocal[2][0];
	*(Diag[3]+indexP) = -betaLocal[3][1] + C1 * betaLocal[3][2] + C3 *betaLocal[3][0];
	*(Diag[4]+indexP) =  betaLocal[4][2] + C2 * betaLocal[4][1] + C3 *betaLocal[4][0];
	*(Diag[5]+indexP) = -betaLocal[5][2] + C2 * betaLocal[5][1] + C3 *betaLocal[5][0];
	
	*(DiaN[0]+indexP)=-1.0*(betaLocal[0][0]+(betaLocal[3][0]+betaLocal[2][0]+betaLocal[4][0]+betaLocal[5][0])*(vol_p/(u[0]+vol_p)));
	
	*(D[0]+indexP)=-1.0*(betaLocal[0][0]);
	*(D[1]+indexP)=0.0;
	*(D[2]+indexP)=-1.0*betaLocal[2][0]*(vol_p/(u[0]+vol_p));
	*(D[3]+indexP)=-1.0*betaLocal[3][0]*(vol_p/(u[0]+vol_p));
	*(D[4]+indexP)=-1.0*betaLocal[4][0]*(vol_p/(u[0]+vol_p));
	*(D[5]+indexP)=-1.0*betaLocal[5][0]*(vol_p/(u[0]+vol_p));
		     
	*(DiaN[1]+indexP)=(betaLocal[1][0]+(betaLocal[3][0]+betaLocal[2][0]+betaLocal[4][0]+betaLocal[5][0])*(vol_p/(u[1]+vol_p)));
	
	*(D[6]+indexP)=0.0;
	*(D[7]+indexP)=betaLocal[1][0];
	*(D[8]+indexP)=betaLocal[2][0]*(vol_p/(u[1]+vol_p));
	*(D[9]+indexP)=betaLocal[3][0]*(vol_p/(u[1]+vol_p));
	*(D[10]+indexP)=betaLocal[4][0]*(vol_p/(u[1]+vol_p));
	*(D[11]+indexP)=betaLocal[5][0]*(vol_p/(u[1]+vol_p));
		     
	*(DiaN[2]+indexP)=-1.0*(betaLocal[2][1]+(betaLocal[0][1]+betaLocal[1][1]+betaLocal[4][1]+betaLocal[5][1])*(vol_p/(u[2]+vol_p)));
	
	
	*(D[12]+indexP)=-1.0*betaLocal[0][1]*(vol_p/(u[2]+vol_p));
	*(D[13]+indexP)=-1.0*betaLocal[1][1]*(vol_p/(u[2]+vol_p));
	*(D[14]+indexP)=-1.0*betaLocal[2][1];
	*(D[15]+indexP)=0.0;
	*(D[16]+indexP)=-1.0*betaLocal[4][1]*(vol_p/(u[2]+vol_p));
	*(D[17]+indexP)=-1.0*betaLocal[5][1]*(vol_p/(u[2]+vol_p));
		     
	*(DiaN[3]+indexP)=(betaLocal[3][1]+(betaLocal[0][1]+betaLocal[1][1]+betaLocal[4][1]+betaLocal[5][1])*(vol_p/(u[3]+vol_p)));
	
	
	*(D[18]+indexP)=betaLocal[0][1]*(vol_p/(u[3]+vol_p));
	*(D[19]+indexP)=betaLocal[1][1]*(vol_p/(u[3]+vol_p));
	*(D[20]+indexP)=0.0;
	*(D[21]+indexP)=betaLocal[3][1];
	*(D[22]+indexP)=betaLocal[4][1]*(vol_p/(u[3]+vol_p));
	*(D[23]+indexP)=betaLocal[5][1]*(vol_p/(u[3]+vol_p));
		      
	*(DiaN[4]+indexP)=-1.0*(betaLocal[4][2]+(betaLocal[0][2]+betaLocal[1][2]+betaLocal[3][2]+betaLocal[2][2])*(vol_p/(u[4]+vol_p)));
	
	
	*(D[24]+indexP)=-1.0*betaLocal[0][2]*(vol_p/(u[4]+vol_p));
	*(D[25]+indexP)=-1.0*betaLocal[1][2]*(vol_p/(u[4]+vol_p));
	*(D[26]+indexP)=-1.0*betaLocal[2][2]*(vol_p/(u[4]+vol_p));
	*(D[27]+indexP)=-1.0*betaLocal[3][2]*(vol_p/(u[4]+vol_p));
	*(D[28]+indexP)=-1.0*(betaLocal[4][2]);
	*(D[29]+indexP)=0.0;
	
	
	*(DiaN[5]+indexP)=(betaLocal[5][2]+(betaLocal[0][2]+betaLocal[1][2]+betaLocal[3][2]+betaLocal[2][2])*(vol_p/(u[5]+vol_p)));
	
	
	*(D[30]+indexP)=betaLocal[0][2]*(vol_p/(u[5]+vol_p));
	*(D[31]+indexP)=betaLocal[1][2]*(vol_p/(u[5]+vol_p));
	*(D[32]+indexP)=betaLocal[2][2]*(vol_p/(u[5]+vol_p));
	*(D[33]+indexP)=betaLocal[3][2]*(vol_p/(u[5]+vol_p));
	*(D[34]+indexP)=0.0;
	*(D[35]+indexP)=betaLocal[5][2];

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
