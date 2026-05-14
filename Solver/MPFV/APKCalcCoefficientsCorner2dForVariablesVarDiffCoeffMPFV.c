#include "APKCalcCoefficientsCorner2dForVariablesVarDiffCoeffMPFV.h"

/**
This function calculates coeffients of diffusion term in the convection-diffusion equation for overallMomentum. This is used for the case of varying diffusion coefficients. It is used for 2D case and it is called only
for corner, edge, slab cells.
*/

double CalcCoefficientsCorner2dForVariablesVarDiffCoeffMPFV(Zone *zone,int index[3],double vol_p,double u[18],int *size,double *beta,double *DiaN[18],int vartag)
{
	int NK,NJ,NI,k,j,i,comp,indexP,iFace,flagrho;
//	double vol_e,vol_w,vol_s,vol_n,vol_t,vol_b,vol_ws,vol_wn,vol_wt,vol_wb,vol_es,vol_en,vol_et,vol_eb,vol_ts,vol_tn,vol_bs,vol_bn;
	double r_p,betaLocal[6][3],diag,*diffusionCoeffdomain,*slabDiffusionCoeff[6],*Domain;
	double diffCoeff_p,diffCoeff_e,diffCoeff_w,diffCoeff_s,diffCoeff_n,diffCoeff_t,diffCoeff_b,*slabRho[6];//These are Cell center values
	double diffCoeff_Eface,diffCoeff_Wface,diffCoeff_Sface,diffCoeff_Nface,diffCoeff_Tface,diffCoeff_Bface;//These are face center values
	double r_e,r_w,r_s,r_n,r_t,r_b,r_ws,r_wn,r_wt,r_wb,r_es,r_en,r_et,r_eb,r_bs,r_bn,r_ts,r_tn;
	//Zone *zone;
// 	FILE *val41;
// 	val41=fopen("Settings/valu1.dat","a");
	NI=*(size+0);
	NJ=*(size+1);
	NK=*(size+2);
	
	flagrho = zone->Flowsol->flagRho;
	diffusionCoeffdomain = MakeDiffusionCoeffDomain(zone,vartag);
	MakeDiffusionCoeffSlabPtrArray(zone,slabDiffusionCoeff,vartag);
	MakeSlabPtrArray(zone,slabRho,710);
	Domain = MakeDomainArray(zone,710);
	
	k=index[0],j=index[1],i=index[2];
	indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
	
	
	diffCoeff_p=*(diffusionCoeffdomain+(NJ-1)*(NI-1)*k+(NI-1)*j+i);
	//printf("\nvartag=%d\t diffcoeff=%e",vartag,diffCoeff_p);
	r_p=1/(*(Domain+(NJ-1)*(NI-1)*k+(NI-1)*j+i));
	
	if(i!=0)
	{
		diffCoeff_w=*(diffusionCoeffdomain+(NJ-1)*(NI-1)*k+(NI-1)*j+(i-1));
		r_w=1/(*(Domain+(NJ-1)*(NI-1)*k+(NI-1)*j+(i-1)));
// 		printf("westcoeff=%e\n",diffCoeff_w);
	}
	else
	{
		diffCoeff_w =*(slabDiffusionCoeff[1]+(NJ-1)*k+j);
		r_w =1/(*(slabRho[1]+(NJ-1)*k+j));
		//printf("westcoeff=%e\n",diffCoeff_w);
	}
	if(i!=(NI-2))
	{
		diffCoeff_e=*(diffusionCoeffdomain+(NJ-1)*(NI-1)*k+(NI-1)*j+(i+1));
		r_e=1/(*(Domain+(NJ-1)*(NI-1)*k+(NI-1)*j+(i+1)));
		//printf("\n1diffcoeffE=%e\ti=%d",diffCoeff_e,i);
	}
	else
	{
		diffCoeff_e=*(slabDiffusionCoeff[0]+(NJ-1)*k+j);
		r_e=1/(*(slabRho[0]+(NJ-1)*k+j));
		//printf("\n2diffcoeffE=%e\ti=%d",diffCoeff_e,i);
;
	}
	if(j!=0)
	{
		diffCoeff_s=*(diffusionCoeffdomain+(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i);
		r_s=1/(*(Domain+(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i));
	}
	else
	{
		diffCoeff_s=*(slabDiffusionCoeff[3]+(NI-1)*k+i);
		r_s=1/(*(slabRho[3]+(NI-1)*k+i));
	}
	if(j!=(NJ-2))
	{
		diffCoeff_n=*(diffusionCoeffdomain+(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i);
		r_n=1/(*(Domain+(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+i));
	}
	else
	{
		diffCoeff_n=*(slabDiffusionCoeff[2]+(NI-1)*k+i);
		r_n=1/(*(slabRho[2]+(NI-1)*k+i));
	}
	if(k!=0)
	{
		diffCoeff_b=*(diffusionCoeffdomain+(NJ-1)*(NI-1)*(k-1)+(NI-1)*j+i);
		r_b=1/(*(Domain+(NJ-1)*(NI-1)*(k-1)+(NI-1)*j+i));
	}
	else
	{
		diffCoeff_b=*(slabDiffusionCoeff[5]+(NI-1)*j+i);
		r_b=1/(*(slabRho[5]+(NI-1)*j+i));
	}
	if(k!=(NK-2))
	{
		diffCoeff_t=*(diffusionCoeffdomain+(NJ-1)*(NI-1)*(k+1)+(NI-1)*j+i);
		
		r_t=1/(*(Domain+(NJ-1)*(NI-1)*(k+1)+(NI-1)*j+i));
	}
	else
	{
		diffCoeff_t=*(slabDiffusionCoeff[4]+(NI-1)*j+i);
		r_t=1/(*(slabRho[4]+(NI-1)*j+i));
	}
	
	
	
	if(!(i==0 && j==0))
	{
		if(i==0)
		{
			
			r_ws=1/(*(slabRho[1]+(NJ-1)*k+(j-1)));
		}
		else
		{
			if(j==0)
			{
				
				r_ws=1/(*(slabRho[3]+(NI-1)*k+(i-1)));
			}
			else
			{
				r_ws=1/(*(Domain+(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+(i-1)));
				
			}
		}
	}
		
	if(!(i==0 && j==(NJ-2)))
	{
		if(i==0)
		{
			
			r_wn=1/(*(slabRho[1]+(NJ-1)*k+(j+1)));
		}
		else
		{
			if(j==(NJ-2))
			{
				
				r_wn=1/(*(slabRho[2]+(NI-1)*k+(i-1)));
			}
			else
			{
				r_wn=1/(*(Domain+(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+(i-1)));
				
			}
		}
	}
	if(!(i==0 && k==0))
	{
		if(i==0)
		{
			
			r_wb=1/(*(slabRho[1]+(NJ-1)*(k-1)+j));
		}
		else
		{
			if(k==0)
			{
				
				r_wb=1/(*(slabRho[5]+(NI-1)*j+(i-1)));
			}
			else
			{
				r_wb=1/(*(Domain+(NJ-1)*(NI-1)*(k-1)+(NI-1)*j+(i-1)));
				
			}
		}
		
	}
	
	if(!(i==0 && k==(NK-2)))
	{
		if(i==0)
		{
			
			r_wt=1/(*(slabRho[1]+(NJ-1)*(k+1)+j));
		}
		else
		{
			if(k==(NK-2))
			{
				
				r_wt=1/(*(slabRho[4]+(NI-1)*j+(i-1)));
			}
			else
			{
				r_wt=1/(*(Domain+(NJ-1)*(NI-1)*(k+1)+(NI-1)*j+(i-1)));
				
			}
		}
		
	}
	if(!(i==(NI-2) && j==0))
	{
		if(i==(NI-2))
		{
			
			r_es=1/(*(slabRho[0]+(NJ-1)*k+(j-1)));
		}
		else
		{
			if(j==0)
			{
				
				r_es=1/(*(slabRho[3]+(NI-1)*k+(i+1)));
			}
			else
			{
				r_es=1/(*(Domain+(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+(i+1)));
				
			}
		}
	}
		
	if(!(i==(NI-2) && j==(NJ-2)))
	{
		if(i==(NI-2))
		{
			
			r_en=1/(*(slabRho[0]+(NJ-1)*k+(j+1)));
		}
		else
		{
			if(j==(NJ-2))
			{
				
				r_en=1/(*(slabRho[2]+(NI-1)*k+(i+1)));
			}
			else
			{
				r_en=1/(*(Domain+(NJ-1)*(NI-1)*k+(NI-1)*(j+1)+(i+1)));
				
			}
		}
	}
	
	if(!(i==(NI-2) && k==0))
	{
		if(i==(NI-2))
		{
			
			r_eb=1/(*(slabRho[0]+(NJ-1)*(k-1)+j));
		}
		else
		{
			if(k==0)
			{
				
				r_eb=1/(*(slabRho[5]+(NI-1)*j+(i+1)));
			}
			else
			{
				r_eb=1/(*(Domain+(NJ-1)*(NI-1)*(k-1)+(NI-1)*j+(i+1)));
				
			}
		}
	}
	
	if(!(i==(NI-2) && k==(NK-2)))
	{
		if(i==(NI-2))
		{
			
			r_et= 1/(*(slabRho[0]+(NJ-1)*(k+1)+j));
		}
		else
		{
			if(k==(NK-2))
			{
				
				r_et=1/(*(slabRho[4]+(NI-1)*j+(i+1)));
			}
			else
			{
				r_et=1/(*(Domain+(NJ-1)*(NI-1)*(k+1)+(NI-1)*j+(i+1)));
				
			}
		}
	}	
	if(!(k==0 && j==0))
	{
		if(k==0)
		{
			
			r_bs=1/(*(slabRho[5]+(NI-1)*(j-1)+i));
		}
		else
		{
			if(j==0)
			{
				
				r_bs=1/(*(slabRho[3]+(NI-1)*(k-1)+i));
			}
			else
			{
				r_bs=1/(*(Domain+(NJ-1)*(NI-1)*(k-1)+(NI-1)*(j-1)+i));
				
			}
		}
		
	}
	if(!(k==0 && j==(NJ-2)))
	{
		if(k==0 )
		{
			
			r_bn=1/(*(slabRho[5]+(NI-1)*(j+1)+i));
		}
		else
		{
			if(j==(NJ-2))
			{
				
				r_bn=1/(*(slabRho[2]+(NI-1)*(k-1)+i));
			}
			else
			{
				r_bn=1/(*(Domain+(NJ-1)*(NI-1)*(k-1)+(NI-1)*(j+1)+i));
				
			}
		}
		
	}
	if(!(k==(NK-2) && j==0))
	{
		if(k==(NK-2))
		{
			
			r_ts=1/(*(slabRho[4]+(NI-1)*(j-1)+i));
		}
		else
		{
			if(j==0)
			{
				
				r_ts=1/(*(slabRho[3]+(NI-1)*(k+1)+i));
			}
			else
			{
				r_ts=1/(*(Domain+(NJ-1)*(NI-1)*(k+1)+(NI-1)*(j-1)+i));
				
			}
		}
	}
	
	if(!(k==(NK-2) && j==(NJ-2)))
	{
		if(k==(NK-2))
		{
			
			r_tn=1/(*(slabRho[4]+(NI-1)*(j+1)+i));
		}
		else
		{
			if(j==(NJ-2))
			{
				
				r_tn=1/(*(slabRho[2]+(NI-1)*(k+1)+i));
			}
			else
			{
				r_tn=1/(*(Domain+(NJ-1)*(NI-1)*(k+1)+(NI-1)*(j+1)+i));
				
			}
		}
	}
	
	
	//end of new added 
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
	
	
	diffCoeff_Eface = (diffCoeff_e*vol_p + diffCoeff_p*u[1])/(vol_p+u[1]);
	diffCoeff_Wface = (diffCoeff_w*vol_p + diffCoeff_p*u[0])/(vol_p+u[0]);
	diffCoeff_Nface = (diffCoeff_n*vol_p + diffCoeff_p*u[3])/(vol_p+u[3]);
	diffCoeff_Sface = (diffCoeff_s*vol_p + diffCoeff_p*u[2])/(vol_p+u[2]);
	diffCoeff_Tface = (diffCoeff_t*vol_p + diffCoeff_p*u[5])/(vol_p+u[5]);
	diffCoeff_Bface = (diffCoeff_b*vol_p + diffCoeff_p*u[4])/(vol_p+u[4]);
// 	if(indexP==0)
//  	printf("\nDiffusion coefficient = %e\t%e\t%e\t%e\t%e\t%e\n",diffCoeff_Eface,diffCoeff_w,diffCoeff_Nface,diffCoeff_Sface,diffCoeff_Tface,diffCoeff_Bface);
	for(iFace=0;iFace<6;iFace++)
		for(comp=0;comp<3;comp++)
			betaLocal[iFace][comp]=*(beta+((NI-1)*(NJ-1)*k+(NI-1)*j+i)*6*3+3*iFace+comp);
			
	indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
	
	
	diag=(diffCoeff_Wface*betaLocal[0][0]-diffCoeff_Eface*betaLocal[1][0]-diffCoeff_Nface*betaLocal[3][1]+diffCoeff_Sface*betaLocal[2][1]+(diffCoeff_Wface*betaLocal[0][1]+diffCoeff_Eface*betaLocal[1][1])*((u[3]/(u[3]+vol_p))-(u[2]/(u[2]+vol_p)))+(diffCoeff_Wface*betaLocal[0][2]+diffCoeff_Eface*betaLocal[1][2]+diffCoeff_Nface*betaLocal[3][2]+diffCoeff_Sface*betaLocal[2][2])*((u[5]/(u[5]+vol_p))-(u[4]/(u[4]+vol_p)))+(diffCoeff_Nface*betaLocal[3][0]+diffCoeff_Sface*betaLocal[2][0])*((u[1]/(u[1]+vol_p))-(u[0]/(u[0]+vol_p))))*r_p;
	//printf("\ndiag=%e",diag);
	
	*(DiaN[0]+indexP)=-1.0*(diffCoeff_Wface*betaLocal[0][0]+(diffCoeff_Nface*betaLocal[3][0]+diffCoeff_Sface*betaLocal[2][0])*(vol_p/(u[0]+vol_p)))*r_w;
		
	*(DiaN[1]+indexP)=(diffCoeff_Eface*betaLocal[1][0]+(diffCoeff_Nface*betaLocal[3][0]+diffCoeff_Sface*betaLocal[2][0])*(vol_p/(u[1]+vol_p)))*r_e;
		
	*(DiaN[2]+indexP)=-(1.0*(diffCoeff_Sface*betaLocal[2][1]+(diffCoeff_Wface*betaLocal[0][1]+diffCoeff_Eface*betaLocal[1][1])*(vol_p/(u[2]+vol_p))))*r_s;
	
	*(DiaN[3]+indexP)=(diffCoeff_Nface*betaLocal[3][1]+(diffCoeff_Wface*betaLocal[0][1]+diffCoeff_Eface*betaLocal[1][1])*(vol_p/(u[3]+vol_p)))*r_n;
	
	*(DiaN[4]+indexP)=-(1.0*(diffCoeff_Wface*betaLocal[0][2]+diffCoeff_Eface*betaLocal[1][2]+diffCoeff_Nface*betaLocal[3][2]+diffCoeff_Sface*betaLocal[2][2])*(vol_p/(u[4]+vol_p)))*r_b;//-1.0*diffCoeff*(betaLocal[4][2]);
	
	*(DiaN[5]+indexP)=((diffCoeff_Wface*betaLocal[0][2]+diffCoeff_Eface*betaLocal[1][2]+diffCoeff_Nface*betaLocal[3][2]+diffCoeff_Sface*betaLocal[2][2])*(vol_p/(u[5]+vol_p)))*r_t;//-1.0*diffCoeff*(-betaLocal[5][2]); 

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

	
	//printf("\nretuned diag=%e",diag);
	return (diag);
	//printf("\noutretuned diag=%e",*(DiaN[5]+indexP));
	
}


