#include "APKFind_Mushy_Source.h" 

/**
This function calculates the extra source term for enthalpy equation incase of mushy region problems.  This function takes the diffusion coefficients of pressure poission equation and calculates the requried term.  For details about the term refer M.Tech thesis on Numerical Simulation of Solidification and Melting Problems.
*/

double Find_Mushy_Source(Base *base,Zone *zone,int INDEX,int index[3])
{

	int NI,NJ,NK,k,j,i,*ptrSize,g,cntDir;
	double *slabVariable[6];
	double Vg,V[6],*CoeffN[6],*diffdia;
	double *variable,value;
	
	ptrSize= zone->Size;
	g = 0;

	diffdia = (zone->grid+g)->Conv_Coeff_P;
	for(i=0;i<6;i++)
		CoeffN[i] = (zone->grid+g)->Conv_Coeff_nb[i];

	
	for(cntDir=0;cntDir<6;cntDir++)
	slabVariable[cntDir] = zone->Flowsol->Slab_Latent_Heat_Content[cntDir];
	
	variable = zone->Flowsol->Latent_Heat_Content;
	
	
	NI=*(ptrSize+0);
	NJ=*(ptrSize+1);
	NK=*(ptrSize+2);
	
	k=index[0];
	j=index[1];
	i=index[2];
	
	if(i!=0)
		V[0]=*(variable+INDEX-1);
	else
		V[0] =*(slabVariable[1]+j);

	if(i!=(NI-2))
		V[1]=*(variable+INDEX+1);
	else
		V[1]= *(slabVariable[0]+j);

	if(j!=0)
		V[2]=*(variable+INDEX-(NI-1));
	else
		V[2]=*(slabVariable[3]+i);

	if(j!=(NJ-2))
		V[3]=*(variable+INDEX+NI-1);
	else
		V[3]=*(slabVariable[2]+i);

	if(k!=0)
		V[4]=*(variable+INDEX-(NI-1)*(NJ-1));
	else
		V[4]=*(slabVariable[5]+(NI-1)*j+i);
		
	if(k!=(NK-2))
		V[5]=*(variable+INDEX+(NI-1)*(NJ-1));
	else
		V[5]=*(slabVariable[4]+(NI-1)*j+i);
		
	Vg = *(variable+INDEX);
	
	value = *(CoeffN[0]+INDEX) *V[0] + *(CoeffN[1]+INDEX) *V[1] + *(CoeffN[2]+INDEX) *V[2] +*(CoeffN[3]+INDEX) *V[3] + *(CoeffN[4]+INDEX) *V[4] + *(CoeffN[5]+INDEX) *V[5]-*(diffdia+INDEX)*Vg;
	return (-value);
}



