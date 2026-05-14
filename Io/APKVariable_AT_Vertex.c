#include "APKVariable_AT_Vertex.h"
/** This function is for finding out the values at the vertices of the cell.  By default the solver is written for cell centered variables.  This function is written in order make the data handling easier for Tecplot 10.0**/


///3D case needs to be checked.

double Variable_AT_Vertex(Base *base,Zone *zone,double *variable,int vartag)
{
	int NI,NJ,NK,k,j,i,*ptrSize,INDEX,indexP,K_Increment,J_Increment;
	double *variable_vertex,*slabVariable[6],*slabVolume[6],*domainVolume;
	double V[12],Vol[12];
	double V_p,Vol_p;
	
	ptrSize= zone->Size;
	MakeSlabPtrArray(zone,slabVariable,vartag);
	MakeSlabPtrArray(zone,slabVolume,5);
	
	domainVolume=MakeDomainArray(zone,5);
	
	variable_vertex = zone->Flowsol->VAR_VERTEX;
	
	if(base->BaseParameters->flag2d==1)
	{
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		for (j=1;j<NJ-1;j++)
		for (i=1;i<NI-1;i++)
		{
			INDEX=(NI-1)*j+i;
			indexP =NI*j+i;
			
			V_p  = *(variable+INDEX);
			Vol_p= *(domainVolume+INDEX);
			
			
			V[0]=*(variable+INDEX-1);
			Vol[0]=*(domainVolume+INDEX-1);
			
			V[1]=*(variable+INDEX-(NI-1));
			Vol[1] =*(domainVolume+INDEX-(NI-1));
			
			V[2]=*(variable+INDEX-(NI-1)-1);
			Vol[2]=*(domainVolume+INDEX-(NI-1)-1);
			
			*(variable_vertex+indexP) = (V_p * Vol[2] + V[2] *Vol_p + Vol[1] * V[0] + Vol[0] *V[1])/(Vol_p+Vol[2]+Vol[0]+Vol[1]);
			//printf("Inter- polated value =%e\ti=%d\tj=%d\n",*(variable_vertex+indexP),i,j);
		}
		
		i=0;
		for(j=1;j<NJ-1;j++)
		{
			INDEX=(NI-1)*j+i;
			indexP =NI*j+i;
			
			V_p  = *(variable+INDEX);
			Vol_p= *(domainVolume+INDEX);
			
			V[0]=*(slabVariable[1]+j);
			Vol[0]=*(slabVolume[1]+j);
			
			V[1]=*(variable+INDEX-(NI-1));
			Vol[1] =*(domainVolume+INDEX-(NI-1));
			
			V[2]=*(slabVariable[1]+j-1);
			Vol[2]=*(slabVolume[1]+j-1);
			
			*(variable_vertex+indexP) = (V_p * Vol[2] + V[2] *Vol_p + Vol[1] * V[0] + Vol[0] *V[1])/(Vol_p+Vol[2]+Vol[0]+Vol[1]);
			//printf("Inter- polated value =%e\ti=%d\tj=%d\n",*(variable_vertex+indexP),i,j);
		}
		
		i = NI-1;
		for(j=1;j<NJ-1;j++)
		{
			INDEX=(NI-1)*j+i-1;
			indexP =NI*j+i;
			
			V_p  = *(variable+INDEX);
			Vol_p= *(domainVolume+INDEX);
			
			V[0]=*(slabVariable[0]+j);
			Vol[0]=*(slabVolume[0]+j);
			
			V[1]=*(variable+INDEX-(NI-1));
			Vol[1] =*(domainVolume+INDEX-(NI-1));
			
			V[2]=*(slabVariable[0]+j-1);
			Vol[2]=*(slabVolume[0]+j-1);
			
			*(variable_vertex+indexP) = (V_p * Vol[2] + V[2] *Vol_p + Vol[1] * V[0] + Vol[0] *V[1])/(Vol_p+Vol[2]+Vol[0]+Vol[1]);
			//printf("Inter- polated value =%e\ti=%d\tj=%d\n",*(variable_vertex+indexP),i,j);
		}
		
		j=0;
		
		for(i=1;i<NI-1;i++)
		{
			INDEX=(NI-1)*j+i;
			indexP =NI*j+i;
			
			V_p  = *(variable+INDEX);
			Vol_p= *(domainVolume+INDEX);
			
			V[0]=*(slabVariable[3]+i);
			Vol[0]=*(slabVolume[3]+i);
			
			V[1]=*(variable+INDEX-1);
			Vol[1] =*(domainVolume+INDEX-1);
			
			V[2]=*(slabVariable[3]+i-1);
			Vol[2]=*(slabVolume[3]+i-1);
			
			*(variable_vertex+indexP) = (V_p * Vol[2] + V[2] *Vol_p + Vol[1] * V[0] + Vol[0] *V[1])/(Vol_p+Vol[2]+Vol[0]+Vol[1]);
			//printf("Inter- polated value =%e\ti=%d\tj=%d\n",*(variable_vertex+indexP),i,j);
		}
		
		j=NJ-1;
		
		for(i=1;i<NI-1;i++)
		{
			INDEX=(NI-1)*(j-1)+i;
			indexP =NI*j+i;
			
			V_p  = *(variable+INDEX);
			Vol_p= *(domainVolume+INDEX);
			
			V[0]=*(slabVariable[2]+i);
			Vol[0]=*(slabVolume[2]+i);
			
			V[1]=*(variable+INDEX-1);
			Vol[1] =*(domainVolume+INDEX-1);
			
			V[2]=*(slabVariable[2]+i-1);
			Vol[2]=*(slabVolume[2]+i-1);
			
			*(variable_vertex+indexP) = (V_p * Vol[2] + V[2] *Vol_p + Vol[1] * V[0] + Vol[0] *V[1])/(Vol_p+Vol[2]+Vol[0]+Vol[1]);
			//printf("Inter- polated value =%e\ti=%d\tj=%d\n",*(variable_vertex+indexP),i,j);
		}
		
		
		
		//First Corner Point formulation is written here
		i=0;
		j=0;
		
		INDEX=0;
		indexP =0;
		
		V_p  = *(variable+INDEX);
		Vol_p= *(domainVolume+INDEX);
		
		V[0]=*(slabVariable[3]+i);
		Vol[0]=*(slabVolume[3]+i);
			
		V[1]=*(slabVariable[1]+j);
		Vol[1] =*(slabVolume[1]+j);
		
		//printf("V[0]=%e\tV[1]=%e\tV_p=%e\nVol[0]=%e\tVol[1]=%e\tVol_p=%e\n\n",V[0],V[1],V_p,Vol[0],Vol[1],Vol_p);
		
		*(variable_vertex+indexP) = ((Vol_p *V[0] + Vol[0] * V_p)/(Vol_p+Vol[0])) + ((Vol_p *V[1] + Vol[1] * V_p)/(Vol_p+Vol[1])) - V_p;
		//printf("Inter- polated value =%e\ti=%d\tj=%d\n",*(variable_vertex+indexP),i,j);
		
		
		//Second Corner Point formulation is written here
		i=0;
		j=NJ-1;
		INDEX=(NI-1)*(j-1)+i;
		indexP =NI*j+i;
		V_p  = *(variable+INDEX);
		Vol_p= *(domainVolume+INDEX);
		
		V[0]=*(slabVariable[2]+i);
		Vol[0]=*(slabVolume[2]+i);
			
		V[1]=*(slabVariable[1]+j-1);
		Vol[1] =*(slabVolume[1]+j-1);
		//printf("V[0]=%e\tV[1]=%e\tV_p=%e\nVol[0]=%e\tVol[1]=%e\tVol_p=%e\n\n",V[0],V[1],V_p,Vol[0],Vol[1],Vol_p);
		*(variable_vertex+indexP) = ((Vol_p *V[0] + Vol[0] * V_p)/(Vol_p+Vol[0])) + ((Vol_p *V[1] + Vol[1] * V_p)/(Vol_p+Vol[1])) - V_p;
		//printf("Inter- polated value =%e\ti=%d\tj=%d\n",*(variable_vertex+indexP),i,j);
		
		//Third Corner Point formulation is written here
		i=NI-1;
		j=0;
		
		INDEX=(NI-1)*(j)+i-1;
		indexP =NI*j+i;
		
		V_p  = *(variable+INDEX);
		Vol_p= *(domainVolume+INDEX);
		
		V[0]=*(slabVariable[3]+i-1);
		Vol[0]=*(slabVolume[3]+i-1);
			
		V[1]=*(slabVariable[0]+j);
		Vol[1] =*(slabVolume[0]+j);
		
		//printf("V[0]=%e\tV[1]=%e\tV_p=%e\nVol[0]=%e\tVol[1]=%e\tVol_p=%e\n\n",V[0],V[1],V_p,Vol[0],Vol[1],Vol_p);
		
		*(variable_vertex+indexP) = ((Vol_p *V[0] + Vol[0] * V_p)/(Vol_p+Vol[0])) + ((Vol_p *V[1] + Vol[1] * V_p)/(Vol_p+Vol[1])) - V_p;
		//printf("Inter- polated value =%e\ti=%d\tj=%d\n",*(variable_vertex+indexP),i,j);
		
		//Fourth Corner Point formulation is written here
		i=NI-1;
		j=NJ-1;
		
		INDEX=(NI-1)*(j-1)+i-1;
		indexP =NI*j+i;
		
		V_p  = *(variable+INDEX);
		Vol_p= *(domainVolume+INDEX);
		
		V[0]=*(slabVariable[2]+i-1);
		Vol[0]=*(slabVolume[2]+i-1);
			
		V[1]=*(slabVariable[0]+j-1);
		Vol[1] =*(slabVolume[0]+j-1);
		//printf("V[0]=%e\tV[1]=%e\tV_p=%e\nVol[0]=%e\tVol[1]=%e\tVol_p=%e\n\n",V[0],V[1],V_p,Vol[0],Vol[1],Vol_p);
		
		*(variable_vertex+indexP) = ((Vol_p *V[0] + Vol[0] * V_p)/(Vol_p+Vol[0])) + ((Vol_p *V[1] + Vol[1] * V_p)/(Vol_p+Vol[1])) - V_p;
		//printf("Inter- polated value =%e\ti=%d\tj=%d\n",*(variable_vertex+indexP),i,j);
		
		
		indexP=NI*NJ;
		INDEX =0;
		for (j=0;j<NJ;j++)
		for (i=0;i<NI;i++)
		{
			*(variable_vertex+indexP) = *(variable_vertex+INDEX);
			indexP++;
			INDEX++;
		}
	}
	
	else
	{
			
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		//printf("NI NJ NK are %d\t%d\t%d\n",NI,NJ,NK);
		K_Increment = ((NI-1)*(NJ-1));
		J_Increment = NI-1;
		/* Domain looping*/
		for (k=1;k<NK-1;k++)
		for (j=1;j<NJ-1;j++)
		for (i=1;i<NI-1;i++)
		{
			INDEX=(NJ-1)*(NI-1)*k+(NI-1)*j+i;
			indexP =NJ*NI*k+NI*j+i;
			
			V_p  = *(variable+INDEX);
			Vol_p= *(domainVolume+INDEX);

			
			V[0]=*(variable+INDEX-1);
			Vol[0]=*(domainVolume+INDEX-1);

			V[1]=*(variable+INDEX-(NI-1));
			Vol[1] =*(domainVolume+INDEX-(NI-1));

			V[2]=*(variable+INDEX-(NI-1)-1);
			Vol[2]=*(domainVolume+INDEX-(NI-1)-1);

			V[3]=*(variable+INDEX-K_Increment);
			Vol[3]=*(domainVolume+INDEX-K_Increment);
			
			V[4]=*(variable+INDEX-K_Increment-1);
			Vol[4]=*(domainVolume+INDEX-K_Increment-1);

			V[5]=*(variable+INDEX-K_Increment-(NI-1));
			Vol[5]=*(domainVolume+INDEX-K_Increment-(NI-1));
			
			V[6]=*(variable+INDEX-K_Increment-(NI-1)-1);
			Vol[6]=*(domainVolume+INDEX-K_Increment-(NI-1)-1);
			
			*(variable_vertex+indexP) = (V_p * Vol[6] + V[6] *Vol_p + Vol[0] * V[5] + Vol[5] *V[0]+Vol[1] * V[4] + Vol[4] *V[1]+Vol[2] * V[3] + Vol[3] *V[2])/(Vol_p+Vol[2]+Vol[0]+Vol[1]+Vol[3]+Vol[4]+Vol[5]+Vol[6]);
//			printf("1Inter- polated value =%e\ti=%d\tj=%d\tk=%d\n",*(variable_vertex+indexP),i,j,k);
		}
	
		
		i=0;/*west*/
		for(k=1;k<NK-1;k++)
		for(j=1;j<NJ-1;j++)
		{
			INDEX=(NJ-1)*(NI-1)*k+(NI-1)*j+i;
			indexP =NJ*NI*k+NI*j+i;
			
			V_p  = *(variable+INDEX);
			Vol_p= *(domainVolume+INDEX);
			
			V[0]  = *(variable+INDEX-(NI-1));
			Vol[0]= *(domainVolume+INDEX-(NI-1));
			
			V[1]  = *(variable+INDEX-K_Increment);
			Vol[1]= *(domainVolume+INDEX-K_Increment);
			
			V[2]  = *(variable+INDEX-K_Increment-(NI-1));
			Vol[2]= *(domainVolume+INDEX-K_Increment-(NI-1));
			
			
			V[3]=*(slabVariable[1]+(NJ-1)*k+j);
			Vol[3]=*(slabVolume[1]+(NJ-1)*k+j);
			
			V[4]=*(slabVariable[1]+(NJ-1)*k+j-1);
			Vol[4]=*(slabVolume[1]+(NJ-1)*k+j-1);
			
			V[5]=*(slabVariable[1]+(NJ-1)*(k-1)+j);
			Vol[5]=*(slabVolume[1]+(NJ-1)*(k-1)+j);
			
			V[6]=*(slabVariable[1]+(NJ-1)*(k-1)+j-1);
			Vol[6]=*(slabVolume[1]+(NJ-1)*(k-1)+j-1);
			
			*(variable_vertex+indexP) = (V_p * Vol[6] + V[6] *Vol_p + Vol[0] * V[5] + Vol[5] *V[0]+Vol[1] * V[4] + Vol[4] *V[1]+Vol[2] * V[3] + Vol[3] *V[2])/(Vol_p+Vol[2]+Vol[0]+Vol[1]+Vol[3]+Vol[4]+Vol[5]+Vol[6]);
//			printf("2Inter- polated value =%e\ti=%d\tj=%d\n",*(variable_vertex+indexP),i,j);
		}
		
		i=NI-1;/*east*/
		for(k=1;k<NK-1;k++)
		for(j=1;j<NJ-1;j++)
		{
			INDEX=(NJ-1)*(NI-1)*k+(NI-1)*j+i-1;
			indexP =NJ*NI*k+NI*j+i;
			
			V_p  = *(variable+INDEX);
			Vol_p= *(domainVolume+INDEX);
			
			V[0]  = *(variable+INDEX-(NI-1));
			Vol[0]= *(domainVolume+INDEX-(NI-1));
			
			V[1]  = *(variable+INDEX-K_Increment);
			Vol[1]= *(domainVolume+INDEX-K_Increment);
			
			V[2]  = *(variable+INDEX-K_Increment-(NI-1));
			Vol[2]= *(domainVolume+INDEX-K_Increment-(NI-1));
			
			
			V[3]=*(slabVariable[0]+(NJ-1)*k+j);
			Vol[3]=*(slabVolume[0]+(NJ-1)*k+j);
			
			V[4]=*(slabVariable[0]+(NJ-1)*k+j-1);
			Vol[4]=*(slabVolume[0]+(NJ-1)*k+j-1);
			
			V[5]=*(slabVariable[0]+(NJ-1)*(k-1)+j);
			Vol[5]=*(slabVolume[0]+(NJ-1)*(k-1)+j);
			
			V[6]=*(slabVariable[0]+(NJ-1)*(k-1)+j-1);
			Vol[6]=*(slabVolume[0]+(NJ-1)*(k-1)+j-1);
			
			*(variable_vertex+indexP) = (V_p * Vol[6] + V[6] *Vol_p + Vol[0] * V[5] + Vol[5] *V[0]+Vol[1] * V[4] + Vol[4] *V[1]+Vol[2] * V[3] + Vol[3] *V[2])/(Vol_p+Vol[2]+Vol[0]+Vol[1]+Vol[3]+Vol[4]+Vol[5]+Vol[6]);
//			printf("3Inter- polated value =%e\ti=%d\tj=%d\n",*(variable_vertex+indexP),i,j);
		}
		
		j=0;/*south*/
		
		for(k=1;k<NK-1;k++)
		for(i=1;i<NI-1;i++)
		{
			INDEX=(NJ-1)*(NI-1)*k+(NI-1)*j+i;
			indexP =NJ*NI*k+NI*j+i;
			
			V_p  = *(variable+INDEX);
			Vol_p= *(domainVolume+INDEX);
			
			V[0]  = *(variable+INDEX-1);
			Vol[0]= *(domainVolume+INDEX-1);
			
			V[1]  = *(variable+INDEX-K_Increment);
			Vol[1]= *(domainVolume+INDEX-K_Increment);
			
			V[2]  = *(variable+INDEX-K_Increment-1);
			Vol[2]= *(domainVolume+INDEX-K_Increment-1);
			
			
			V[3]=*(slabVariable[3]+(NI-1)*k+i);
			Vol[3]=*(slabVolume[3]+(NI-1)*k+i);
			
			V[4]=*(slabVariable[3]+(NI-1)*k+i-1);
			Vol[4]=*(slabVolume[3]+(NI-1)*k+i-1);
			
			V[5]=*(slabVariable[3]+(NI-1)*(k-1)+i);
			Vol[5]=*(slabVolume[3]+(NI-1)*(k-1)+i);
			
			V[6]=*(slabVariable[3]+(NI-1)*(k-1)+i-1);
			Vol[6]=*(slabVolume[3]+(NI-1)*(k-1)+i-1);
			
			*(variable_vertex+indexP) = (V_p * Vol[6] + V[6] *Vol_p + Vol[0] * V[5] + Vol[5] *V[0]+Vol[1] * V[4] + Vol[4] *V[1]+Vol[2] * V[3] + Vol[3] *V[2])/(Vol_p+Vol[2]+Vol[0]+Vol[1]+Vol[3]+Vol[4]+Vol[5]+Vol[6]);
//			printf("4Inter- polated value =%e\ti=%d\tj=%d\n",*(variable_vertex+indexP),i,j);
		}
		
		j=NJ-1;/*north*/
		for(k=1;k<NK-1;k++)
		for(i=1;i<NI-1;i++)
		{
			INDEX=(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i;
			indexP =NJ*NI*k+NI*j+i;
			
			V_p  = *(variable+INDEX);
			Vol_p= *(domainVolume+INDEX);
			
			V[0]  = *(variable+INDEX-1);
			Vol[0]= *(domainVolume+INDEX-1);
			
			V[1]  = *(variable+INDEX-K_Increment);
			Vol[1]= *(domainVolume+INDEX-K_Increment);
			
			V[2]  = *(variable+INDEX-K_Increment-1);
			Vol[2]= *(domainVolume+INDEX-K_Increment-1);
			
			
			V[3]=*(slabVariable[2]+(NI-1)*k+i);
			Vol[3]=*(slabVolume[2]+(NI-1)*k+i);
			
			V[4]=*(slabVariable[2]+(NI-1)*k+i-1);
			Vol[4]=*(slabVolume[2]+(NI-1)*k+i-1);
			
			V[5]=*(slabVariable[2]+(NI-1)*(k-1)+i);
			Vol[5]=*(slabVolume[2]+(NI-1)*(k-1)+i);
			
			V[6]=*(slabVariable[2]+(NI-1)*(k-1)+i-1);
			Vol[6]=*(slabVolume[2]+(NI-1)*(k-1)+i-1);
			
			*(variable_vertex+indexP) = (V_p * Vol[6] + V[6] *Vol_p + Vol[0] * V[5] + Vol[5] *V[0]+Vol[1] * V[4] + Vol[4] *V[1]+Vol[2] * V[3] + Vol[3] *V[2])/(Vol_p+Vol[2]+Vol[0]+Vol[1]+Vol[3]+Vol[4]+Vol[5]+Vol[6]);
//			printf("5Inter- polated value =%e\ti=%d\tj=%d\n",*(variable_vertex+indexP),i,j);
		}
		
		k=0;/*bott*/
		for(j=1;j<NJ-1;j++)
		for(i=1;i<NI-1;i++)
		{
			INDEX=(NJ-1)*(NI-1)*k+(NI-1)*j+i;
			indexP =NJ*NI*k+NI*j+i;
			
			V_p  = *(variable+INDEX);
			Vol_p= *(domainVolume+INDEX);
			
			V[0]  = *(variable+INDEX-1);
			Vol[0]= *(domainVolume+INDEX-1);
			
			V[1]  = *(variable+INDEX-(NI-1));
			Vol[1]= *(domainVolume+INDEX-(NI-1));
			
			V[2]  = *(variable+INDEX-(NI-1)-1);
			Vol[2]= *(domainVolume+INDEX-(NI-1)-1);
			
			
			V[3]=*(slabVariable[5]+(NI-1)*j+i);
			Vol[3]=*(slabVolume[5]+(NI-1)*j+i);
			
			V[4]=*(slabVariable[5]+(NI-1)*j+i-1);
			Vol[4]=*(slabVolume[5]+(NI-1)*j+i-1);
			
			V[5]=*(slabVariable[5]+(NI-1)*(j-1)+i);
			Vol[5]=*(slabVolume[5]+(NI-1)*(j-1)+i);
			
			V[6]=*(slabVariable[5]+(NI-1)*(j-1)+i-1);
			Vol[6]=*(slabVolume[5]+(NI-1)*(j-1)+i-1);
			
			*(variable_vertex+indexP) = (V_p * Vol[6] + V[6] *Vol_p + Vol[0] * V[5] + Vol[5] *V[0]+Vol[1] * V[4] + Vol[4] *V[1]+Vol[2] * V[3] + Vol[3] *V[2])/(Vol_p+Vol[2]+Vol[0]+Vol[1]+Vol[3]+Vol[4]+Vol[5]+Vol[6]);
//			printf("6Inter- polated value =%e\ti=%d\tj=%d\n",*(variable_vertex+indexP),i,j);
		}
		
		k=NK-1;/*top*/
		for(j=1;j<NJ-1;j++)
		for(i=1;i<NI-1;i++)
		{
			INDEX=(NJ-1)*(NI-1)*(k-1)+(NI-1)*j+i;
			indexP =NJ*NI*k+NI*j+i;
			
			V_p  = *(variable+INDEX);
			Vol_p= *(domainVolume+INDEX);
			
			V[0]  = *(variable+INDEX-1);
			Vol[0]= *(domainVolume+INDEX-1);
			
			V[1]  = *(variable+INDEX-(NI-1));
			Vol[1]= *(domainVolume+INDEX-(NI-1));
			
			V[2]  = *(variable+INDEX-(NI-1)-1);
			Vol[2]= *(domainVolume+INDEX-(NI-1)-1);
			
			
			V[3]=*(slabVariable[4]+(NI-1)*j+i);
			Vol[3]=*(slabVolume[4]+(NI-1)*j+i);
			
			V[4]=*(slabVariable[4]+(NI-1)*j+i-1);
			Vol[4]=*(slabVolume[4]+(NI-1)*j+i-1);
			
			V[5]=*(slabVariable[4]+(NI-1)*(j-1)+i);
			Vol[5]=*(slabVolume[4]+(NI-1)*(j-1)+i);
			
			V[6]=*(slabVariable[4]+(NI-1)*(j-1)+i-1);
			Vol[6]=*(slabVolume[4]+(NI-1)*(j-1)+i-1);
			
			*(variable_vertex+indexP) = (V_p * Vol[6] + V[6] *Vol_p + Vol[0] * V[5] + Vol[5] *V[0]+Vol[1] * V[4] + Vol[4] *V[1]+Vol[2] * V[3] + Vol[3] *V[2])/(Vol_p+Vol[2]+Vol[0]+Vol[1]+Vol[3]+Vol[4]+Vol[5]+Vol[6]);
			//printf("Values possible are %e\t%e\n",(Vol_p+Vol[2]+Vol[0]+Vol[1]+Vol[3]+Vol[4]+Vol[5]+Vol[6]),(V_p * Vol[6] + V[6] *Vol_p + Vol[0] * V[5] + Vol[5] *V[0]+Vol[1] * V[4] + Vol[4] *V[1]+Vol[2] * V[3] + Vol[3] *V[2]));
			//printf("Volumes are %e\t%e\t%e\t%e\t%e\t%e\t%e\t\n",Vol_p,Vol[0],Vol[1],Vol[2],Vol[3],Vol[4],Vol[5]);
//			printf("7Inter- polated value =%e\ti=%d\tj=%d\n",*(variable_vertex+indexP),i,j);
			
		}
		//exit(12);
		i=0;j=0;
		for(k=1;k<NK-1;k++)
		{
			INDEX=(NJ-1)*(NI-1)*k+(NI-1)*j+i;
			indexP =NJ*NI*k+NI*j+i;
			
			V_p  = *(variable+INDEX);
			Vol_p= *(domainVolume+INDEX);
			
			V[0]=*(slabVariable[1]+(NJ-1)*k+j);
			Vol[0]=*(slabVolume[1]+(NJ-1)*k+j);
			
			V[1]=*(slabVariable[3]+(NI-1)*k+i);
			Vol[1]=*(slabVolume[3]+(NI-1)*k+i);
			
			V[2]=0.0;
			Vol[2]=0.0;
			
			V[3]  = *(variable+INDEX-K_Increment);
			Vol[3]= *(domainVolume+INDEX-K_Increment);
			
			
			V[4]=*(slabVariable[1]+(NJ-1)*(k-1)+j);
			Vol[4]=*(slabVolume[1]+(NJ-1)*(k-1)+j);
			
			V[5]=*(slabVariable[3]+(NI-1)*(k-1)+i);
			Vol[5]=*(slabVolume[3]+(NI-1)*(k-1)+i);
			
			V[6]=0.0;
			Vol[6]=0.0;
			
			*(variable_vertex+indexP) = (V_p * Vol[6] + V[6] *Vol_p + Vol[0] * V[5] + Vol[5] *V[0]+Vol[1] * V[4] + Vol[4] *V[1]+Vol[2] * V[3] + Vol[3] *V[2])/(Vol_p+Vol[2]+Vol[0]+Vol[1]+Vol[3]+Vol[4]+Vol[5]+Vol[6]);
//			printf("8Inter- polated value =%e\ti=%d\tj=%d\n",*(variable_vertex+indexP),i,j);
		}
		
		i=NI-1;j=0;
		for(k=1;k<NK-1;k++)
		{
			INDEX=(NJ-1)*(NI-1)*k+(NI-1)*j+i-1;
			indexP =NJ*NI*k+NI*j+i;
			
			V_p  = *(variable+INDEX);
			Vol_p= *(domainVolume+INDEX);
			
			V[0]=*(slabVariable[0]+(NJ-1)*k+j);
			Vol[0]=*(slabVolume[0]+(NJ-1)*k+j);
			
			V[1]=*(slabVariable[3]+(NI-1)*k+i-1);
			Vol[1]=*(slabVolume[3]+(NI-1)*k+i-1);
			
			V[2]=0.0;
			Vol[2]=0.0;
			
			V[3]  = *(variable+INDEX-K_Increment);
			Vol[3]= *(domainVolume+INDEX-K_Increment);
			
			
			V[4]=*(slabVariable[0]+(NJ-1)*(k-1)+j);
			Vol[4]=*(slabVolume[0]+(NJ-1)*(k-1)+j);
			
			V[5]=*(slabVariable[3]+(NI-1)*(k-1)+i-1);
			Vol[5]=*(slabVolume[3]+(NI-1)*(k-1)+i-1);
			
			V[6]=0.0;
			Vol[6]=0.0;
			
			*(variable_vertex+indexP) = (V_p * Vol[6] + V[6] *Vol_p + Vol[0] * V[5] + Vol[5] *V[0]+Vol[1] * V[4] + Vol[4] *V[1]+Vol[2] * V[3] + Vol[3] *V[2])/(Vol_p+Vol[2]+Vol[0]+Vol[1]+Vol[3]+Vol[4]+Vol[5]+Vol[6]);
//			printf("9Inter- polated value =%e\ti=%d\tj=%d\n",*(variable_vertex+indexP),i,j);
		}
		
		i=0;j=NJ-1;
		
		for(k=1;k<NK-1;k++)
		{
			INDEX=(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i;
			indexP =NJ*NI*k+NI*j+i;
			
			V_p  = *(variable+INDEX);
			Vol_p= *(domainVolume+INDEX);
			
			V[0]=*(slabVariable[1]+(NJ-1)*k+j-1);
			Vol[0]=*(slabVolume[1]+(NJ-1)*k+j-1);
			
			V[1]=*(slabVariable[2]+(NI-1)*k+i);
			Vol[1]=*(slabVolume[2]+(NI-1)*k+i);
			
			V[2]=0.0;
			Vol[2]=0.0;
			
			V[3]  = *(variable+INDEX-K_Increment);
			Vol[3]= *(domainVolume+INDEX-K_Increment);
			
			
			V[4]=*(slabVariable[1]+(NJ-1)*(k-1)+j-1);
			Vol[4]=*(slabVolume[1]+(NJ-1)*(k-1)+j-1);
			
			V[5]=*(slabVariable[2]+(NI-1)*(k-1)+i);
			Vol[5]=*(slabVolume[2]+(NI-1)*(k-1)+i);
			
			V[6]=0.0;
			Vol[6]=0.0;
			
			*(variable_vertex+indexP) = (V_p * Vol[6] + V[6] *Vol_p + Vol[0] * V[5] + Vol[5] *V[0]+Vol[1] * V[4] + Vol[4] *V[1]+Vol[2] * V[3] + Vol[3] *V[2])/(Vol_p+Vol[2]+Vol[0]+Vol[1]+Vol[3]+Vol[4]+Vol[5]+Vol[6]);
//			printf("10Inter- polated value =%e\ti=%d\tj=%d\n",*(variable_vertex+indexP),i,j);
		}
		i=NI-1;j=NJ-1;
		for(k=1;k<NK-1;k++)
		{
			INDEX=(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i-1;
			indexP =NJ*NI*k+NI*j+i;
			
			V_p  = *(variable+INDEX);
			Vol_p= *(domainVolume+INDEX);
			
			V[0]=*(slabVariable[0]+(NJ-1)*k+j-1);
			Vol[0]=*(slabVolume[0]+(NJ-1)*k+j-1);
			
			V[1]=*(slabVariable[2]+(NI-1)*k+i-1);
			Vol[1]=*(slabVolume[2]+(NI-1)*k+i-1);
			
			V[2]=0.0;
			Vol[2]=0.0;
			
			V[3]  = *(variable+INDEX-K_Increment);
			Vol[3]= *(domainVolume+INDEX-K_Increment);
			
			
			V[4]=*(slabVariable[0]+(NJ-1)*(k-1)+j-1);
			Vol[4]=*(slabVolume[0]+(NJ-1)*(k-1)+j-1);
			
			V[5]=*(slabVariable[2]+(NI-1)*(k-1)+i-1);
			Vol[5]=*(slabVolume[2]+(NI-1)*(k-1)+i-1);
			
			V[6]=0.0;
			Vol[6]=0.0;
			
			*(variable_vertex+indexP) = (V_p * Vol[6] + V[6] *Vol_p + Vol[0] * V[5] + Vol[5] *V[0]+Vol[1] * V[4] + Vol[4] *V[1]+Vol[2] * V[3] + Vol[3] *V[2])/(Vol_p+Vol[2]+Vol[0]+Vol[1]+Vol[3]+Vol[4]+Vol[5]+Vol[6]);
//			printf("11Inter- polated value =%e\ti=%d\tj=%d\n",*(variable_vertex+indexP),i,j);
		}
		
		i=0;k=0;
		
		for(j=1;j<NJ-1;j++)
		{
			INDEX=(NJ-1)*(NI-1)*k+(NI-1)*j+i;
			indexP =NJ*NI*k+NI*j+i;
			
			V_p  = *(variable+INDEX);
			Vol_p= *(domainVolume+INDEX);
			
			V[0]=*(slabVariable[1]+(NJ-1)*k+j);
			Vol[0]=*(slabVolume[1]+(NJ-1)*k+j);
			
			V[1]=*(slabVariable[5]+(NI-1)*j+i);
			Vol[1]=*(slabVolume[5]+(NI-1)*j+i);
			
			V[2]=0.0;
			Vol[2]=0.0;
			
			V[3]  = *(variable+INDEX-(NI-1));
			Vol[3]= *(domainVolume+INDEX-(NI-1));
			
			
			V[4]=*(slabVariable[1]+(NJ-1)*k+j-1);
			Vol[4]=*(slabVolume[1]+(NJ-1)*k+j-1);
			
			V[5]=*(slabVariable[5]+(NI-1)*(j-1)+i);
			Vol[5]=*(slabVolume[5]+(NI-1)*(j-1)+i);
			
			V[6]=0.0;
			Vol[6]=0.0;
			
			*(variable_vertex+indexP) = (V_p * Vol[6] + V[6] *Vol_p + Vol[0] * V[5] + Vol[5] *V[0]+Vol[1] * V[4] + Vol[4] *V[1]+Vol[2] * V[3] + Vol[3] *V[2])/(Vol_p+Vol[2]+Vol[0]+Vol[1]+Vol[3]+Vol[4]+Vol[5]+Vol[6]);
//			printf("12Inter- polated value =%e\ti=%d\tj=%d\n",*(variable_vertex+indexP),i,j);
		}
		i=NI-1;k=0;
		
		for(j=1;j<NJ-1;j++)
		{
			INDEX=(NJ-1)*(NI-1)*k+(NI-1)*j+i-1;
			indexP =NJ*NI*k+NI*j+i;
			
			V_p  = *(variable+INDEX);
			Vol_p= *(domainVolume+INDEX);
			
			V[0]=*(slabVariable[0]+(NJ-1)*k+j);
			Vol[0]=*(slabVolume[0]+(NJ-1)*k+j);
			
			V[1]=*(slabVariable[5]+(NI-1)*j+i-1);
			Vol[1]=*(slabVolume[5]+(NI-1)*j+i-1);
			
			V[2]=0.0;
			Vol[2]=0.0;
			
			V[3]  = *(variable+INDEX-(NI-1));
			Vol[3]= *(domainVolume+INDEX-(NI-1));
			
			
			V[4]=*(slabVariable[0]+(NJ-1)*k+j-1);
			Vol[4]=*(slabVolume[0]+(NJ-1)*k+j-1);
			
			V[5]=*(slabVariable[5]+(NI-1)*(j-1)+i-1);
			Vol[5]=*(slabVolume[5]+(NI-1)*(j-1)+i-1);
			
			V[6]=0.0;
			Vol[6]=0.0;
			
			*(variable_vertex+indexP) = (V_p * Vol[6] + V[6] *Vol_p + Vol[0] * V[5] + Vol[5] *V[0]+Vol[1] * V[4] + Vol[4] *V[1]+Vol[2] * V[3] + Vol[3] *V[2])/(Vol_p+Vol[2]+Vol[0]+Vol[1]+Vol[3]+Vol[4]+Vol[5]+Vol[6]);
//			printf("13Inter- polated value =%e\ti=%d\tj=%d\n",*(variable_vertex+indexP),i,j);
		}
		
		i=0;k=NK-1;
		
		for(j=1;j<NJ-1;j++)
		{
			INDEX=(NJ-1)*(NI-1)*(k-1)+(NI-1)*j+i;
			indexP =NJ*NI*k+NI*j+i;
			
			V_p  = *(variable+INDEX);
			Vol_p= *(domainVolume+INDEX);
			
			V[0]=*(slabVariable[1]+(NJ-1)*(k-1)+j);
			Vol[0]=*(slabVolume[1]+(NJ-1)*(k-1)+j);
			
			V[1]=*(slabVariable[4]+(NI-1)*j+i);
			Vol[1]=*(slabVolume[4]+(NI-1)*j+i);
			
			V[2]=0.0;
			Vol[2]=0.0;
			
			V[3]  = *(variable+INDEX-(NI-1));
			Vol[3]= *(domainVolume+INDEX-(NI-1));
			
			
			V[4]=*(slabVariable[1]+(NJ-1)*(k-1)+j-1);
			Vol[4]=*(slabVolume[1]+(NJ-1)*(k-1)+j-1);
			
			V[5]=*(slabVariable[4]+(NI-1)*(j-1)+i);
			Vol[5]=*(slabVolume[4]+(NI-1)*(j-1)+i);
			
			V[6]=0.0;
			Vol[6]=0.0;
			
			*(variable_vertex+indexP) = (V_p * Vol[6] + V[6] *Vol_p + Vol[0] * V[5] + Vol[5] *V[0]+Vol[1] * V[4] + Vol[4] *V[1]+Vol[2] * V[3] + Vol[3] *V[2])/(Vol_p+Vol[2]+Vol[0]+Vol[1]+Vol[3]+Vol[4]+Vol[5]+Vol[6]);
//			printf("14Inter- polated value =%e\ti=%d\tj=%d\n",*(variable_vertex+indexP),i,j);
		}
		i=NI-1;k=NK-1;
		
		for(j=1;j<NJ-1;j++)
		{
			INDEX=(NJ-1)*(NI-1)*(k-1)+(NI-1)*j+i-1;
			indexP =NJ*NI*k+NI*j+i;
			
			V_p  = *(variable+INDEX);
			Vol_p= *(domainVolume+INDEX);
			
			V[0]=*(slabVariable[0]+(NJ-1)*(k-1)+j);
			Vol[0]=*(slabVolume[0]+(NJ-1)*(k-1)+j);
			
			V[1]=*(slabVariable[4]+(NI-1)*j+i-1);
			Vol[1]=*(slabVolume[4]+(NI-1)*j+i-1);
			
			V[2]=0.0;
			Vol[2]=0.0;
			
			V[3]  = *(variable+INDEX-(NI-1));
			Vol[3]= *(domainVolume+INDEX-(NI-1));
			
			
			V[4]=*(slabVariable[0]+(NJ-1)*(k-1)+j-1);
			Vol[4]=*(slabVolume[0]+(NJ-1)*(k-1)+j-1);
			
			V[5]=*(slabVariable[4]+(NI-1)*(j-1)+i-1);
			Vol[5]=*(slabVolume[4]+(NI-1)*(j-1)+i-1);
			
			V[6]=0.0;
			Vol[6]=0.0;
			
			*(variable_vertex+indexP) = (V_p * Vol[6] + V[6] *Vol_p + Vol[0] * V[5] + Vol[5] *V[0]+Vol[1] * V[4] + Vol[4] *V[1]+Vol[2] * V[3] + Vol[3] *V[2])/(Vol_p+Vol[2]+Vol[0]+Vol[1]+Vol[3]+Vol[4]+Vol[5]+Vol[6]);
//			printf("15Inter- polated value =%e\ti=%d\tj=%d\n",*(variable_vertex+indexP),i,j);
		}
		
		j=0;k=0;
		for(i=1;i<NI-1;i++)
		{
			INDEX=(NJ-1)*(NI-1)*k+(NI-1)*j+i;
			indexP =NJ*NI*k+NI*j+i;
			
			V_p  = *(variable+INDEX);
			Vol_p= *(domainVolume+INDEX);
			
			V[0]=*(slabVariable[3]+(NI-1)*k+i);
			Vol[0]=*(slabVolume[3]+(NI-1)*k+i);
			
			V[1]=*(slabVariable[5]+(NI-1)*j+i);
			Vol[1]=*(slabVolume[5]+(NI-1)*j+i);
			
			V[2]=0.0;
			Vol[2]=0.0;
			
			V[3]  = *(variable+INDEX-(NI-1));
			Vol[3]= *(domainVolume+INDEX-(NI-1));
			
			
			V[4]=*(slabVariable[3]+(NI-1)*k+i-1);
			Vol[4]=*(slabVolume[3]+(NI-1)*k+i-1);
			
			V[5]=*(slabVariable[5]+(NI-1)*j+i-1);
			Vol[5]=*(slabVolume[5]+(NI-1)*j+i-1);
			
			V[6]=0.0;
			Vol[6]=0.0;
			
			*(variable_vertex+indexP) = (V_p * Vol[6] + V[6] *Vol_p + Vol[0] * V[5] + Vol[5] *V[0]+Vol[1] * V[4] + Vol[4] *V[1]+Vol[2] * V[3] + Vol[3] *V[2])/(Vol_p+Vol[2]+Vol[0]+Vol[1]+Vol[3]+Vol[4]+Vol[5]+Vol[6]);
//			printf("16Inter- polated value =%e\ti=%d\tj=%d\n",*(variable_vertex+indexP),i,j);
		}
		
		
		j=NJ-1;k=0;
		for(i=1;i<NI-1;i++)
		{
			INDEX=(NJ-1)*(NI-1)*k+(NI-1)*(j-1)+i;
			indexP =NJ*NI*k+NI*j+i;
			
			V_p  = *(variable+INDEX);
			Vol_p= *(domainVolume+INDEX);
			
			V[0]=*(slabVariable[2]+(NI-1)*k+i);
			Vol[0]=*(slabVolume[2]+(NI-1)*k+i);
			
			V[1]=*(slabVariable[5]+(NI-1)*(j-1)+i);
			Vol[1]=*(slabVolume[5]+(NI-1)*(j-1)+i);
			
			V[2]=0.0;
			Vol[2]=0.0;
			
			V[3]  = *(variable+INDEX-(NI-1));
			Vol[3]= *(domainVolume+INDEX-(NI-1));
			
			
			V[4]=*(slabVariable[2]+(NI-1)*k+i-1);
			Vol[4]=*(slabVolume[2]+(NI-1)*k+i-1);
			
			V[5]=*(slabVariable[5]+(NI-1)*(j-1)+i-1);
			Vol[5]=*(slabVolume[5]+(NI-1)*(j-1)+i-1);
			
			V[6]=0.0;
			Vol[6]=0.0;
			
			*(variable_vertex+indexP) = (V_p * Vol[6] + V[6] *Vol_p + Vol[0] * V[5] + Vol[5] *V[0]+Vol[1] * V[4] + Vol[4] *V[1]+Vol[2] * V[3] + Vol[3] *V[2])/(Vol_p+Vol[2]+Vol[0]+Vol[1]+Vol[3]+Vol[4]+Vol[5]+Vol[6]);
//			printf("17Inter- polated value =%e\ti=%d\tj=%d\n",*(variable_vertex+indexP),i,j);
		}
		
		j=0;k=NK-1;
		for(i=1;i<NI-1;i++)
		{
			INDEX=(NJ-1)*(NI-1)*(k-1)+(NI-1)*j+i;
			indexP =NJ*NI*k+NI*j+i;
			
			V_p  = *(variable+INDEX);
			Vol_p= *(domainVolume+INDEX);
			
			V[0]=*(slabVariable[3]+(NI-1)*(k-1)+i);
			Vol[0]=*(slabVolume[3]+(NI-1)*(k-1)+i);
			
			V[1]=*(slabVariable[4]+(NI-1)*j+i);
			Vol[1]=*(slabVolume[4]+(NI-1)*j+i);
			
			V[2]=0.0;
			Vol[2]=0.0;
			
			V[3]  = *(variable+INDEX-(NI-1));
			Vol[3]= *(domainVolume+INDEX-(NI-1));
			
			
			V[4]=*(slabVariable[3]+(NI-1)*(k-1)+i-1);
			Vol[4]=*(slabVolume[3]+(NI-1)*(k-1)+i-1);
			
			V[5]=*(slabVariable[4]+(NI-1)*j+i-1);
			Vol[5]=*(slabVolume[4]+(NI-1)*j+i-1);
			
			V[6]=0.0;
			Vol[6]=0.0;
			
			*(variable_vertex+indexP) = (V_p * Vol[6] + V[6] *Vol_p + Vol[0] * V[5] + Vol[5] *V[0]+Vol[1] * V[4] + Vol[4] *V[1]+Vol[2] * V[3] + Vol[3] *V[2])/(Vol_p+Vol[2]+Vol[0]+Vol[1]+Vol[3]+Vol[4]+Vol[5]+Vol[6]);
//			printf("18Inter- polated value =%e\ti=%d\tj=%d\n",*(variable_vertex+indexP),i,j);
		}
		j=NJ-1;k=NK-1;
		for(i=1;i<NI-1;i++)
		{
			INDEX=(NJ-1)*(NI-1)*(k-1)+(NI-1)*(j-1)+i;
			indexP =NJ*NI*k+NI*j+i;
			
			V_p  = *(variable+INDEX);
			Vol_p= *(domainVolume+INDEX);
			
			V[0]=*(slabVariable[2]+(NI-1)*(k-1)+i);
			Vol[0]=*(slabVolume[2]+(NI-1)*(k-1)+i);
			
			V[1]=*(slabVariable[4]+(NI-1)*(j-1)+i);
			Vol[1]=*(slabVolume[4]+(NI-1)*(j-1)+i);
			
			V[2]=0.0;
			Vol[2]=0.0;
			
			V[3]  = *(variable+INDEX-(NI-1));
			Vol[3]= *(domainVolume+INDEX-(NI-1));
			
			
			V[4]=*(slabVariable[2]+(NI-1)*(k-1)+i-1);
			Vol[4]=*(slabVolume[2]+(NI-1)*(k-1)+i-1);
			
			V[5]=*(slabVariable[4]+(NI-1)*(j-1)+i-1);
			Vol[5]=*(slabVolume[4]+(NI-1)*(j-1)+i-1);
			
			V[6]=0.0;
			Vol[6]=0.0;
			
			*(variable_vertex+indexP) = (V_p * Vol[6] + V[6] *Vol_p + Vol[0] * V[5] + Vol[5] *V[0]+Vol[1] * V[4] + Vol[4] *V[1]+Vol[2] * V[3] + Vol[3] *V[2])/(Vol_p+Vol[2]+Vol[0]+Vol[1]+Vol[3]+Vol[4]+Vol[5]+Vol[6]);
//			printf("19Inter- polated value =%e\ti=%d\tj=%d\n",*(variable_vertex+indexP),i,j);
		}
		
		i=0;j=0;k=0;
		INDEX=(NJ-2)*(NI-2)*k+(NI-2)*j+i;
		indexP =(NJ-1)*(NI-1)*k+(NI-1)*j+i;
		*(variable_vertex+indexP)= *(variable+INDEX);
		
		i=NI-1;j=0;k=0;
		INDEX=(NJ-2)*(NI-2)*k+(NI-2)*j+i-1;
		indexP =(NJ-1)*(NI-1)*k+(NI-1)*j+i;
		*(variable_vertex+indexP)= *(variable+INDEX);
		
		i=0;j=NJ-1;k=0;
		INDEX=(NJ-2)*(NI-2)*k+(NI-2)*(j-1)+i;
		indexP =(NJ-1)*(NI-1)*k+(NI-1)*j+i;
		*(variable_vertex+indexP)= *(variable+INDEX);
		
		i=0;j=0;k=NK-1;
		INDEX=(NJ-2)*(NI-2)*(k-1)+(NI-2)*j+i;
		indexP =(NJ-1)*(NI-1)*k+(NI-1)*j+i;
		*(variable_vertex+indexP)= *(variable+INDEX);
		
		i=NI-1;j=NJ-1;k=0;
		INDEX=(NJ-2)*(NI-2)*k+(NI-2)*(j-1)+i-1;
		indexP =(NJ-1)*(NI-1)*k+(NI-1)*j+i;
		*(variable_vertex+indexP)= *(variable+INDEX);
		
		i=NI-1;j=0;k=NK-1;
		INDEX=(NJ-2)*(NI-2)*(k-1)+(NI-2)*j+i-1;
		indexP =(NJ-1)*(NI-1)*k+(NI-1)*j+i;
		*(variable_vertex+indexP)= *(variable+INDEX);
		
		i=0;j=NJ-1;k=NK-1;
		INDEX=(NJ-2)*(NI-2)*(k-1)+(NI-2)*(j-1)+i;
		indexP =(NJ-1)*(NI-1)*k+(NI-1)*j+i;
		*(variable_vertex+indexP)= *(variable+INDEX);
		
		i=NI-1;j=NJ-1;k=NK-1;
		INDEX=(NJ-2)*(NI-2)*(k-1)+(NI-2)*(j-1)+i-1;
		indexP =(NJ-1)*(NI-1)*k+(NI-1)*j+i;
		*(variable_vertex+indexP)= *(variable+INDEX);	
		
		/*indexP=0;
		for (k=1;k<NK-1;k++)
		for (j=1;j<NJ-1;j++)
		for (i=1;i<NI-1;i++)
		{
			printf("%.2lf\ti=%d\tj=%d\tk=%d\n",*(variable_vertex+indexP),i,j,k);	
			indexP++;
		}
		exit(12);*/	
	}
	return (0.0);
}

