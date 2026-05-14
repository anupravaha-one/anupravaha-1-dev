#include "APKconvfluxC3D.h"
/**
The function convfluxC3D() is the main control file for calculating convective fluxes when the spatial discretisation
is Central Difference and domain is 3 Dimensional. This function passes the control to the functions ConvfluxIntC3D(),
ConvfluxExtCorner3D() and ConvfluxExtCother3D() depending on the cell location within the zone. If the cells are completely
interior cells then ConvfluxIntC3D() is called. If a cell has only one boundary or an interface among any of its six 
faces(since 3D), then ConvfluxExtCother3D() is called. If a cell has more than one boundary or an interface among any of its
six faces(since 3D), then ConvfluxExtCorner3D() is called (This condition occurs for completely corner cells).
*/
int convfluxC3D(double *phicellold,double *phicellnew,double alpha,Zone *zone,
		double *Fc_n,double *Fn_d,double *convdia,int vartag,int scheme,int IndexConv[3])
{
	//printf("I am in convflux\n");
	
	int NI,NJ,NK,NcellI,NcellJ,NcellK,k,j,i;   
	int face,index[3];
	int *ptrSize;
	double *volume,*F;
	double  *slabPhi[6],*slabvolume[6],NCoeff[6],Diagonal[6];
	double diagonal,ddiag;
		
	ptrSize=zone->Size;
	
		NI=*(ptrSize+0);
		NJ=*(ptrSize+1);
		NK=*(ptrSize+2);
		NcellI=NI-1;
		NcellJ=NJ-1;
		NcellK=NK-1;
		volume=zone->ZGeom->Volume;
		/*Fc_n=(double*)malloc(6*NcellK*NcellJ*NcellI*sizeof(double));
		convdia=(double *)malloc(NcellK*NcellJ*NcellI*sizeof(double));
		Fn_d=(double*)malloc(6*NcellK*NcellJ*NcellI*sizeof(double));*/
		if(scheme==0)
		{
			if(vartag == 9)
				F=zone->Flowsol->Electricflux;//Electric Field Flux
			else
				F=zone->Flowsol->Massflux;
		}
		else
		{
			if(vartag == 9)
				F=zone->Flowsol->Electricflux;//Electric Field Flux
			else
				F=zone->PrevFlowsol->Massflux;
		}
		//printf("i am outside face loop");
		MakeSlabPtrArray(zone,slabPhi,vartag);/*Called for making slabs of u */
		MakeSlabPtrArray(zone,slabvolume,5);
		//printf("pointer=%u\n",slabvolume[face]);
		
			k=IndexConv[0],j=IndexConv[1],i=IndexConv[2];
			index[0]=k,index[1]=j,index[2]=i;
			//printf("\n real routine called ");
			
		if((i>0 && i<NI-2) && (j>0 && j<NJ-2) && (k>0 && k<NK-2))
		{
		
		ConvfluxIntC3D(zone,F,phicellnew,phicellold,volume,Fc_n,Fn_d,index,NCoeff,Diagonal,slabvolume,slabPhi,alpha);
			
			diagonal=0.0;
			
			for(face=0;face<6;face++)
			{
				diagonal=diagonal+Diagonal[face];
				//*(Ndia+((NI-1)*(NJ-1)*k+(NI-1)*j+i)*6+face)=NCoeff[face];
			}
			*(convdia+(NI-1)*(NJ-1)*k+(NI-1)*j+i)=diagonal;
			//printf("Convdia=%lf\ti=%d\tj=%d\tk=%d\n",*(convdia+(NI-1)*(NJ-1)*k+(NI-1)*j+i),i,j,k);
			
			/*if(*(convdia+(NI-1)*(NJ-1)*k+(NI-1)*j+i)<-10)
			{	exit(1);
			}*/
		}
		else
		{
			
			diagonal=0.0;
		for(face=0;face<6;face++)
		{

			switch(face)
			{   
				case 0:/*Imin->West Face*/
					if(i==0)
					{
					ddiag=ConvfluxExtCorner3D(zone,F,phicellnew,phicellold,volume,Fc_n,Fn_d,index,NCoeff,slabvolume,slabPhi,face,alpha);
					Diagonal[face]=ddiag;
					}/*allocate memory convdia and Fn_d*/
					else
					{
					ddiag=ConvfluxExtCother3D(zone,F,phicellnew,phicellold,volume,Fc_n,Fn_d,index,NCoeff,slabvolume,slabPhi,face,alpha);
					Diagonal[face]=ddiag;
					//diagonal=diagonal+Diagonal[face];
					}
							
					break;

				case 1:/*Imax->East Face*/
					
					if(i==NI-2)
					{
					ddiag=ConvfluxExtCorner3D(zone,F,phicellnew,phicellold,volume,Fc_n,Fn_d,index,NCoeff,slabvolume,slabPhi,face,alpha);
					Diagonal[face]=ddiag;
					//diagonal=diagonal+Diagonal[face];
					//printf("diagonal=%lf\n",Diagonal[face]);
					}
					else
					{
					ddiag=ConvfluxExtCother3D(zone,F,phicellnew,phicellold,volume,Fc_n,Fn_d,index,NCoeff,slabvolume,slabPhi,face,alpha);
					Diagonal[face]=ddiag;
					//diagonal=diagonal+Diagonal[face];
					}	
						
					break;

				case 2:/*Jmin->South face*/
					
					if(j==0)
					{
					ddiag=ConvfluxExtCorner3D(zone,F,phicellnew,phicellold,volume,Fc_n,Fn_d,index,NCoeff,slabvolume,slabPhi,face,alpha);	
					Diagonal[face]=ddiag;
					//diagonal=diagonal+Diagonal[face];
					
					}
					else
					{
					ddiag=ConvfluxExtCother3D(zone,F,phicellnew,phicellold,volume,Fc_n,Fn_d,index,NCoeff,slabvolume,slabPhi,face,alpha);
					Diagonal[face]=ddiag;
					//diagonal=diagonal+Diagonal[face];
					}
								
					break;

					case 3:/*Jmax->North Face*/
					
					if(j==NJ-2)
					{
					ddiag=ConvfluxExtCorner3D(zone,F,phicellnew,phicellold,volume,Fc_n,Fn_d,index,NCoeff,slabvolume,slabPhi,face,alpha);	
					Diagonal[face]=ddiag;
					}
					else
					{
					ddiag=ConvfluxExtCother3D(zone,F,phicellnew,phicellold,volume,Fc_n,Fn_d,index,NCoeff,slabvolume,slabPhi,face,alpha);			
					Diagonal[face]=ddiag;
					//diagonal=diagonal+Diagonal[face];
					}
					
					break;

				case 4:/*Kmin->Bottam Face*/
					
					if(k==0)
					{
					ddiag=ConvfluxExtCorner3D(zone,F,phicellnew,phicellold,volume,Fc_n,Fn_d,index,NCoeff,slabvolume,slabPhi,face,alpha);
					Diagonal[face]=ddiag;
					//diagonal=diagonal+Diagonal[face];
					}
					else
					{
					ddiag=ConvfluxExtCother3D(zone,F,phicellnew,phicellold,volume,Fc_n,Fn_d,index,NCoeff,slabvolume,slabPhi,face,alpha);			
					Diagonal[face]=ddiag;
					//diagonal=diagonal+Diagonal[face];
					}
							
					break;
				case 5:/*Kmax->Top face*/
					
					if(k==NK-2)
					{
					ddiag=ConvfluxExtCorner3D(zone,F,phicellnew,phicellold,volume,Fc_n,Fn_d,index,NCoeff,slabvolume,slabPhi,face,alpha);
					Diagonal[face]=ddiag;
					}
					else
					{
					ddiag=ConvfluxExtCother3D(zone,F,phicellnew,phicellold,volume,Fc_n,Fn_d,index,NCoeff,slabvolume,slabPhi,face,alpha);
					Diagonal[face]=ddiag;
					}
							
					break;
			}
			//*(CoeffN+((NJ-1)*(NI-1)*k+(NI-1)*j+i)*6+face)=NCoeff[face];
			
			diagonal=diagonal+Diagonal[face];
			//*(Ndia+((NI-1)*(NJ-1)*k+(NI-1)*j+i)*6+face)=NCoeff[face];	
				
		}
		*(convdia+(NJ-1)*(NI-1)*k+(NI-1)*j+i)=diagonal;
		//printf("Convdia=%lf\ti=%d\tj=%d\tk=%d\n",*(convdia+(NI-1)*(NJ-1)*k+(NI-1)*j+i),i,j,k);
		
		
		}
		
	
	
		return(0);
}
