#include "APKCallCalcCoefficientsForVariables.h"

/**
This function calls respective CalcCoefficients functions depending on constant and varying diffusion coeffients cases 
 */

double CallCalcCoefficientsForVariables(Base *base,int g,int vartag)
{
	int NI,NJ,NK,i,j,k,indexP,*ptrSize,flagrho,flagdiffusionCoeff,flagRhoDisperse,che;
	double *diffdia,*CoeffN[18],*slabVolume[6],*slabRho[6],*Domain;
	double *beta,*volumeDomain;
	double Vol[18];
	double *diffdia_GEOM,*CoeffN_GEOM[18];
	
	double diff_coeff;
	double r[18],r_p;
	double *diffusionCoeffdomain,*slabDiffusionCoeff[6];
	double diffCoeff_p,diffCoeff_e,diffCoeff_w,diffCoeff_s,diffCoeff_n,diffCoeff_t,diffCoeff_b;//These are Cell center values
	//double diffCoeff_Eface,diffCoeff_Wface,diffCoeff_Sface,diffCoeff_Nface,diffCoeff_Tface,diffCoeff_Bface;
	double diffCoeff_face[6];
	double vol_p,*Diag[6],*D[36],*DN[24];
	int face,K_Increment,J_Increment;
	//chandrasekhar-Turbulent gas-particle
	double *Domain_Vartheta,*slabVarthetaDisperse[6],v[18],v_p,rd[18],rd_p;
	double *Domain_RhoDisperse,*slabRhoDisperse[6];
	double EPSILONU;
	EPSILONU = base->BaseParameters->EPSILONDrop;
	Zone *zone;
	int iZone,numZones;
	numZones=base->Number_of_zones;
	if(base->flagMPFV==0)
	{	
		for (iZone=0;iZone<numZones;iZone++)
		{		
		
			zone = base->CgnsZone+iZone;	
			if((base->flagMPFV==1)&&(vartag==4))//Is there any use of this???--chandrasekhar
			{	
				diffdia = (zone->grid+g)->DiffusionCoeff_P;
				for(i=0;i<18;i++)
				{
					CoeffN[i] = (zone->grid+g)->DiffusionCoeff_nb[i];
				}
				
			}
			else
			{
				diffdia = (zone->grid+g)->DiffusionCoeffPhi_P;
				
				for(i=0;i<18;i++)
				{
					CoeffN[i] = (zone->grid+g)->DiffusionCoeffPhi_nb[i];
				}
				for(i=0;i<6;i++)
				{
					Diag[i] = (zone->grid+g)->Diffusion_Diag[i];
				}
				for(i=0;i<36;i++)
				{
					D[i] = (zone->grid+g)->Diffusion_nb[i];
				}
				for(i=0;i<24;i++)
				{
					DN[i] = (zone->grid+g)->Diffusion_Nonortho[i];
				}
			}
			
			
			ptrSize = zone->Size ;
			beta=zone->ZGeom->Beta;
			volumeDomain = zone->ZGeom->Volume;
			
			NI=*(ptrSize+0);
			NJ=*(ptrSize+1);
			NK=*(ptrSize+2);
			
			K_Increment = (NI-1) *(NJ-1);
			J_Increment = NI-1;
			flagdiffusionCoeff = MakeDiffusionCoeffFlag(zone,vartag);
			//printf("\nZone-%s$$Vartag=%d $$flagdiffusionCoeff=%d",zone->zonename,vartag,flagdiffusionCoeff);
			MakeSlabPtrArray(zone,slabVolume,5);/*Called for making slabs of Volumes*/
			flagrho = zone->Flowsol->flagRho;
			
			if(base->flagParticle==1)//chandrasekhar--Turbulent gas-particle
			{
				Domain_Vartheta=MakeDomainArray(zone,650);
				MakeSlabPtrArray(zone,slabVarthetaDisperse,650);
				Domain_RhoDisperse=MakeDomainArray(zone,66);
				MakeSlabPtrArray(zone,slabRhoDisperse,66);
				flagRhoDisperse=base->BaseParameters->flagRhoDisperse;
			}
				
			///This loop is added to avoid unnecessary accessing the volumes and other multiplications in case of a variable with Constant diffusion coefficient.  Still this extension is to be made to variable with variable diffusion coefficient.  This is not been done right now as the extension needs more time and memory.  This needs to be discussed with sir before finalizing this subroutine.
	
			if(flagdiffusionCoeff==0)
			{
				diffusionCoeffdomain = MakeDiffusionCoeffDomain(zone,vartag);
				diff_coeff = *diffusionCoeffdomain;
				
				diffdia_GEOM = (zone->grid+g)->DiffusionCoeff_P;
				
				for(i=0;i<18;i++)
				{
					CoeffN_GEOM[i] = (zone->grid+g)->DiffusionCoeff_nb[i];
				}
				indexP =0;
				for(k=0;k<NK-1;k++)
					for(j=0;j<NJ-1;j++)
						for(i=0;i<NI-1;i++)
				{
					*(diffdia  +indexP) = diff_coeff * *(diffdia_GEOM  +indexP);
					*(CoeffN[0]+indexP) = diff_coeff * *(CoeffN_GEOM[0]+indexP);
					*(CoeffN[1]+indexP) = diff_coeff * *(CoeffN_GEOM[1]+indexP);
					*(CoeffN[2]+indexP) = diff_coeff * *(CoeffN_GEOM[2]+indexP);
					*(CoeffN[3]+indexP) = diff_coeff * *(CoeffN_GEOM[3]+indexP);
					*(CoeffN[4]+indexP) = diff_coeff * *(CoeffN_GEOM[4]+indexP);
					*(CoeffN[5]+indexP) = diff_coeff * *(CoeffN_GEOM[5]+indexP);
					*(CoeffN[6]+indexP) = diff_coeff * *(CoeffN_GEOM[6]+indexP);
					*(CoeffN[7]+indexP) = diff_coeff * *(CoeffN_GEOM[7]+indexP);
					*(CoeffN[8]+indexP) = diff_coeff * *(CoeffN_GEOM[8]+indexP);
					*(CoeffN[9]+indexP) = diff_coeff * *(CoeffN_GEOM[9]+indexP);
					*(CoeffN[10]+indexP) = diff_coeff * *(CoeffN_GEOM[10]+indexP);
					*(CoeffN[11]+indexP) = diff_coeff * *(CoeffN_GEOM[11]+indexP);
					*(CoeffN[12]+indexP) = diff_coeff * *(CoeffN_GEOM[12]+indexP);
					*(CoeffN[13]+indexP) = diff_coeff * *(CoeffN_GEOM[13]+indexP);
					*(CoeffN[14]+indexP) = diff_coeff * *(CoeffN_GEOM[14]+indexP);
					*(CoeffN[15]+indexP) = diff_coeff * *(CoeffN_GEOM[15]+indexP);
					*(CoeffN[16]+indexP) = diff_coeff * *(CoeffN_GEOM[16]+indexP);
					*(CoeffN[17]+indexP) = diff_coeff * *(CoeffN_GEOM[17]+indexP);
					indexP++;
				}
			}
			else
			{
				diffusionCoeffdomain = MakeDiffusionCoeffDomain(zone,vartag);
				MakeDiffusionCoeffSlabPtrArray(zone,slabDiffusionCoeff,vartag);
// 				printf("The [zone--%s]--[diff value --%e]--[vartag--%d]\n",zone->zonename,*diffusionCoeffdomain,vartag);
// 				printf("The slabDiffusionCoeff at 0th slab in zone %s is %e\n",zone->zonename,*slabDiffusionCoeff[0]);
// 				printf("The slabDiffusionCoeff at 1st slab in zone %s is %e\n",zone->zonename,*slabDiffusionCoeff[1]);
// 				printf("The slabDiffusionCoeff at 2nd slab in zone %s is %e\n",zone->zonename,*slabDiffusionCoeff[2]);
// 				printf("The slabDiffusionCoeff at 3rd slab in zone %s is %e\n",zone->zonename,*slabDiffusionCoeff[3]);
// 				printf("The slabDiffusionCoeff at 4th slab in zone %s is %e\n",zone->zonename,*slabDiffusionCoeff[4]);
// 				printf("The slabDiffusionCoeff at 5th slab in zone %s is %e\n",zone->zonename,*slabDiffusionCoeff[5]);
				if((base->BaseParameters->flag2d == 1)&&((vartag==62)||(vartag==63)||(vartag==64)))//chandrasekhar--turbulent gas-particle
				{	
					indexP=0;//printf("\nI am entering 1 and vartag=%d",vartag);
					for(k=0;k<NK-1;k++)
					for(j=0;j<NJ-1;j++)
					for(i=0;i<NI-1;i++)
					{
						if(i!=0)
						{
							Vol[0]=*(volumeDomain+indexP-1);
							diffCoeff_w=*(diffusionCoeffdomain+indexP-1);
							
							if((*(Domain_Vartheta+indexP-1))>EPSILONU)
								v[0] = 1.0/(*(Domain_Vartheta+indexP-1));
							else
								v[0]=0.0;
                            
							rd[0] = 1.0/(*(Domain_RhoDisperse+(indexP-1)*flagRhoDisperse));	
						}
						else
						{
							Vol[0] =*(slabVolume[1]+j);
							diffCoeff_w =*(slabDiffusionCoeff[1]+j);
							
							if((*(slabVarthetaDisperse[1]+j))>EPSILONU)
								v[0] =1.0/(*(slabVarthetaDisperse[1]+j));
							else
								v[0]=0.0;
							
							rd[0] =1.0/(*(slabRhoDisperse[1]+j*flagRhoDisperse));	
						}
						if(i!=(NI-2))
						{										
							Vol[1]=*(volumeDomain+indexP+1);
							diffCoeff_e=*(diffusionCoeffdomain+indexP+1);
							
							if((*(Domain_Vartheta+indexP+1))>EPSILONU)
								v[1]=1.0/(*(Domain_Vartheta+indexP+1));
							else
								v[1]=0.0;
							
							rd[1] = 1.0/(*(Domain_RhoDisperse+(indexP+1)*flagRhoDisperse));
						}
						else
						{
							Vol[1]= *(slabVolume[0]+j);
							diffCoeff_e=*(slabDiffusionCoeff[0]+j);
							
							if((*(slabVarthetaDisperse[0]+j))>EPSILONU)
								v[1]=1.0/(*(slabVarthetaDisperse[0]+j));
							else
								v[1]=0.0;
							
							rd[1]=1.0/(*(slabRhoDisperse[0]+j*flagRhoDisperse));
						}

						if(j!=0)
						{
							Vol[2]=*(volumeDomain+indexP-J_Increment);
							diffCoeff_s=*(diffusionCoeffdomain+indexP-J_Increment);
							
							if((*(Domain_Vartheta+indexP-J_Increment))>EPSILONU)
								v[2]=1.0/(*(Domain_Vartheta+indexP-J_Increment));
							else
								v[2]=0.0;
							
							rd[2]=1.0/(*(Domain_RhoDisperse+(indexP-J_Increment)*flagRhoDisperse));
						}
						else
						{
							Vol[2]=*(slabVolume[3]+i);
							diffCoeff_s=*(slabDiffusionCoeff[3]+i);
							
							if((*(slabVarthetaDisperse[3]+i))>EPSILONU)
								v[2]=1.0/(*(slabVarthetaDisperse[3]+i));
							else
								v[2]=0.0;
							
							rd[2]=1.0/(*(slabRhoDisperse[3]+i*flagRhoDisperse));
						}

						if(j!=(NJ-2))
						{
							Vol[3]=*(volumeDomain+indexP+J_Increment);
							diffCoeff_n=*(diffusionCoeffdomain+indexP+J_Increment);
							
							if((*(Domain_Vartheta+indexP+J_Increment))>EPSILONU)
								v[3]=1.0/(*(Domain_Vartheta+indexP+J_Increment));
							else
								v[3]=0.0;
							
							rd[3]=1.0/(*(Domain_RhoDisperse+(indexP+J_Increment)*flagRhoDisperse));	
						}
						else
						{
							Vol[3]=*(slabVolume[2]+i);
							diffCoeff_n=*(slabDiffusionCoeff[2]+i);
							
							if((*(slabVarthetaDisperse[2]+i))>EPSILONU)
								v[3]=1.0/(*(slabVarthetaDisperse[2]+i));
							else
								v[3]=0.0;
							
							rd[3]=1.0/(*(slabRhoDisperse[2]+i*flagRhoDisperse));
						}
						
						if(!(i==0 && j==0))
						{
							if(i==0)
							{
			
								if((*(slabVarthetaDisperse[1]+(NJ-1)*k+(j-1)))>EPSILONU)
									v[7]=1.0/(*(slabVarthetaDisperse[1]+(NJ-1)*k+(j-1)));
								else
									v[7]=0.0;
								
								rd[7]=1.0/(*(slabRhoDisperse[1]+((NJ-1)*k+(j-1))*flagRhoDisperse));	
							}
							else
							{
								if(j==0)
								{
				
									if((*(slabVarthetaDisperse[3]+(NI-1)*k+(i-1)))>EPSILONU)
										v[7]=1.0/(*(slabVarthetaDisperse[3]+(NI-1)*k+(i-1)));
									else
										v[7]=0.0;
									
									rd[7]=1.0/(*(slabRhoDisperse[3]+((NI-1)*k+(i-1))*flagRhoDisperse));	
								}
								else
								{
									if((*(Domain_Vartheta+indexP-J_Increment-1))>EPSILONU)
										v[7]=1.0/(*(Domain_Vartheta+indexP-J_Increment-1));
									else
										v[7]=0.0;
									
									rd[7]=1.0/(*(Domain_RhoDisperse+(indexP-J_Increment-1)*flagRhoDisperse));
								}
							}
						}
						else
						{
							v[7]=0.0;
							rd[7]=0.0;
						}
						if(!(i==0 && j==(NJ-2)))
						{
							if(i==0)
							{
			
								if((*(slabVarthetaDisperse[1]+(NJ-1)*k+(j+1)))>EPSILONU)
									v[6]=1.0/(*(slabVarthetaDisperse[1]+(NJ-1)*k+(j+1)));
								else
									v[6]=0.0;
								
								rd[6]=1.0/(*(slabRhoDisperse[1]+((NJ-1)*k+(j+1))*flagRhoDisperse));
							}
							else
							{
								if(j==(NJ-2))
								{
									if((*(slabVarthetaDisperse[2]+(NI-1)*k+(i-1)))>EPSILONU)
										v[6]=1.0/(*(slabVarthetaDisperse[2]+(NI-1)*k+(i-1)));
									else
										v[6]=0.0;
									
									rd[6]=1.0/(*(slabRhoDisperse[2]+((NI-1)*k+(i-1))*flagRhoDisperse));
								}
								else
								{
									if((*(Domain_Vartheta+indexP+J_Increment-1))>EPSILONU)
										v[6]=1.0/(*(Domain_Vartheta+indexP+J_Increment-1));
									else
										v[6]=0.0;
									
									rd[6]=1.0/(*(Domain_RhoDisperse+(indexP+J_Increment-1)*flagRhoDisperse));
								}
							}
						}
						else
						{
							v[6]=0.0;
							rd[6]=0.0;
						}
						if(!(i==0 && k==0))
						{
							if(i==0)
							{
								if((*(slabVarthetaDisperse[1]+(NJ-1)*(k-1)+j))>EPSILONU)
									v[9]=1.0/(*(slabVarthetaDisperse[1]+(NJ-1)*(k-1)+j));
								else
									v[9]=0.0;
								
								rd[9]=1.0/(*(slabRhoDisperse[1]+((NJ-1)*(k-1)+j)*flagRhoDisperse));
							}
							else
							{
								if(k==0)
								{
									if((*(slabVarthetaDisperse[5]+(NI-1)*j+(i-1)))>EPSILONU)
										v[9]=1.0/(*(slabVarthetaDisperse[5]+(NI-1)*j+(i-1)));
									else
										v[9]=0.0;
									
									rd[9]=1.0/(*(slabRhoDisperse[5]+((NI-1)*j+(i-1))*flagRhoDisperse));
								}
								else
								{
									if((*(Domain_Vartheta+indexP-K_Increment-1))>EPSILONU)
										v[9]=1.0/(*(Domain_Vartheta+indexP-K_Increment-1));
									else
										v[9]=0.0;
									
									rd[9]=1.0/(*(Domain_RhoDisperse+(indexP-K_Increment-1)*flagRhoDisperse));
								}
							}
		
						}
						else
						{
							v[9]=0.0;
							rd[9]=0.0;
						}
						if(!(i==0 && k==(NK-2)))
						{
							if(i==0)
							{
			
								if((*(slabVarthetaDisperse[1]+(NJ-1)*(k+1)+j))>EPSILONU)
									v[8]=1.0/(*(slabVarthetaDisperse[1]+(NJ-1)*(k+1)+j));
								else
									v[8]=0.0;
								
								rd[8]=1.0/(*(slabRhoDisperse[1]+((NJ-1)*(k+1)+j)*flagRhoDisperse));
							}
							else
							{
								if(k==(NK-2))
								{
									if((*(slabVarthetaDisperse[4]+(NI-1)*j+(i-1)))>EPSILONU)
										v[8]=1.0/(*(slabVarthetaDisperse[4]+(NI-1)*j+(i-1)));
									else
										v[8]=0.0;
									
									rd[8]=1.0/(*(slabRhoDisperse[4]+((NI-1)*j+(i-1))*flagRhoDisperse));
								}
								else
								{
									if((*(Domain_Vartheta+indexP+K_Increment-1))>EPSILONU)
										v[8]=1.0/(*(Domain_Vartheta+indexP+K_Increment-1));
									else
										v[8]=0.0;
									
									rd[8]=1.0/(*(Domain_RhoDisperse+(indexP+K_Increment-1)*flagRhoDisperse));
								}
							}
		
						}
						else
						{
							v[8]=0.0;
							rd[8]=0.0;
						}
						if(!(i==(NI-2) && j==0))
						{
							if(i==(NI-2))
							{
			
								if((*(slabVarthetaDisperse[0]+(NJ-1)*k+(j-1)))>EPSILONU)
									v[11]=1.0/(*(slabVarthetaDisperse[0]+(NJ-1)*k+(j-1)));
								else
									v[11]=0.0;
								
								rd[11]=1.0/(*(slabRhoDisperse[0]+((NJ-1)*k+(j-1))*flagRhoDisperse));
							}
							else
							{
								if(j==0)
								{
									if((*(slabVarthetaDisperse[3]+(NI-1)*k+(i+1)))>EPSILONU)
										v[11]=1.0/(*(slabVarthetaDisperse[3]+(NI-1)*k+(i+1)));
									else
										v[11]=0.0;
									
									rd[11]=1.0/(*(slabRhoDisperse[3]+((NI-1)*k+(i+1))*flagRhoDisperse));
								}
								else
								{
									if((*(Domain_Vartheta+indexP-J_Increment+1))>EPSILONU)
										v[11]=1.0/(*(Domain_Vartheta+indexP-J_Increment+1));
									else
										v[11]=0.0;
									
									rd[11]=1.0/(*(Domain_RhoDisperse+(indexP-J_Increment+1)*flagRhoDisperse));
								}
							}
						}
						else
						{
							v[11]=0.0;
							rd[11]=0.0;
						}
						if(!(i==(NI-2) && j==(NJ-2)))
						{
							if(i==(NI-2))
							{
			
								if((*(slabVarthetaDisperse[0]+(NJ-1)*k+(j+1)))>EPSILONU)
									v[10]=1.0/(*(slabVarthetaDisperse[0]+(NJ-1)*k+(j+1)));
								else
									v[10]=0.0;
								
								rd[10]=1.0/(*(slabRhoDisperse[0]+((NJ-1)*k+(j+1))*flagRhoDisperse));
							}
							else
							{
								if(j==(NJ-2))
								{
				
									if((*(slabVarthetaDisperse[2]+(NI-1)*k+(i+1)))>EPSILONU)
										v[10]=1.0/(*(slabVarthetaDisperse[2]+(NI-1)*k+(i+1)));
									else
										v[10]=0.0;
									
									rd[10]=1.0/(*(slabRhoDisperse[2]+((NI-1)*k+(i+1))*flagRhoDisperse));
								}
								else
								{
									if((*(Domain_Vartheta+indexP+J_Increment+1))>EPSILONU)
										v[10]=1.0/(*(Domain_Vartheta+indexP+J_Increment+1));
									else
										v[10]=0.0;
									
									rd[10]=1.0/(*(Domain_RhoDisperse+(indexP+J_Increment+1)*flagRhoDisperse));
								}
							}
						}
						else
						{
							v[10]=0.0;
							rd[10]=0.0;
						}
						if(!(i==(NI-2) && k==0))
						{
							if(i==(NI-2))
							{
			
								if((*(slabVarthetaDisperse[0]+(NJ-1)*(k-1)+j))>EPSILONU)
									v[13]=1.0/(*(slabVarthetaDisperse[0]+(NJ-1)*(k-1)+j));
								else
									v[13]=0.0;
								
								rd[13]=1.0/(*(slabRhoDisperse[0]+((NJ-1)*(k-1)+j)*flagRhoDisperse));
							}
							else
							{
								if(k==0)
								{
				
									if((*(slabVarthetaDisperse[5]+(NI-1)*j+(i+1)))>EPSILONU)
										v[13]=1.0/(*(slabVarthetaDisperse[5]+(NI-1)*j+(i+1)));
									else
										v[13]=0.0;
									
									rd[13]=1.0/(*(slabRhoDisperse[5]+((NI-1)*j+(i+1))*flagRhoDisperse));
								}
								else
								{
									if((*(Domain_Vartheta+indexP-K_Increment+1))>EPSILONU)
										v[13]=1.0/(*(Domain_Vartheta+indexP-K_Increment+1));
									else
										v[13]=0.0;
									
									rd[13]=1.0/(*(Domain_RhoDisperse+(indexP-K_Increment+1)*flagRhoDisperse));
								}
							}
						}
						else
						{
							v[13]=0.0;
							rd[13]=0.0;
						}
									
						if(!(i==(NI-2) && k==(NK-2)))
						{
							if(i==(NI-2))
							{
			
								if((*(slabVarthetaDisperse[0]+(NJ-1)*(k+1)+j))>EPSILONU)
									v[12]= 1.0/(*(slabVarthetaDisperse[0]+(NJ-1)*(k+1)+j));
								else
									v[12]=0.0;
								
								rd[12]= 1.0/(*(slabRhoDisperse[0]+((NJ-1)*(k+1)+j)*flagRhoDisperse));
							}
							else
							{
								if(k==(NK-2))
								{
				
									if((*(slabVarthetaDisperse[4]+(NI-1)*j+(i+1)))>EPSILONU)
										v[12]=1.0/(*(slabVarthetaDisperse[4]+(NI-1)*j+(i+1)));
									else
										v[12]=0.0;
									
									rd[12]=1.0/(*(slabRhoDisperse[4]+((NI-1)*j+(i+1))*flagRhoDisperse));
								}
								else
								{
									if((*(Domain_Vartheta+indexP+K_Increment+1))>EPSILONU)
										v[12]=1.0/(*(Domain_Vartheta+indexP+K_Increment+1));
									else
										v[12]=0.0;
									
									rd[12]=1.0/(*(Domain_RhoDisperse+(indexP+K_Increment+1)*flagRhoDisperse));
								}
							}
						}
						else
						{
							v[12]=0.0;
							rd[12]=0.0;
						}
							
						if(!(k==0 && j==0))
						{
							if(k==0)
							{
			
								if((*(slabVarthetaDisperse[5]+(NI-1)*(j-1)+i))>EPSILONU)
									v[17]=1.0/(*(slabVarthetaDisperse[5]+(NI-1)*(j-1)+i));
								else
									v[17]=0.0;
								
								rd[17]=1.0/(*(slabRhoDisperse[5]+((NI-1)*(j-1)+i)*flagRhoDisperse));
							}
							else
							{
								if(j==0)
								{
				
									if((*(slabVarthetaDisperse[3]+(NI-1)*(k-1)+i))>EPSILONU)
										v[17]=1.0/(*(slabVarthetaDisperse[3]+(NI-1)*(k-1)+i));
									else
										v[17]=0.0;
									
									rd[17]=1.0/(*(slabRhoDisperse[3]+((NI-1)*(k-1)+i)*flagRhoDisperse));
								}
								else
								{
									if((*(Domain_Vartheta+indexP-J_Increment-K_Increment))>EPSILONU)
										v[17]=1.0/(*(Domain_Vartheta+indexP-J_Increment-K_Increment));
									else
										v[17]=0.0;
									
									rd[17]=1.0/(*(Domain_RhoDisperse+(indexP-J_Increment-K_Increment)*flagRhoDisperse));
								}
							}
		
						}
						else
						{
							v[17]=0.0;
							rd[17]=0.0;
						}
						if(!(k==0 && j==(NJ-2)))
						{
							if(k==0 )
							{
			
								if((*(slabVarthetaDisperse[5]+(NI-1)*(j+1)+i))>EPSILONU)
									v[16]=1.0/(*(slabVarthetaDisperse[5]+(NI-1)*(j+1)+i));
								else
									v[16]=0.0;
								
								rd[16]=1.0/(*(slabRhoDisperse[5]+((NI-1)*(j+1)+i)*flagRhoDisperse));
							}
							else
							{
								if(j==(NJ-2))
								{
				
									if((*(slabVarthetaDisperse[2]+(NI-1)*(k-1)+i))>EPSILONU)
										v[16]=1.0/(*(slabVarthetaDisperse[2]+(NI-1)*(k-1)+i));
									else
										v[16]=0.0;
									
									rd[16]=1.0/(*(slabRhoDisperse[2]+((NI-1)*(k-1)+i)*flagRhoDisperse));
								}
								else
								{
									if((*(Domain_Vartheta+indexP+J_Increment-K_Increment))>EPSILONU)
										v[16]=1.0/(*(Domain_Vartheta+indexP+J_Increment-K_Increment));
									else
										v[16]=0.0;
									
									rd[16]=1.0/(*(Domain_RhoDisperse+(indexP+J_Increment-K_Increment)*flagRhoDisperse));
								}
							}
		
						}
						else
						{
							v[16]=0.0;
							rd[16]=0.0;
						}
						if(!(k==(NK-2) && j==0))
						{
							if(k==(NK-2))
							{
			
								if((*(slabVarthetaDisperse[4]+(NI-1)*(j-1)+i))>EPSILONU)
									v[15]=1.0/(*(slabVarthetaDisperse[4]+(NI-1)*(j-1)+i));
								else
									v[15]=0.0;
								
								rd[15]=1.0/(*(slabRhoDisperse[4]+((NI-1)*(j-1)+i)*flagRhoDisperse));
							}
							else
							{
								if(j==0)
								{
				
									if((*(slabVarthetaDisperse[3]+(NI-1)*(k+1)+i))>EPSILONU)
										v[15]=1.0/(*(slabVarthetaDisperse[3]+(NI-1)*(k+1)+i));
									else
										v[15]=0.0;
									
									rd[15]=1.0/(*(slabRhoDisperse[3]+((NI-1)*(k+1)+i)*flagRhoDisperse));
								}
								else
								{
									if((*(Domain_Vartheta+indexP-J_Increment+K_Increment))>EPSILONU)
										v[15]=1.0/(*(Domain_Vartheta+indexP-J_Increment+K_Increment));
									else
										v[15]=0.0;
									
									rd[15]=1.0/(*(Domain_RhoDisperse+(indexP-J_Increment+K_Increment)*flagRhoDisperse));
								}
							}
						}
						else
						{
							v[15]=0.0;
							rd[15]=0.0;
						}
						if(!(k==(NK-2) && j==(NJ-2)))
						{
							if(k==(NK-2))
							{
			
								if((*(slabVarthetaDisperse[4]+(NI-1)*(j+1)+i))>EPSILONU)
									v[14]=1.0/(*(slabVarthetaDisperse[4]+(NI-1)*(j+1)+i));
								else
									v[14]=0.0;
								
								rd[14]=1.0/(*(slabRhoDisperse[4]+((NI-1)*(j+1)+i)*flagRhoDisperse));
							}
							else
							{
								if(j==(NJ-2))
								{
				
									if((*(slabVarthetaDisperse[2]+(NI-1)*(k+1)+i))>EPSILONU)
										v[14]=1.0/(*(slabVarthetaDisperse[2]+(NI-1)*(k+1)+i));
									else
										v[14]=0.0;
									
									rd[14]=1.0/(*(slabRhoDisperse[2]+((NI-1)*(k+1)+i)*flagRhoDisperse));
								}
								else
								{
									if((*(Domain_Vartheta+indexP+J_Increment+K_Increment))>EPSILONU)
										v[14]=1.0/(*(Domain_Vartheta+indexP+J_Increment+K_Increment));
									else
										v[14]=0.0;
									
									rd[14]=1.0/(*(Domain_RhoDisperse+(indexP+J_Increment+K_Increment)*flagRhoDisperse));
								}
							}
						}
						else
						{
							v[14]=0.0;
							rd[14]=0.0;		
						}
						diffCoeff_p=*(diffusionCoeffdomain+indexP);
						
						vol_p=*(volumeDomain+indexP);
// 						for(che=0;che<18;che++)
// 							printf("\nV[%d]=%e",che,rd[che]);
						
						
						if(*(Domain_Vartheta+indexP)>EPSILONU)
						v_p=1.0/(*(Domain_Vartheta+indexP));
						else
						v_p=0.0;
								
						rd_p=1.0/(*(Domain_RhoDisperse+indexP*flagRhoDisperse));
						
						diffCoeff_face[0] = (diffCoeff_w*vol_p + diffCoeff_p*Vol[0])/(vol_p+Vol[0]);
						diffCoeff_face[1] = (diffCoeff_e*vol_p + diffCoeff_p*Vol[1])/(vol_p+Vol[1]);
						diffCoeff_face[2] = (diffCoeff_s*vol_p + diffCoeff_p*Vol[2])/(vol_p+Vol[2]);
						diffCoeff_face[3] = (diffCoeff_n*vol_p + diffCoeff_p*Vol[3])/(vol_p+Vol[3]);
				
						*(diffdia  +indexP) = 0.0;
						*(CoeffN[0]+indexP) = 0.0;
						*(CoeffN[1]+indexP) = 0.0;
						*(CoeffN[2]+indexP) = 0.0;
						*(CoeffN[3]+indexP) = 0.0;
					
						for (face=0;face<4;face++)
						{
							*(diffdia  +indexP) = *(diffdia  +indexP) + ((*(Diag[face]+indexP)) * diffCoeff_face[face])*v_p*rd_p;
							*(CoeffN[0]+indexP) = *(CoeffN[0]+indexP) + ((*(D[face]+indexP))    * diffCoeff_face[face])*v[0]*rd[0];
							*(CoeffN[1]+indexP) = *(CoeffN[1]+indexP) + ((*(D[face+6]+indexP))    * diffCoeff_face[face])*v[1]*rd[1];
							*(CoeffN[2]+indexP) = *(CoeffN[2]+indexP) +( (*(D[face+12]+indexP))    * diffCoeff_face[face])*v[2]*rd[2];
						
							*(CoeffN[3]+indexP) = *(CoeffN[3]+indexP) + ((*(D[face+18]+indexP))    * diffCoeff_face[face])*v[3]*rd[3];
						}
						*(CoeffN[6]+indexP) = ((*(DN[0]+indexP))    * diffCoeff_face[0] + (*(DN[1]+indexP))    * diffCoeff_face[3])*v[6]*rd[6];
						*(CoeffN[7]+indexP) = ((*(DN[2]+indexP))    * diffCoeff_face[0] + (*(DN[3]+indexP))    * diffCoeff_face[2])*v[7]*rd[7];
						*(CoeffN[10]+indexP) =((*(DN[8]+indexP))    * diffCoeff_face[1] + (*(DN[9]+indexP))    * diffCoeff_face[3])*v[10]*rd[10];
						*(CoeffN[11]+indexP) =((*(DN[10]+indexP))   * diffCoeff_face[1] +(*(DN[11]+indexP))    * diffCoeff_face[2])*v[11]*rd[11];
					
						indexP++;
					}
				}
				else if(base->BaseParameters->flag2d == 1)
				{
					indexP=0;
					for(k=0;k<NK-1;k++)
						for(j=0;j<NJ-1;j++)
							for(i=0;i<NI-1;i++)
					{
						if(i!=0)
						{
							Vol[0]=*(volumeDomain+indexP-1);
							diffCoeff_w=*(diffusionCoeffdomain+indexP-1);
						}
						else
						{
							Vol[0] =*(slabVolume[1]+j);
							diffCoeff_w =*(slabDiffusionCoeff[1]+j);
						}
						if(i!=(NI-2))
						{
							Vol[1]=*(volumeDomain+indexP+1);
							diffCoeff_e=*(diffusionCoeffdomain+indexP+1);
						}
						else
						{
							Vol[1]= *(slabVolume[0]+j);
							diffCoeff_e=*(slabDiffusionCoeff[0]+j);
						}
	
						if(j!=0)
						{
							Vol[2]=*(volumeDomain+indexP-J_Increment);
							diffCoeff_s=*(diffusionCoeffdomain+indexP-J_Increment);
						}
						else
						{
							Vol[2]=*(slabVolume[3]+i);
							diffCoeff_s=*(slabDiffusionCoeff[3]+i);
						}
	
						if(j!=(NJ-2))
						{
							Vol[3]=*(volumeDomain+indexP+J_Increment);
							diffCoeff_n=*(diffusionCoeffdomain+indexP+J_Increment);
						}
						else
						{
							Vol[3]=*(slabVolume[2]+i);
							diffCoeff_n=*(slabDiffusionCoeff[2]+i);
						}
						diffCoeff_p=*(diffusionCoeffdomain+indexP);
						vol_p=*(volumeDomain+indexP);
						
						
						diffCoeff_face[0] = (diffCoeff_w*vol_p + diffCoeff_p*Vol[0])/(vol_p+Vol[0]);
						diffCoeff_face[1] = (diffCoeff_e*vol_p + diffCoeff_p*Vol[1])/(vol_p+Vol[1]);
						diffCoeff_face[2] = (diffCoeff_s*vol_p + diffCoeff_p*Vol[2])/(vol_p+Vol[2]);
						diffCoeff_face[3] = (diffCoeff_n*vol_p + diffCoeff_p*Vol[3])/(vol_p+Vol[3]);
						
						*(diffdia  +indexP) = 0.0;
						*(CoeffN[0]+indexP) = 0.0;
						*(CoeffN[1]+indexP) = 0.0;
						*(CoeffN[2]+indexP) = 0.0;
						*(CoeffN[3]+indexP) = 0.0;
						
						for (face=0;face<4;face++)
						{
							*(diffdia  +indexP) = *(diffdia  +indexP) + (*(Diag[face]+indexP)) * diffCoeff_face[face];
							*(CoeffN[0]+indexP) = *(CoeffN[0]+indexP) + (*(D[face]+indexP))    * diffCoeff_face[face];
							*(CoeffN[1]+indexP) = *(CoeffN[1]+indexP) + (*(D[face+6]+indexP))    * diffCoeff_face[face];
							*(CoeffN[2]+indexP) = *(CoeffN[2]+indexP) + (*(D[face+12]+indexP))    * diffCoeff_face[face];
							
							*(CoeffN[3]+indexP) = *(CoeffN[3]+indexP) + (*(D[face+18]+indexP))    * diffCoeff_face[face];
						}
						*(CoeffN[6]+indexP) = (*(DN[0]+indexP))    * diffCoeff_face[0] + (*(DN[1]+indexP))    * diffCoeff_face[3];
						*(CoeffN[7]+indexP) = (*(DN[2]+indexP))    * diffCoeff_face[0] + (*(DN[3]+indexP))    * diffCoeff_face[2];
						*(CoeffN[10]+indexP) =(*(DN[8]+indexP))    * diffCoeff_face[1] + (*(DN[9]+indexP))    * diffCoeff_face[3];
						*(CoeffN[11]+indexP) =(*(DN[10]+indexP))   * diffCoeff_face[1] +(*(DN[11]+indexP))    * diffCoeff_face[2];
	
						indexP++;
					}
		
				}
				else if((base->BaseParameters->flag2d == 0)&&((vartag==62)||(vartag==63)||(vartag==64)))
				{
					indexP=0;
					for(k=0;k<NK-1;k++)
					for(j=0;j<NJ-1;j++)
					for(i=0;i<NI-1;i++)
					{
						if(i!=0)
						{
							Vol[0]=*(volumeDomain+indexP-1);
							diffCoeff_w=*(diffusionCoeffdomain+indexP-1);
							v[0] = 1/(*(Domain_Vartheta+indexP-1));
							rd[0] = 1/(*(Domain_RhoDisperse+(indexP-1)*flagRhoDisperse));
						}
						else
						{
							Vol[0] =*(slabVolume[1]+j);
							diffCoeff_w =*(slabDiffusionCoeff[1]+j);
							v[0] =1/(*(slabVarthetaDisperse[1]+j));
							rd[0] =1/(*(slabRhoDisperse[1]+j*flagRhoDisperse));
						}
						
						if(i!=(NI-2))
						{
							Vol[1]=*(volumeDomain+indexP+1);
							diffCoeff_e=*(diffusionCoeffdomain+indexP+1);
							v[1]=1/(*(Domain_Vartheta+indexP+1));
							rd[1] = 1/(*(Domain_RhoDisperse+(indexP+1)*flagRhoDisperse));
						}
						else
						{
							Vol[1]= *(slabVolume[0]+j);
							diffCoeff_e=*(slabDiffusionCoeff[0]+j);
							v[1]=1/(*(slabVarthetaDisperse[0]+j));
							rd[1]=1/(*(slabRhoDisperse[0]+j*flagRhoDisperse));
						}
	
						if(j!=0)
						{
							Vol[2]=*(volumeDomain+indexP-J_Increment);
							diffCoeff_s=*(diffusionCoeffdomain+indexP-J_Increment);
							v[2]=1/(*(Domain_Vartheta+indexP-J_Increment));
							rd[2]=1/(*(Domain_RhoDisperse+(indexP-J_Increment)*flagRhoDisperse));
						}
						else
						{
							Vol[2]=*(slabVolume[3]+i);
							diffCoeff_s=*(slabDiffusionCoeff[3]+i);
							v[2]=1/(*(slabVarthetaDisperse[3]+i));
							rd[2]=1/(*(slabRhoDisperse[3]+i*flagRhoDisperse));
						}
	
						if(j!=(NJ-2))
						{
							Vol[3]=*(volumeDomain+indexP+J_Increment);
							diffCoeff_n=*(diffusionCoeffdomain+indexP+J_Increment);
							v[3]=1/(*(Domain_Vartheta+indexP+J_Increment));
							rd[3]=1/(*(Domain_RhoDisperse+(indexP+J_Increment)*flagRhoDisperse));	
						}
						else
						{
							Vol[3]=*(slabVolume[2]+i);
							diffCoeff_n=*(slabDiffusionCoeff[2]+i);
							v[3]=1/(*(slabVarthetaDisperse[2]+i));
							rd[3]=1/(*(slabRhoDisperse[2]+i*flagRhoDisperse));
						}
						if(k!=0)
						{
							Vol[4]=*(volumeDomain+indexP-K_Increment);
							diffCoeff_b=*(diffusionCoeffdomain+indexP-K_Increment);
							v[4]=1/(*(Domain_Vartheta+indexP-K_Increment));
							rd[4]=1/(*(Domain_RhoDisperse+(indexP-K_Increment)*flagRhoDisperse));
						}
						else
						{
							Vol[4]=*(slabVolume[5]+(NI-1)*j+i);
							diffCoeff_b=*(slabDiffusionCoeff[5]+(NI-1)*j+i);
							v[4]=1/(*(slabVarthetaDisperse[5]+(NI-1)*j+i));
							rd[4]=1/(*(slabRhoDisperse[5]+((NI-1)*j+i)*flagRhoDisperse));
						}
						
						if(k!=(NK-2))
						{
							Vol[5] = *(volumeDomain+indexP+K_Increment);
							diffCoeff_t=*(diffusionCoeffdomain+indexP+K_Increment);
							v[5] = 1/(*(Domain_Vartheta+indexP+K_Increment));
							rd[5] = 1/(*(Domain_RhoDisperse+(indexP+K_Increment)*flagRhoDisperse));
						}
						else
						{
							Vol[5] = *(slabVolume[4]+(NI-1)*j+i);
							diffCoeff_t=*(slabDiffusionCoeff[4]+(NI-1)*j+i);
							v[5] = 1/(*(slabVarthetaDisperse[4]+(NI-1)*j+i));
							rd[5] = 1/(*(slabRhoDisperse[4]+((NI-1)*j+i)*flagRhoDisperse));
						}
						if(!(i==0 && j==0))
						{
							if(i==0)
							{
			
								v[7]=1/(*(slabVarthetaDisperse[1]+(NJ-1)*k+(j-1)));
								rd[7]=1/(*(slabRhoDisperse[1]+((NJ-1)*k+(j-1))*flagRhoDisperse));	
							}
							else
							{
								if(j==0)
								{
				
									v[7]=1/(*(slabVarthetaDisperse[3]+(NI-1)*k+(i-1)));
									rd[7]=1/(*(slabRhoDisperse[3]+((NI-1)*k+(i-1))*flagRhoDisperse));	
								}
								else
								{
									v[7]=1/(*(Domain_Vartheta+indexP-J_Increment-1));
									rd[7]=1/(*(Domain_RhoDisperse+(indexP-J_Increment-1)*flagRhoDisperse));
								}
							}
						}
						else
						{
							v[7]=0.0;
							rd[7]=0.0;
						}
						if(!(i==0 && j==(NJ-2)))
						{
							if(i==0)
							{
			
								v[6]=1/(*(slabVarthetaDisperse[1]+(NJ-1)*k+(j+1)));
								rd[6]=1/(*(slabRhoDisperse[1]+((NJ-1)*k+(j+1))*flagRhoDisperse));
							}
							else
							{
								if(j==(NJ-2))
								{
				
									v[6]=1/(*(slabVarthetaDisperse[2]+(NI-1)*k+(i-1)));
									rd[6]=1/(*(slabRhoDisperse[2]+((NI-1)*k+(i-1))*flagRhoDisperse));
								}
								else
								{
									v[6]=1/(*(Domain_Vartheta+indexP+J_Increment-1));
									rd[6]=1/(*(Domain_RhoDisperse+(indexP+J_Increment-1)*flagRhoDisperse));
								}
							}
						}
						else
						{
							v[6]=0.0;
							rd[6]=0.0;
						}
						if(!(i==0 && k==0))
						{
							if(i==0)
							{
		
								v[9]=1/(*(slabVarthetaDisperse[1]+(NJ-1)*(k-1)+j));
								rd[9]=1/(*(slabRhoDisperse[1]+((NJ-1)*(k-1)+j)*flagRhoDisperse));
							}
							else
							{
								if(k==0)
								{
			
									v[9]=1/(*(slabVarthetaDisperse[5]+(NI-1)*j+(i-1)));
									rd[9]=1/(*(slabRhoDisperse[5]+((NI-1)*j+(i-1))*flagRhoDisperse));
								}
								else
								{
									v[9]=1/(*(Domain_Vartheta+indexP-K_Increment-1));
									rd[9]=1/(*(Domain_RhoDisperse+(indexP-K_Increment-1)*flagRhoDisperse));
								}
							}
	
						}
						else
						{
							v[9]=0.0;
							rd[9]=0.0;
						}
						if(!(i==0 && k==(NK-2)))
						{
							if(i==0)
							{
		
								v[8]=1/(*(slabVarthetaDisperse[1]+(NJ-1)*(k+1)+j));
								rd[8]=1/(*(slabRhoDisperse[1]+((NJ-1)*(k+1)+j)*flagRhoDisperse));
							}
							else
							{
								if(k==(NK-2))
								{
			
									v[8]=1/(*(slabVarthetaDisperse[4]+(NI-1)*j+(i-1)));
									rd[8]=1/(*(slabRhoDisperse[4]+((NI-1)*j+(i-1))*flagRhoDisperse));
								}
								else
								{
									v[8]=1/(*(Domain_Vartheta+indexP+K_Increment-1));
									rd[8]=1/(*(Domain_RhoDisperse+(indexP+K_Increment-1)*flagRhoDisperse));
								}
							}
	
						}
						else
						{
							v[8]=0.0;
							rd[8]=0.0;
						}
						if(!(i==(NI-2) && j==0))
						{
							if(i==(NI-2))
							{
			
								v[11]=1/(*(slabVarthetaDisperse[0]+(NJ-1)*k+(j-1)));
								rd[11]=1/(*(slabRhoDisperse[0]+((NJ-1)*k+(j-1))*flagRhoDisperse));
							}
							else
							{
								if(j==0)
								{
				
									v[11]=1/(*(slabVarthetaDisperse[3]+(NI-1)*k+(i+1)));
									rd[11]=1/(*(slabRhoDisperse[3]+((NI-1)*k+(i+1))*flagRhoDisperse));
								}
								else
								{
									v[11]=1/(*(Domain_Vartheta+indexP-J_Increment+1));
									rd[11]=1/(*(Domain_RhoDisperse+(indexP-J_Increment+1)*flagRhoDisperse));
								}
							}
						}
						else
						{
							v[11]=0.0;
							rd[11]=0.0;
						}
						if(!(i==(NI-2) && j==(NJ-2)))
						{
							if(i==(NI-2))
							{
			
								v[10]=1/(*(slabVarthetaDisperse[0]+(NJ-1)*k+(j+1)));
								rd[10]=1/(*(slabRhoDisperse[0]+((NJ-1)*k+(j+1))*flagRhoDisperse));
							}
							else
							{
								if(j==(NJ-2))
								{
				
									v[10]=1/(*(slabVarthetaDisperse[2]+(NI-1)*k+(i+1)));
									rd[10]=1/(*(slabRhoDisperse[2]+((NI-1)*k+(i+1))*flagRhoDisperse));
								}
								else
								{
									v[10]=1/(*(Domain_Vartheta+indexP+J_Increment+1));
									rd[10]=1/(*(Domain_RhoDisperse+(indexP+J_Increment+1)*flagRhoDisperse));
								}
							}
						}
						else
						{
							v[10]=0.0;
							rd[10]=0.0;
						}
						if(!(i==(NI-2) && k==0))
						{
							if(i==(NI-2))
							{
			
								v[13]=1/(*(slabVarthetaDisperse[0]+(NJ-1)*(k-1)+j));
								rd[13]=1/(*(slabRhoDisperse[0]+((NJ-1)*(k-1)+j)*flagRhoDisperse));
							}
							else
							{
								if(k==0)
								{
				
									v[13]=1/(*(slabVarthetaDisperse[5]+(NI-1)*j+(i+1)));
									rd[13]=1/(*(slabRhoDisperse[5]+((NI-1)*j+(i+1))*flagRhoDisperse));
								}
								else
								{
									v[13]=1/(*(Domain_Vartheta+indexP-K_Increment+1));
									rd[13]=1/(*(Domain_RhoDisperse+(indexP-K_Increment+1)*flagRhoDisperse));
								}
							}
						}
						else
						{
							v[13]=0.0;
							rd[13]=0.0;
						}
						if(!(i==(NI-2) && k==(NK-2)))
						{
							if(i==(NI-2))
							{
			
								v[12]= 1/(*(slabVarthetaDisperse[0]+(NJ-1)*(k+1)+j));
								rd[12]= 1/(*(slabRhoDisperse[0]+((NJ-1)*(k+1)+j)*flagRhoDisperse));
							}
							else
							{
								if(k==(NK-2))
								{
				
									v[12]=1/(*(slabVarthetaDisperse[4]+(NI-1)*j+(i+1)));
									rd[12]=1/(*(slabRhoDisperse[4]+((NI-1)*j+(i+1))*flagRhoDisperse));
								}
								else
								{
									v[12]=1/(*(Domain_Vartheta+indexP+K_Increment+1));
									rd[12]=1/(*(Domain_RhoDisperse+(indexP+K_Increment+1)*flagRhoDisperse));
								}
							}
						}
						else
						{
							v[12]=0.0;
							rd[12]=0.0;
						}	
						if(!(k==0 && j==0))
						{
							if(k==0)
							{
			
								v[17]=1/(*(slabVarthetaDisperse[5]+(NI-1)*(j-1)+i));
								rd[17]=1/(*(slabRhoDisperse[5]+((NI-1)*(j-1)+i)*flagRhoDisperse));
							}
							else
							{
								if(j==0)
								{
				
									v[17]=1/(*(slabVarthetaDisperse[3]+(NI-1)*(k-1)+i));
									rd[17]=1/(*(slabRhoDisperse[3]+((NI-1)*(k-1)+i)*flagRhoDisperse));
								}
								else
								{
									v[17]=1/(*(Domain_Vartheta+indexP-J_Increment-K_Increment));
									rd[17]=1/(*(Domain_RhoDisperse+(indexP-J_Increment-K_Increment)*flagRhoDisperse));
								}
							}
		
						}
						else
						{
							v[17]=0.0;
							rd[17]=0.0;
						}
						if(!(k==0 && j==(NJ-2)))
						{
							if(k==0 )
							{
			
								v[16]=1/(*(slabVarthetaDisperse[5]+(NI-1)*(j+1)+i));
								rd[16]=1/(*(slabRhoDisperse[5]+((NI-1)*(j+1)+i)*flagRhoDisperse));
							}
							else
							{
								if(j==(NJ-2))
								{
				
									v[16]=1/(*(slabVarthetaDisperse[2]+(NI-1)*(k-1)+i));
									rd[16]=1/(*(slabRhoDisperse[2]+((NI-1)*(k-1)+i)*flagRhoDisperse));
								}
								else
								{
									v[16]=1/(*(Domain_Vartheta+indexP+J_Increment-K_Increment));
									rd[16]=1/(*(Domain_RhoDisperse+(indexP+J_Increment-K_Increment)*flagRhoDisperse));
								}
							}
		
						}
						else
						{
							v[16]=0.0;
							rd[16]=0.0;
						}
						if(!(k==(NK-2) && j==0))
						{
							if(k==(NK-2))
							{
			
								v[15]=1/(*(slabVarthetaDisperse[4]+(NI-1)*(j-1)+i));
								rd[15]=1/(*(slabRhoDisperse[4]+((NI-1)*(j-1)+i)*flagRhoDisperse));
							}
							else
							{
								if(j==0)
								{
				
									v[15]=1/(*(slabVarthetaDisperse[3]+(NI-1)*(k+1)+i));
									rd[15]=1/(*(slabRhoDisperse[3]+((NI-1)*(k+1)+i)*flagRhoDisperse));
								}
								else
								{
									v[15]=1/(*(Domain_Vartheta+indexP-J_Increment+K_Increment));
									rd[15]=1/(*(Domain_RhoDisperse+(indexP-J_Increment+K_Increment)*flagRhoDisperse));
								}
							}
						}
						else
						{
							v[15]=0.0;
							rd[15]=0.0;
						}
						if(!(k==(NK-2) && j==(NJ-2)))
						{
							if(k==(NK-2))
							{
			
								v[14]=1/(*(slabVarthetaDisperse[4]+(NI-1)*(j+1)+i));
								rd[14]=1/(*(slabRhoDisperse[4]+((NI-1)*(j+1)+i)*flagRhoDisperse));
							}
							else
							{
								if(j==(NJ-2))
								{
				
									v[14]=1/(*(slabVarthetaDisperse[2]+(NI-1)*(k+1)+i));
									rd[14]=1/(*(slabRhoDisperse[2]+((NI-1)*(k+1)+i)*flagRhoDisperse));
								}
								else
								{
									v[14]=1/(*(Domain_Vartheta+indexP+J_Increment+K_Increment));
									rd[14]=1/(*(Domain_RhoDisperse+(indexP+J_Increment+K_Increment)*flagRhoDisperse));
								}
							}
						}
						else
						{
							v[14]=0.0;
							rd[14]=0.0;		
						}
						diffCoeff_p=*(diffusionCoeffdomain+indexP);
						vol_p=*(volumeDomain+indexP);
						v_p=1/(*(Domain_Vartheta+indexP));
						rd_p=1/(*(Domain_RhoDisperse+indexP*flagRhoDisperse));
						
						diffCoeff_face[0] = (diffCoeff_w*vol_p + diffCoeff_p*Vol[0])/(vol_p+Vol[0]);
						diffCoeff_face[1] = (diffCoeff_e*vol_p + diffCoeff_p*Vol[1])/(vol_p+Vol[1]);
						diffCoeff_face[2] = (diffCoeff_s*vol_p + diffCoeff_p*Vol[2])/(vol_p+Vol[2]);
						diffCoeff_face[3] = (diffCoeff_n*vol_p + diffCoeff_p*Vol[3])/(vol_p+Vol[3]);
						diffCoeff_face[4] = (diffCoeff_b*vol_p + diffCoeff_p*Vol[4])/(vol_p+Vol[4]);
						diffCoeff_face[5] = (diffCoeff_t*vol_p + diffCoeff_p*Vol[5])/(vol_p+Vol[5]);
						
						*(diffdia  +indexP) = 0.0;
						*(CoeffN[0]+indexP) = 0.0;
						*(CoeffN[1]+indexP) = 0.0;
						*(CoeffN[2]+indexP) = 0.0;
						*(CoeffN[3]+indexP) = 0.0;
						*(CoeffN[4]+indexP) = 0.0;
						*(CoeffN[5]+indexP) = 0.0;
						
						for (face=0;face<6;face++)
						{
							*(diffdia  +indexP) = *(diffdia  +indexP) + ((*(Diag[face]+indexP)) * diffCoeff_face[face])*v_p*rd_p;
							*(CoeffN[0]+indexP) = *(CoeffN[0]+indexP) + ((*(D[face]+indexP))    * diffCoeff_face[face])*v[0]*rd[0];
							*(CoeffN[1]+indexP) = *(CoeffN[1]+indexP) + ((*(D[face+6]+indexP))    * diffCoeff_face[face])*v[1]*rd[1];
							*(CoeffN[2]+indexP) = *(CoeffN[2]+indexP) + ((*(D[face+12]+indexP))    * diffCoeff_face[face])*v[2]*rd[2];
							*(CoeffN[3]+indexP) = *(CoeffN[3]+indexP) + ((*(D[face+18]+indexP))    * diffCoeff_face[face])*v[3]*rd[3];
							*(CoeffN[4]+indexP) = *(CoeffN[4]+indexP) + ((*(D[face+24]+indexP))    * diffCoeff_face[face])*v[4]*rd[4];
							*(CoeffN[5]+indexP) = *(CoeffN[5]+indexP) + ((*(D[face+30]+indexP))    * diffCoeff_face[face])*v[5]*rd[5];
						}
	
						*(CoeffN[6]+indexP) = ((*(DN[0]+indexP))    * diffCoeff_face[0] + (*(DN[1]+indexP))    * diffCoeff_face[3])*v[6]*rd[6];
						*(CoeffN[7]+indexP) = ((*(DN[2]+indexP))    * diffCoeff_face[0] + (*(DN[3]+indexP))    * diffCoeff_face[2])*v[7]*rd[7];
						
						*(CoeffN[8]+indexP) = ((*(DN[4]+indexP))    * diffCoeff_face[0] + (*(DN[5]+indexP))    * diffCoeff_face[5])*v[8]*rd[8];
						*(CoeffN[9]+indexP) = ((*(DN[6]+indexP))    * diffCoeff_face[0] + (*(DN[7]+indexP))    * diffCoeff_face[4])*v[9]*rd[9];
						
						*(CoeffN[10]+indexP) =((*(DN[8]+indexP))    * diffCoeff_face[1] + (*(DN[9]+indexP))    * diffCoeff_face[3])*v[10]*rd[10];
						*(CoeffN[11]+indexP) =((*(DN[10]+indexP))   * diffCoeff_face[1] +(*(DN[11]+indexP))    * diffCoeff_face[2])*v[11]*rd[11];
						
						*(CoeffN[12]+indexP) = ((*(DN[12]+indexP))    * diffCoeff_face[1] + (*(DN[13]+indexP))    * diffCoeff_face[5])*v[12]*rd[12];
						*(CoeffN[13]+indexP) = ((*(DN[14]+indexP))    * diffCoeff_face[1] + (*(DN[15]+indexP))    * diffCoeff_face[4])*v[13]*rd[13];
						
						*(CoeffN[14]+indexP) = ((*(DN[16]+indexP))    * diffCoeff_face[5] + (*(DN[17]+indexP))    * diffCoeff_face[3])*v[14]*rd[14];
						*(CoeffN[15]+indexP) = ((*(DN[18]+indexP))    * diffCoeff_face[5] + (*(DN[19]+indexP))    * diffCoeff_face[2])*v[15]*rd[15];
						
						*(CoeffN[16]+indexP) = ((*(DN[20]+indexP))    * diffCoeff_face[4] + (*(DN[21]+indexP))    * diffCoeff_face[3])*v[16]*rd[16];
						*(CoeffN[17]+indexP) = ((*(DN[22]+indexP))    * diffCoeff_face[4] + (*(DN[23]+indexP))    * diffCoeff_face[2])*v[17]*rd[17];
						
						indexP++;
					}
				}
				else if(base->BaseParameters->flag2d == 0)
				{
					indexP=0;
					for(k=0;k<NK-1;k++)
						for(j=0;j<NJ-1;j++)
							for(i=0;i<NI-1;i++)
					{
						if(i!=0)
						{
							Vol[0]=*(volumeDomain+indexP-1);
							diffCoeff_w=*(diffusionCoeffdomain+indexP-1);
						}
						else
						{
							Vol[0] =*(slabVolume[1]+(NJ-1)*k+j);
							diffCoeff_w =*(slabDiffusionCoeff[1]+(NJ-1)*k+j);
						}
						
						if(i!=(NI-2))
						{
							Vol[1]=*(volumeDomain+indexP+1);
							diffCoeff_e=*(diffusionCoeffdomain+indexP+1);
						}
						else
						{
							Vol[1]= *(slabVolume[0]+(NJ-1)*k+j);
							diffCoeff_e=*(slabDiffusionCoeff[0]+(NJ-1)*k+j);
						}
	
						if(j!=0)
						{
							Vol[2]=*(volumeDomain+indexP-J_Increment);
							diffCoeff_s=*(diffusionCoeffdomain+indexP-J_Increment);
						}
						else
						{
							Vol[2]=*(slabVolume[3]+i);
							diffCoeff_s=*(slabDiffusionCoeff[3]+(NI-1)*k+i);
						}
	
						if(j!=(NJ-2))
						{
							Vol[3]=*(volumeDomain+indexP+J_Increment);
							diffCoeff_n=*(diffusionCoeffdomain+indexP+J_Increment);
						}
						else
						{
							Vol[3]=*(slabVolume[2]+i);
							diffCoeff_n=*(slabDiffusionCoeff[2]+(NI-1)*k+i);
						}
						if(k!=0)
						{
							Vol[4]=*(volumeDomain+indexP-K_Increment);
							diffCoeff_b=*(diffusionCoeffdomain+indexP-K_Increment);
						}
						else
						{
							Vol[4]=*(slabVolume[5]+(NI-1)*j+i);
							diffCoeff_b=*(slabDiffusionCoeff[5]+(NI-1)*j+i);
						}
						
						if(k!=(NK-2))
						{
							Vol[5] = *(volumeDomain+indexP+K_Increment);
							diffCoeff_t=*(diffusionCoeffdomain+indexP+K_Increment);
						}
						else
						{
							Vol[5] = *(slabVolume[4]+(NI-1)*j+i);
							diffCoeff_t=*(slabDiffusionCoeff[4]+(NI-1)*j+i);
						}
						
						diffCoeff_p=*(diffusionCoeffdomain+indexP);
						vol_p=*(volumeDomain+indexP);
						
						diffCoeff_face[0] = (diffCoeff_w*vol_p + diffCoeff_p*Vol[0])/(vol_p+Vol[0]);
						diffCoeff_face[1] = (diffCoeff_e*vol_p + diffCoeff_p*Vol[1])/(vol_p+Vol[1]);
						diffCoeff_face[2] = (diffCoeff_s*vol_p + diffCoeff_p*Vol[2])/(vol_p+Vol[2]);
						diffCoeff_face[3] = (diffCoeff_n*vol_p + diffCoeff_p*Vol[3])/(vol_p+Vol[3]);
						diffCoeff_face[4] = (diffCoeff_b*vol_p + diffCoeff_p*Vol[4])/(vol_p+Vol[4]);
						diffCoeff_face[5] = (diffCoeff_t*vol_p + diffCoeff_p*Vol[5])/(vol_p+Vol[5]);
						
						*(diffdia  +indexP) = 0.0;
						*(CoeffN[0]+indexP) = 0.0;
						*(CoeffN[1]+indexP) = 0.0;
						*(CoeffN[2]+indexP) = 0.0;
						*(CoeffN[3]+indexP) = 0.0;
						*(CoeffN[4]+indexP) = 0.0;
						*(CoeffN[5]+indexP) = 0.0;
						
						for (face=0;face<6;face++)
						{
							*(diffdia  +indexP) = *(diffdia  +indexP) + (*(Diag[face]+indexP)) * diffCoeff_face[face];
							*(CoeffN[0]+indexP) = *(CoeffN[0]+indexP) + (*(D[face]+indexP))    * diffCoeff_face[face];
							*(CoeffN[1]+indexP) = *(CoeffN[1]+indexP) + (*(D[face+6]+indexP))    * diffCoeff_face[face];
							*(CoeffN[2]+indexP) = *(CoeffN[2]+indexP) + (*(D[face+12]+indexP))    * diffCoeff_face[face];
							*(CoeffN[3]+indexP) = *(CoeffN[3]+indexP) + (*(D[face+18]+indexP))    * diffCoeff_face[face];
							*(CoeffN[4]+indexP) = *(CoeffN[4]+indexP) + (*(D[face+24]+indexP))    * diffCoeff_face[face];
							*(CoeffN[5]+indexP) = *(CoeffN[5]+indexP) + (*(D[face+30]+indexP))    * diffCoeff_face[face];
						}
	
						*(CoeffN[6]+indexP) = (*(DN[0]+indexP))    * diffCoeff_face[0] + (*(DN[1]+indexP))    * diffCoeff_face[3];
						*(CoeffN[7]+indexP) = (*(DN[2]+indexP))    * diffCoeff_face[0] + (*(DN[3]+indexP))    * diffCoeff_face[2];
						
						*(CoeffN[8]+indexP) = (*(DN[4]+indexP))    * diffCoeff_face[0] + (*(DN[5]+indexP))    * diffCoeff_face[5];
						*(CoeffN[9]+indexP) = (*(DN[6]+indexP))    * diffCoeff_face[0] + (*(DN[7]+indexP))    * diffCoeff_face[4];
						
						*(CoeffN[10]+indexP) =(*(DN[8]+indexP))    * diffCoeff_face[1] + (*(DN[9]+indexP))    * diffCoeff_face[3];
						*(CoeffN[11]+indexP) =(*(DN[10]+indexP))   * diffCoeff_face[1] +(*(DN[11]+indexP))    * diffCoeff_face[2];
						
						*(CoeffN[12]+indexP) = (*(DN[12]+indexP))    * diffCoeff_face[1] + (*(DN[13]+indexP))    * diffCoeff_face[5];
						*(CoeffN[13]+indexP) = (*(DN[14]+indexP))    * diffCoeff_face[1] + (*(DN[15]+indexP))    * diffCoeff_face[4];
						
						*(CoeffN[14]+indexP) = (*(DN[16]+indexP))    * diffCoeff_face[5] + (*(DN[17]+indexP))    * diffCoeff_face[3];
						*(CoeffN[15]+indexP) = (*(DN[18]+indexP))    * diffCoeff_face[5] + (*(DN[19]+indexP))    * diffCoeff_face[2];
						
						*(CoeffN[16]+indexP) = (*(DN[20]+indexP))    * diffCoeff_face[4] + (*(DN[21]+indexP))    * diffCoeff_face[3];
						*(CoeffN[17]+indexP) = (*(DN[22]+indexP))    * diffCoeff_face[4] + (*(DN[23]+indexP))    * diffCoeff_face[2];
						
						indexP++;
					}
				}
			}
		}
	}
	else if(base->flagMPFV==1)
	{	
		for (iZone=0;iZone<numZones;iZone++)
		{
	
			zone = base->CgnsZone+iZone;	
			if(vartag==4)
			{	
				diffdia = (zone->grid+g)->DiffusionCoeff_P;
				for(i=0;i<18;i++)
				{
					CoeffN[i] = (zone->grid+g)->DiffusionCoeff_nb[i];
				}
				for(i=0;i<6;i++)
				{
					Diag[i] = (zone->grid+g)->Diffusion_Diag[i];
				}
				for(i=0;i<36;i++)
				{
					D[i] = (zone->grid+g)->Diffusion_nb[i];
				}
				for(i=0;i<24;i++)
				{
					DN[i] = (zone->grid+g)->Diffusion_Nonortho[i];
				}
				
			}
			else
			{
				diffdia = (zone->grid+g)->DiffusionCoeffPhi_P;
				
				for(i=0;i<18;i++)
				{
					CoeffN[i] = (zone->grid+g)->DiffusionCoeffPhi_nb[i];
				}
				for(i=0;i<6;i++)
				{
					Diag[i] = (zone->grid+g)->Diffusion_Diag[i];
				}
				for(i=0;i<36;i++)
				{
					D[i] = (zone->grid+g)->Diffusion_nb[i];
				}
				for(i=0;i<24;i++)
				{
					DN[i] = (zone->grid+g)->Diffusion_Nonortho[i];
				}
			}
			
			
			ptrSize = zone->Size;
			beta=zone->ZGeom->Beta;
			volumeDomain = zone->ZGeom->Volume;
			
			NI=*(ptrSize+0);
			NJ=*(ptrSize+1);
			NK=*(ptrSize+2);
			
			K_Increment = (NI-1) *(NJ-1);
			J_Increment = NI-1;
			flagdiffusionCoeff = MakeDiffusionCoeffFlag(zone,vartag);
			//printf("\nVartag=%d and flagdiffcoeff=%d",vartag,flagdiffusionCoeff);
			MakeSlabPtrArray(zone ,slabRho,710);
			Domain=MakeDomainArray(zone,710);
			
			MakeSlabPtrArray(zone,slabVolume,5);/*Called for making slabs of Volumes*/
			flagrho = zone->Flowsol->flagRho;
			
			///This loop is added to avoid unnecessary accessing the volumes and other multiplications in case of a variable with Constant diffusion coefficient.  Still this extension is to be made to variable with variable diffusion coefficient.  This is not been done right now as the extension needs more time and memory.  This needs to be discussed with sir before finalizing this subroutine.
	
			if(flagdiffusionCoeff==0)
			{
				diffusionCoeffdomain = MakeDiffusionCoeffDomain(zone,vartag);
				diff_coeff = *diffusionCoeffdomain;
				
				diffdia_GEOM = (zone->grid+g)->DiffusionCoeff_P;
				
				for(i=0;i<18;i++)
				{
					CoeffN_GEOM[i] = (zone->grid+g)->DiffusionCoeff_nb[i];
				}
				indexP =0;
				for(k=0;k<NK-1;k++)
					for(j=0;j<NJ-1;j++)
						for(i=0;i<NI-1;i++)
				{
					*(diffdia  +indexP) = diff_coeff * *(diffdia_GEOM  +indexP);
					*(CoeffN[0]+indexP) = diff_coeff * *(CoeffN_GEOM[0]+indexP);
					*(CoeffN[1]+indexP) = diff_coeff * *(CoeffN_GEOM[1]+indexP);
					*(CoeffN[2]+indexP) = diff_coeff * *(CoeffN_GEOM[2]+indexP);
					*(CoeffN[3]+indexP) = diff_coeff * *(CoeffN_GEOM[3]+indexP);
					*(CoeffN[4]+indexP) = diff_coeff * *(CoeffN_GEOM[4]+indexP);
					*(CoeffN[5]+indexP) = diff_coeff * *(CoeffN_GEOM[5]+indexP);
					*(CoeffN[6]+indexP) = diff_coeff * *(CoeffN_GEOM[6]+indexP);
					*(CoeffN[7]+indexP) = diff_coeff * *(CoeffN_GEOM[7]+indexP);
					*(CoeffN[8]+indexP) = diff_coeff * *(CoeffN_GEOM[8]+indexP);
					*(CoeffN[9]+indexP) = diff_coeff * *(CoeffN_GEOM[9]+indexP);
					*(CoeffN[10]+indexP) = diff_coeff * *(CoeffN_GEOM[10]+indexP);
					*(CoeffN[11]+indexP) = diff_coeff * *(CoeffN_GEOM[11]+indexP);
					*(CoeffN[12]+indexP) = diff_coeff * *(CoeffN_GEOM[12]+indexP);
					*(CoeffN[13]+indexP) = diff_coeff * *(CoeffN_GEOM[13]+indexP);
					*(CoeffN[14]+indexP) = diff_coeff * *(CoeffN_GEOM[14]+indexP);
					*(CoeffN[15]+indexP) = diff_coeff * *(CoeffN_GEOM[15]+indexP);
					*(CoeffN[16]+indexP) = diff_coeff * *(CoeffN_GEOM[16]+indexP);
					*(CoeffN[17]+indexP) = diff_coeff * *(CoeffN_GEOM[17]+indexP);
					indexP++;
				}
			}
			else
			{
				diffusionCoeffdomain = MakeDiffusionCoeffDomain(zone,vartag);
				MakeDiffusionCoeffSlabPtrArray(zone,slabDiffusionCoeff,vartag);
// 				printf("The [zone--%s]--[diff value --%lf]--[vartag--%d]\n",zone->zonename,*diffusionCoeffdomain,vartag);
// 				printf("The slabDiffusionCoeff at 0th slab in zone %s is %lf\n",zone->zonename,*slabDiffusionCoeff[0]);
// 				printf("The slabDiffusionCoeff at 1st slab in zone %s is %lf\n",zone->zonename,*slabDiffusionCoeff[1]);
// 				printf("The slabDiffusionCoeff at 2nd slab in zone %s is %lf\n",zone->zonename,*slabDiffusionCoeff[2]);
// 				printf("The slabDiffusionCoeff at 3rd slab in zone %s is %lf\n",zone->zonename,*slabDiffusionCoeff[3]);
// 				printf("The slabDiffusionCoeff at 4th slab in zone %s is %lf\n",zone->zonename,*slabDiffusionCoeff[4]);
// 				printf("The slabDiffusionCoeff at 5th slab in zone %s is %lf\n",zone->zonename,*slabDiffusionCoeff[5]);
				if((base->BaseParameters->flag2d == 1)&&(vartag==4))
				{
					indexP=0;
					for(k=0;k<NK-1;k++)
					for(j=0;j<NJ-1;j++)
					for(i=0;i<NI-1;i++)
					{
						if(i!=0)
						{
							Vol[0]=*(volumeDomain+indexP-1);
							diffCoeff_w=*(diffusionCoeffdomain+indexP-1);
							
							
						}
						else
						{
							Vol[0] =*(slabVolume[1]+j);
							diffCoeff_w =*(slabDiffusionCoeff[1]+j);
							
						}
						if(i!=(NI-2))
						{
							Vol[1]=*(volumeDomain+indexP+1);
							diffCoeff_e=*(diffusionCoeffdomain+indexP+1);
							
						}
						else
						{
							Vol[1]= *(slabVolume[0]+j);
							diffCoeff_e=*(slabDiffusionCoeff[0]+j);
							
						}
	
						if(j!=0)
						{
							Vol[2]=*(volumeDomain+indexP-J_Increment);
							diffCoeff_s=*(diffusionCoeffdomain+indexP-J_Increment);
							
						}
						else
						{
							Vol[2]=*(slabVolume[3]+i);
							diffCoeff_s=*(slabDiffusionCoeff[3]+i);
							
						}
	
						if(j!=(NJ-2))
						{
							Vol[3]=*(volumeDomain+indexP+J_Increment);
							diffCoeff_n=*(diffusionCoeffdomain+indexP+J_Increment);
							
						}
						else
						{
							Vol[3]=*(slabVolume[2]+i);
							diffCoeff_n=*(slabDiffusionCoeff[2]+i);
							
						}
						diffCoeff_p=*(diffusionCoeffdomain+indexP);
						vol_p=*(volumeDomain+indexP);
						
						
						
						diffCoeff_face[0] = 1/((diffCoeff_w*vol_p + diffCoeff_p*Vol[0])/(vol_p+Vol[0]));
						diffCoeff_face[1] = 1/((diffCoeff_e*vol_p + diffCoeff_p*Vol[1])/(vol_p+Vol[1]));
						diffCoeff_face[2] = 1/((diffCoeff_s*vol_p + diffCoeff_p*Vol[2])/(vol_p+Vol[2]));
						diffCoeff_face[3] = 1/((diffCoeff_n*vol_p + diffCoeff_p*Vol[3])/(vol_p+Vol[3]));
					
						*(diffdia  +indexP) = 0.0;
						*(CoeffN[0]+indexP) = 0.0;
						*(CoeffN[1]+indexP) = 0.0;
						*(CoeffN[2]+indexP) = 0.0;
						*(CoeffN[3]+indexP) = 0.0;
						
						for (face=0;face<4;face++)
						{
							*(diffdia  +indexP) = *(diffdia  +indexP) + (*(Diag[face]+indexP)) * diffCoeff_face[face];
							*(CoeffN[0]+indexP) = *(CoeffN[0]+indexP) + (*(D[face]+indexP))    * diffCoeff_face[face];
							*(CoeffN[1]+indexP) = *(CoeffN[1]+indexP) + (*(D[face+6]+indexP))    * diffCoeff_face[face];
							*(CoeffN[2]+indexP) = *(CoeffN[2]+indexP) + (*(D[face+12]+indexP))    * diffCoeff_face[face];
							
							*(CoeffN[3]+indexP) = *(CoeffN[3]+indexP) + (*(D[face+18]+indexP))    * diffCoeff_face[face];
						}
						*(CoeffN[6]+indexP) = (*(DN[0]+indexP))    * diffCoeff_face[0] + (*(DN[1]+indexP))    * diffCoeff_face[3];
						*(CoeffN[7]+indexP) = (*(DN[2]+indexP))    * diffCoeff_face[0] + (*(DN[3]+indexP))    * diffCoeff_face[2];
						*(CoeffN[10]+indexP) =(*(DN[8]+indexP))    * diffCoeff_face[1] + (*(DN[9]+indexP))    * diffCoeff_face[3];
						*(CoeffN[11]+indexP) =(*(DN[10]+indexP))   * diffCoeff_face[1] +(*(DN[11]+indexP))    * diffCoeff_face[2];
	
						indexP++;
					}
		
				}
				else if((base->BaseParameters->flag2d == 1)&&(vartag==600||vartag==601||vartag==35))
				{
					indexP=0;
					for(k=0;k<NK-1;k++)
					for(j=0;j<NJ-1;j++)
					for(i=0;i<NI-1;i++)
					{
						//printf("Diff Domain=%e\n",*(diffusionCoeffdomain+indexP));
						if(i!=0)
						{
							Vol[0]=*(volumeDomain+indexP-1);
							diffCoeff_w=*(diffusionCoeffdomain+indexP-1);
							
							
						}
						else
						{
							Vol[0] =*(slabVolume[1]+j);
							diffCoeff_w =*(slabDiffusionCoeff[1]+j);
							
						}
						if(i!=(NI-2))
						{
							Vol[1]=*(volumeDomain+indexP+1);
							diffCoeff_e=*(diffusionCoeffdomain+indexP+1);
							
						}
						else
						{
							Vol[1]= *(slabVolume[0]+j);
							diffCoeff_e=*(slabDiffusionCoeff[0]+j);
							
						}
	
						if(j!=0)
						{
							Vol[2]=*(volumeDomain+indexP-J_Increment);
							diffCoeff_s=*(diffusionCoeffdomain+indexP-J_Increment);
							
						}
						else
						{
							Vol[2]=*(slabVolume[3]+i);
							diffCoeff_s=*(slabDiffusionCoeff[3]+i);
							
						}
	
						if(j!=(NJ-2))
						{
							Vol[3]=*(volumeDomain+indexP+J_Increment);
							diffCoeff_n=*(diffusionCoeffdomain+indexP+J_Increment);
							
						}
						else
						{
							Vol[3]=*(slabVolume[2]+i);
							diffCoeff_n=*(slabDiffusionCoeff[2]+i);
							
						}
						diffCoeff_p=*(diffusionCoeffdomain+indexP);
						vol_p=*(volumeDomain+indexP);
						
						if(vartag==35)
							r_p=1;
						else 
							r_p=1/(*(Domain+indexP));//printf("\nrho p=%lf",r_p);
						
						
						/*printf("d_p=%e\n",diffCoeff_p);
						printf("d_e=%e\n",diffCoeff_e);
						printf("d_w=%e\n",diffCoeff_w);
						printf("d_n=%e\n",diffCoeff_n);
						printf("d_s=%e\n",diffCoeff_s);	*/					
						
						diffCoeff_face[0] = ((diffCoeff_w*vol_p + diffCoeff_p*Vol[0])/(vol_p+Vol[0]));
						diffCoeff_face[1] = ((diffCoeff_e*vol_p + diffCoeff_p*Vol[1])/(vol_p+Vol[1]));
						diffCoeff_face[2] = ((diffCoeff_s*vol_p + diffCoeff_p*Vol[2])/(vol_p+Vol[2]));
						diffCoeff_face[3] = ((diffCoeff_n*vol_p + diffCoeff_p*Vol[3])/(vol_p+Vol[3]));
						//printf("0=%e,1=%e,2=%e,3=%e\n",diffCoeff_face[0],diffCoeff_face[1],diffCoeff_face[2],diffCoeff_face[3]);
						*(diffdia  +indexP) = 0.0;
						*(CoeffN[0]+indexP) = 0.0;
						*(CoeffN[1]+indexP) = 0.0;
						*(CoeffN[2]+indexP) = 0.0;
						*(CoeffN[3]+indexP) = 0.0;
						
						for (face=0;face<4;face++)
						{
							*(diffdia  +indexP) = *(diffdia  +indexP) + (*(Diag[face]+indexP)) * diffCoeff_face[face]*r_p;
							*(CoeffN[0]+indexP) = *(CoeffN[0]+indexP) + (*(D[face]+indexP))    * diffCoeff_face[face]*r_p;
							*(CoeffN[1]+indexP) = *(CoeffN[1]+indexP) + (*(D[face+6]+indexP))    * diffCoeff_face[face]*r_p;
							*(CoeffN[2]+indexP) = *(CoeffN[2]+indexP) + (*(D[face+12]+indexP))    * diffCoeff_face[face]*r_p;
							
							*(CoeffN[3]+indexP) = *(CoeffN[3]+indexP) + (*(D[face+18]+indexP))    * diffCoeff_face[face]*r_p;
						}
						*(CoeffN[6]+indexP) = (*(DN[0]+indexP))    * diffCoeff_face[0] + (*(DN[1]+indexP))    * diffCoeff_face[3]*r_p;
						*(CoeffN[7]+indexP) = (*(DN[2]+indexP))    * diffCoeff_face[0] + (*(DN[3]+indexP))    * diffCoeff_face[2]*r_p;
						*(CoeffN[10]+indexP) =(*(DN[8]+indexP))    * diffCoeff_face[1] + (*(DN[9]+indexP))    * diffCoeff_face[3]*r_p;
						*(CoeffN[11]+indexP) =(*(DN[10]+indexP))   * diffCoeff_face[1] +(*(DN[11]+indexP))    * diffCoeff_face[2]*r_p;
						
						indexP++;
					}
		
				}
				else if(base->BaseParameters->flag2d == 1)
				{	
					
				
					indexP=0;
					for(k=0;k<NK-1;k++)
					for(j=0;j<NJ-1;j++)
					for(i=0;i<NI-1;i++)
					{
						if(i!=0)
						{
							Vol[0]=*(volumeDomain+indexP-1);
							diffCoeff_w=*(diffusionCoeffdomain+indexP-1);
							r[0] = 1/(*(Domain+indexP-1));
						}
						else
						{
							Vol[0] =*(slabVolume[1]+j);
							diffCoeff_w =*(slabDiffusionCoeff[1]+j);
							r[0] =1/(*(slabRho[1]+j));
						}
						if(i!=(NI-2))
						{
							Vol[1]=*(volumeDomain+indexP+1);
							diffCoeff_e=*(diffusionCoeffdomain+indexP+1);
							r[1]=1/(*(Domain+indexP+1));
						}
						else
						{
							Vol[1]= *(slabVolume[0]+j);
							diffCoeff_e=*(slabDiffusionCoeff[0]+j);
							r[1]=1/(*(slabRho[0]+j));
						}

						if(j!=0)
						{
							Vol[2]=*(volumeDomain+indexP-J_Increment);
							diffCoeff_s=*(diffusionCoeffdomain+indexP-J_Increment);
							r[2]=1/(*(Domain+indexP-J_Increment));
						}
						else
						{
							Vol[2]=*(slabVolume[3]+i);
							diffCoeff_s=*(slabDiffusionCoeff[3]+i);
							r[2]=1/(*(slabRho[3]+i));
						}

						if(j!=(NJ-2))
						{
							Vol[3]=*(volumeDomain+indexP+J_Increment);
							diffCoeff_n=*(diffusionCoeffdomain+indexP+J_Increment);
							r[3]=1/(*(Domain+indexP+J_Increment));
							
						}
						else
						{
							Vol[3]=*(slabVolume[2]+i);
							diffCoeff_n=*(slabDiffusionCoeff[2]+i);
							r[3]=1/(*(slabRho[2]+i));
						}
						
						if(!(i==0 && j==0))
						{
							if(i==0)
							{
			
								r[7]=1/(*(slabRho[1]+(NJ-1)*k+(j-1)));
							}
							else
							{
								if(j==0)
								{
				
									r[7]=1/(*(slabRho[3]+(NI-1)*k+(i-1)));
								}
								else
								{
									r[7]=1/(*(Domain+indexP-J_Increment-1));
				
								}
							}
						}
						else
							r[7]=0.0;
						if(!(i==0 && j==(NJ-2)))
						{
							if(i==0)
							{
			
								r[6]=1/(*(slabRho[1]+(NJ-1)*k+(j+1)));
							}
							else
							{
								if(j==(NJ-2))
								{
				
									r[6]=1/(*(slabRho[2]+(NI-1)*k+(i-1)));
								}
								else
								{
									r[6]=1/(*(Domain+indexP+J_Increment-1));
				
								}
							}
						}
						else
							r[6]=0.0;
						if(!(i==0 && k==0))
						{
							if(i==0)
							{
			
								r[9]=1/(*(slabRho[1]+(NJ-1)*(k-1)+j));
							}
							else
							{
								if(k==0)
								{
				
									r[9]=1/(*(slabRho[5]+(NI-1)*j+(i-1)));
								}
								else
								{
									r[9]=1/(*(Domain+indexP-K_Increment-1));
				
								}
							}
		
						}
						else
							r[9]=0.0;
						if(!(i==0 && k==(NK-2)))
						{
							if(i==0)
							{
			
								r[8]=1/(*(slabRho[1]+(NJ-1)*(k+1)+j));
							}
							else
							{
								if(k==(NK-2))
								{
				
									r[8]=1/(*(slabRho[4]+(NI-1)*j+(i-1)));
								}
								else
								{
									r[8]=1/(*(Domain+indexP+K_Increment-1));
				
								}
							}
		
						}
						else
							r[8]=0.0;
						if(!(i==(NI-2) && j==0))
						{
							if(i==(NI-2))
							{
			
								r[11]=1/(*(slabRho[0]+(NJ-1)*k+(j-1)));
							}
							else
							{
								if(j==0)
								{
				
									r[11]=1/(*(slabRho[3]+(NI-1)*k+(i+1)));
								}
								else
								{
									r[11]=1/(*(Domain+indexP-J_Increment+1));
				
								}
							}
						}
						else
							r[11]=0.0;
						if(!(i==(NI-2) && j==(NJ-2)))
						{
							if(i==(NI-2))
							{
			
								r[10]=1/(*(slabRho[0]+(NJ-1)*k+(j+1)));
							}
							else
							{
								if(j==(NJ-2))
								{
				
									r[10]=1/(*(slabRho[2]+(NI-1)*k+(i+1)));
								}
								else
								{
									r[10]=1/(*(Domain+indexP+J_Increment+1));
				
								}
							}
						}
						else
							r[10]=0.0;
						if(!(i==(NI-2) && k==0))
						{
							if(i==(NI-2))
							{
			
								r[13]=1/(*(slabRho[0]+(NJ-1)*(k-1)+j));
							}
							else
							{
								if(k==0)
								{
				
									r[13]=1/(*(slabRho[5]+(NI-1)*j+(i+1)));
								}
								else
								{
									r[13]=1/(*(Domain+indexP-K_Increment+1));
				
								}
							}
						}
						else
							r[13]=0.0;
						if(!(i==(NI-2) && k==(NK-2)))
						{
							if(i==(NI-2))
							{
			
								r[12]= 1/(*(slabRho[0]+(NJ-1)*(k+1)+j));
							}
							else
							{
								if(k==(NK-2))
								{
				
									r[12]=1/(*(slabRho[4]+(NI-1)*j+(i+1)));
								}
								else
								{
									r[12]=1/(*(Domain+indexP+K_Increment+1));
				
								}
							}
						}
						else
							r[12]=0.0;
							
						if(!(k==0 && j==0))
						{
							if(k==0)
							{
			
								r[17]=1/(*(slabRho[5]+(NI-1)*(j-1)+i));
							}
							else
							{
								if(j==0)
								{
				
									r[17]=1/(*(slabRho[3]+(NI-1)*(k-1)+i));
								}
								else
								{
									r[17]=1/(*(Domain+indexP-J_Increment-K_Increment));
				
								}
							}
		
						}
						else
							r[17]=0.0;
						if(!(k==0 && j==(NJ-2)))
						{
							if(k==0 )
							{
			
								r[16]=1/(*(slabRho[5]+(NI-1)*(j+1)+i));
							}
							else
							{
								if(j==(NJ-2))
								{
				
									r[16]=1/(*(slabRho[2]+(NI-1)*(k-1)+i));
								}
								else
								{
									r[16]=1/(*(Domain+indexP+J_Increment-K_Increment));
				
								}
							}
		
						}
						else
							r[16]=0.0;
						if(!(k==(NK-2) && j==0))
						{
							if(k==(NK-2))
							{
			
								r[15]=1/(*(slabRho[4]+(NI-1)*(j-1)+i));
							}
							else
							{
								if(j==0)
								{
				
									r[15]=1/(*(slabRho[3]+(NI-1)*(k+1)+i));
								}
								else
								{
									r[15]=1/(*(Domain+indexP-J_Increment+K_Increment));
				
								}
							}
						}
						else
							r[15]=0.0;
						if(!(k==(NK-2) && j==(NJ-2)))
						{
							if(k==(NK-2))
							{
			
								r[14]=1/(*(slabRho[4]+(NI-1)*(j+1)+i));
							}
							else
							{
								if(j==(NJ-2))
								{
				
									r[14]=1/(*(slabRho[2]+(NI-1)*(k+1)+i));
								}
								else
								{
									r[14]=1/(*(Domain+indexP+J_Increment+K_Increment));
				
								}
							}
						}
						else
							r[14]=0.0;
						diffCoeff_p=*(diffusionCoeffdomain+indexP);
						
						vol_p=*(volumeDomain+indexP);
						
						r_p=1/(*(Domain+indexP));
					
						diffCoeff_face[0] = (diffCoeff_w*vol_p + diffCoeff_p*Vol[0])/(vol_p+Vol[0]);
						diffCoeff_face[1] = (diffCoeff_e*vol_p + diffCoeff_p*Vol[1])/(vol_p+Vol[1]);
						diffCoeff_face[2] = (diffCoeff_s*vol_p + diffCoeff_p*Vol[2])/(vol_p+Vol[2]);
						diffCoeff_face[3] = (diffCoeff_n*vol_p + diffCoeff_p*Vol[3])/(vol_p+Vol[3]);
				
						*(diffdia  +indexP) = 0.0;
						*(CoeffN[0]+indexP) = 0.0;
						*(CoeffN[1]+indexP) = 0.0;
						*(CoeffN[2]+indexP) = 0.0;
						*(CoeffN[3]+indexP) = 0.0;
					
						for (face=0;face<4;face++)
						{
							*(diffdia  +indexP) = *(diffdia  +indexP) + ((*(Diag[face]+indexP)) * diffCoeff_face[face])*r_p;
							*(CoeffN[0]+indexP) = *(CoeffN[0]+indexP) + ((*(D[face]+indexP))    * diffCoeff_face[face])*r[0];
							*(CoeffN[1]+indexP) = *(CoeffN[1]+indexP) + ((*(D[face+6]+indexP))    * diffCoeff_face[face])*r[1];
							*(CoeffN[2]+indexP) = *(CoeffN[2]+indexP) +( (*(D[face+12]+indexP))    * diffCoeff_face[face])*r[2];
						
							*(CoeffN[3]+indexP) = *(CoeffN[3]+indexP) + ((*(D[face+18]+indexP))    * diffCoeff_face[face])*r[3];
						}
						*(CoeffN[6]+indexP) = ((*(DN[0]+indexP))    * diffCoeff_face[0] + (*(DN[1]+indexP))    * diffCoeff_face[3])*r[6];
						*(CoeffN[7]+indexP) = ((*(DN[2]+indexP))    * diffCoeff_face[0] + (*(DN[3]+indexP))    * diffCoeff_face[2])*r[7];
						*(CoeffN[10]+indexP) =((*(DN[8]+indexP))    * diffCoeff_face[1] + (*(DN[9]+indexP))    * diffCoeff_face[3])*r[10];
						*(CoeffN[11]+indexP) =((*(DN[10]+indexP))   * diffCoeff_face[1] +(*(DN[11]+indexP))    * diffCoeff_face[2])*r[11];

						indexP++;
					}
	
					
					
				}
				else if((base->BaseParameters->flag2d == 0)&&(vartag==4))
				{
					indexP=0;
					for(k=0;k<NK-1;k++)
					for(j=0;j<NJ-1;j++)
					for(i=0;i<NI-1;i++)
					{
						if(i!=0)
						{
							Vol[0]=*(volumeDomain+indexP-1);
							diffCoeff_w=*(diffusionCoeffdomain+indexP-1);
						}
						else
						{
							Vol[0] =*(slabVolume[1]+(NJ-1)*k+j);
							diffCoeff_w =*(slabDiffusionCoeff[1]+(NJ-1)*k+j);
						}
						
						if(i!=(NI-2))
						{
							Vol[1]=*(volumeDomain+indexP+1);
							diffCoeff_e=*(diffusionCoeffdomain+indexP+1);
						}
						else
						{
							Vol[1]= *(slabVolume[0]+(NJ-1)*k+j);
							diffCoeff_e=*(slabDiffusionCoeff[0]+(NJ-1)*k+j);
						}
	
						if(j!=0)
						{
							Vol[2]=*(volumeDomain+indexP-J_Increment);
							diffCoeff_s=*(diffusionCoeffdomain+indexP-J_Increment);
						}
						else
						{
							Vol[2]=*(slabVolume[3]+i);
							diffCoeff_s=*(slabDiffusionCoeff[3]+(NI-1)*k+i);
						}
	
						if(j!=(NJ-2))
						{
							Vol[3]=*(volumeDomain+indexP+J_Increment);
							diffCoeff_n=*(diffusionCoeffdomain+indexP+J_Increment);
						}
						else
						{
							Vol[3]=*(slabVolume[2]+i);
							diffCoeff_n=*(slabDiffusionCoeff[2]+(NI-1)*k+i);
						}
						if(k!=0)
						{
							Vol[4]=*(volumeDomain+indexP-K_Increment);
							diffCoeff_b=*(diffusionCoeffdomain+indexP-K_Increment);
						}
						else
						{
							Vol[4]=*(slabVolume[5]+(NI-1)*j+i);
							diffCoeff_b=*(slabDiffusionCoeff[5]+(NI-1)*j+i);
						}
						
						if(k!=(NK-2))
						{
							Vol[5] = *(volumeDomain+indexP+K_Increment);
							diffCoeff_t=*(diffusionCoeffdomain+indexP+K_Increment);
						}
						else
						{
							Vol[5] = *(slabVolume[4]+(NI-1)*j+i);
							diffCoeff_t=*(slabDiffusionCoeff[4]+(NI-1)*j+i);
						}
						
						diffCoeff_p=*(diffusionCoeffdomain+indexP);
						vol_p=*(volumeDomain+indexP);
						
						diffCoeff_face[0] = 1/((diffCoeff_w*vol_p + diffCoeff_p*Vol[0])/(vol_p+Vol[0]));
						diffCoeff_face[1] = 1/((diffCoeff_e*vol_p + diffCoeff_p*Vol[1])/(vol_p+Vol[1]));
						diffCoeff_face[2] = 1/((diffCoeff_s*vol_p + diffCoeff_p*Vol[2])/(vol_p+Vol[2]));
						diffCoeff_face[3] = 1/((diffCoeff_n*vol_p + diffCoeff_p*Vol[3])/(vol_p+Vol[3]));
						diffCoeff_face[4] = 1/((diffCoeff_b*vol_p + diffCoeff_p*Vol[4])/(vol_p+Vol[4]));
						diffCoeff_face[5] = 1/((diffCoeff_t*vol_p + diffCoeff_p*Vol[5])/(vol_p+Vol[5]));
						
						*(diffdia  +indexP) = 0.0;
						*(CoeffN[0]+indexP) = 0.0;
						*(CoeffN[1]+indexP) = 0.0;
						*(CoeffN[2]+indexP) = 0.0;
						*(CoeffN[3]+indexP) = 0.0;
						*(CoeffN[4]+indexP) = 0.0;
						*(CoeffN[5]+indexP) = 0.0;
						
						for (face=0;face<6;face++)
						{
							*(diffdia  +indexP) = *(diffdia  +indexP) + (*(Diag[face]+indexP)) * diffCoeff_face[face];
							*(CoeffN[0]+indexP) = *(CoeffN[0]+indexP) + (*(D[face]+indexP))    * diffCoeff_face[face];
							*(CoeffN[1]+indexP) = *(CoeffN[1]+indexP) + (*(D[face+6]+indexP))    * diffCoeff_face[face];
							*(CoeffN[2]+indexP) = *(CoeffN[2]+indexP) + (*(D[face+12]+indexP))    * diffCoeff_face[face];
							*(CoeffN[3]+indexP) = *(CoeffN[3]+indexP) + (*(D[face+18]+indexP))    * diffCoeff_face[face];
							*(CoeffN[4]+indexP) = *(CoeffN[4]+indexP) + (*(D[face+24]+indexP))    * diffCoeff_face[face];
							*(CoeffN[5]+indexP) = *(CoeffN[5]+indexP) + (*(D[face+30]+indexP))    * diffCoeff_face[face];
						}
	
						*(CoeffN[6]+indexP) = (*(DN[0]+indexP))    * diffCoeff_face[0] + (*(DN[1]+indexP))    * diffCoeff_face[3];
						*(CoeffN[7]+indexP) = (*(DN[2]+indexP))    * diffCoeff_face[0] + (*(DN[3]+indexP))    * diffCoeff_face[2];
						
						*(CoeffN[8]+indexP) = (*(DN[4]+indexP))    * diffCoeff_face[0] + (*(DN[5]+indexP))    * diffCoeff_face[5];
						*(CoeffN[9]+indexP) = (*(DN[6]+indexP))    * diffCoeff_face[0] + (*(DN[7]+indexP))    * diffCoeff_face[4];
						
						*(CoeffN[10]+indexP) =(*(DN[8]+indexP))    * diffCoeff_face[1] + (*(DN[9]+indexP))    * diffCoeff_face[3];
						*(CoeffN[11]+indexP) =(*(DN[10]+indexP))   * diffCoeff_face[1] +(*(DN[11]+indexP))    * diffCoeff_face[2];
						
						*(CoeffN[12]+indexP) = (*(DN[12]+indexP))    * diffCoeff_face[1] + (*(DN[13]+indexP))    * diffCoeff_face[5];
						*(CoeffN[13]+indexP) = (*(DN[14]+indexP))    * diffCoeff_face[1] + (*(DN[15]+indexP))    * diffCoeff_face[4];
						
						*(CoeffN[14]+indexP) = (*(DN[16]+indexP))    * diffCoeff_face[5] + (*(DN[17]+indexP))    * diffCoeff_face[3];
						*(CoeffN[15]+indexP) = (*(DN[18]+indexP))    * diffCoeff_face[5] + (*(DN[19]+indexP))    * diffCoeff_face[2];
						
						*(CoeffN[16]+indexP) = (*(DN[20]+indexP))    * diffCoeff_face[4] + (*(DN[21]+indexP))    * diffCoeff_face[3];
						*(CoeffN[17]+indexP) = (*(DN[22]+indexP))    * diffCoeff_face[4] + (*(DN[23]+indexP))    * diffCoeff_face[2];
						
						indexP++;
					}
				}
				else if((base->BaseParameters->flag2d == 0)&&(vartag==600||vartag==601||vartag==35))
				{
					indexP=0;
					for(k=0;k<NK-1;k++)
					for(j=0;j<NJ-1;j++)
					for(i=0;i<NI-1;i++)
					{
						if(i!=0)
						{
							Vol[0]=*(volumeDomain+indexP-1);
							diffCoeff_w=*(diffusionCoeffdomain+indexP-1);
						}
						else
						{
							Vol[0] =*(slabVolume[1]+(NJ-1)*k+j);
							diffCoeff_w =*(slabDiffusionCoeff[1]+(NJ-1)*k+j);
						}
						
						if(i!=(NI-2))
						{
							Vol[1]=*(volumeDomain+indexP+1);
							diffCoeff_e=*(diffusionCoeffdomain+indexP+1);
						}
						else
						{
							Vol[1]= *(slabVolume[0]+(NJ-1)*k+j);
							diffCoeff_e=*(slabDiffusionCoeff[0]+(NJ-1)*k+j);
						}
	
						if(j!=0)
						{
							Vol[2]=*(volumeDomain+indexP-J_Increment);
							diffCoeff_s=*(diffusionCoeffdomain+indexP-J_Increment);
						}
						else
						{
							Vol[2]=*(slabVolume[3]+i);
							diffCoeff_s=*(slabDiffusionCoeff[3]+(NI-1)*k+i);
						}
	
						if(j!=(NJ-2))
						{
							Vol[3]=*(volumeDomain+indexP+J_Increment);
							diffCoeff_n=*(diffusionCoeffdomain+indexP+J_Increment);
						}
						else
						{
							Vol[3]=*(slabVolume[2]+i);
							diffCoeff_n=*(slabDiffusionCoeff[2]+(NI-1)*k+i);
						}
						if(k!=0)
						{
							Vol[4]=*(volumeDomain+indexP-K_Increment);
							diffCoeff_b=*(diffusionCoeffdomain+indexP-K_Increment);
						}
						else
						{
							Vol[4]=*(slabVolume[5]+(NI-1)*j+i);
							diffCoeff_b=*(slabDiffusionCoeff[5]+(NI-1)*j+i);
						}
						
						if(k!=(NK-2))
						{
							Vol[5] = *(volumeDomain+indexP+K_Increment);
							diffCoeff_t=*(diffusionCoeffdomain+indexP+K_Increment);
						}
						else
						{
							Vol[5] = *(slabVolume[4]+(NI-1)*j+i);
							diffCoeff_t=*(slabDiffusionCoeff[4]+(NI-1)*j+i);
						}
						
						diffCoeff_p=*(diffusionCoeffdomain+indexP);
						vol_p=*(volumeDomain+indexP);
						if(vartag==35)
							r_p=1;
						else
							r_p=1/(*(Domain+indexP));
						diffCoeff_face[0] = 1/((diffCoeff_w*vol_p + diffCoeff_p*Vol[0])/(vol_p+Vol[0]));
						diffCoeff_face[1] = 1/((diffCoeff_e*vol_p + diffCoeff_p*Vol[1])/(vol_p+Vol[1]));
						diffCoeff_face[2] = 1/((diffCoeff_s*vol_p + diffCoeff_p*Vol[2])/(vol_p+Vol[2]));
						diffCoeff_face[3] = 1/((diffCoeff_n*vol_p + diffCoeff_p*Vol[3])/(vol_p+Vol[3]));
						diffCoeff_face[4] = 1/((diffCoeff_b*vol_p + diffCoeff_p*Vol[4])/(vol_p+Vol[4]));
						diffCoeff_face[5] = 1/((diffCoeff_t*vol_p + diffCoeff_p*Vol[5])/(vol_p+Vol[5]));
						
						*(diffdia  +indexP) = 0.0;
						*(CoeffN[0]+indexP) = 0.0;
						*(CoeffN[1]+indexP) = 0.0;
						*(CoeffN[2]+indexP) = 0.0;
						*(CoeffN[3]+indexP) = 0.0;
						*(CoeffN[4]+indexP) = 0.0;
						*(CoeffN[5]+indexP) = 0.0;
						
						for (face=0;face<6;face++)
						{
							*(diffdia  +indexP) = *(diffdia  +indexP) + (*(Diag[face]+indexP)) * diffCoeff_face[face]*r_p;
							*(CoeffN[0]+indexP) = *(CoeffN[0]+indexP) + (*(D[face]+indexP))    * diffCoeff_face[face]*r_p;
							*(CoeffN[1]+indexP) = *(CoeffN[1]+indexP) + (*(D[face+6]+indexP))    * diffCoeff_face[face]*r_p;
							*(CoeffN[2]+indexP) = *(CoeffN[2]+indexP) + (*(D[face+12]+indexP))    * diffCoeff_face[face]*r_p;
							*(CoeffN[3]+indexP) = *(CoeffN[3]+indexP) + (*(D[face+18]+indexP))    * diffCoeff_face[face]*r_p;
							*(CoeffN[4]+indexP) = *(CoeffN[4]+indexP) + (*(D[face+24]+indexP))    * diffCoeff_face[face]*r_p;
							*(CoeffN[5]+indexP) = *(CoeffN[5]+indexP) + (*(D[face+30]+indexP))    * diffCoeff_face[face]*r_p;
						}
	
						*(CoeffN[6]+indexP) = (*(DN[0]+indexP))    * diffCoeff_face[0] + (*(DN[1]+indexP))    * diffCoeff_face[3]*r_p;
						*(CoeffN[7]+indexP) = (*(DN[2]+indexP))    * diffCoeff_face[0] + (*(DN[3]+indexP))    * diffCoeff_face[2]*r_p;
						
						*(CoeffN[8]+indexP) = (*(DN[4]+indexP))    * diffCoeff_face[0] + (*(DN[5]+indexP))    * diffCoeff_face[5]*r_p;
						*(CoeffN[9]+indexP) = (*(DN[6]+indexP))    * diffCoeff_face[0] + (*(DN[7]+indexP))    * diffCoeff_face[4]*r_p;
						
						*(CoeffN[10]+indexP) =(*(DN[8]+indexP))    * diffCoeff_face[1] + (*(DN[9]+indexP))    * diffCoeff_face[3]*r_p;
						*(CoeffN[11]+indexP) =(*(DN[10]+indexP))   * diffCoeff_face[1] +(*(DN[11]+indexP))    * diffCoeff_face[2]*r_p;
						
						*(CoeffN[12]+indexP) = (*(DN[12]+indexP))    * diffCoeff_face[1] + (*(DN[13]+indexP))    * diffCoeff_face[5]*r_p;
						*(CoeffN[13]+indexP) = (*(DN[14]+indexP))    * diffCoeff_face[1] + (*(DN[15]+indexP))    * diffCoeff_face[4]*r_p;
						
						*(CoeffN[14]+indexP) = (*(DN[16]+indexP))    * diffCoeff_face[5] + (*(DN[17]+indexP))    * diffCoeff_face[3]*r_p;
						*(CoeffN[15]+indexP) = (*(DN[18]+indexP))    * diffCoeff_face[5] + (*(DN[19]+indexP))    * diffCoeff_face[2]*r_p;
						
						*(CoeffN[16]+indexP) = (*(DN[20]+indexP))    * diffCoeff_face[4] + (*(DN[21]+indexP))    * diffCoeff_face[3]*r_p;
						*(CoeffN[17]+indexP) = (*(DN[22]+indexP))    * diffCoeff_face[4] + (*(DN[23]+indexP))    * diffCoeff_face[2]*r_p;
						
						indexP++;
					}
				}
				
				else if(base->BaseParameters->flag2d == 0)
				{
					indexP=0;
					for(k=0;k<NK-1;k++)
					for(j=0;j<NJ-1;j++)
					for(i=0;i<NI-1;i++)
					{
						if(i!=0)
						{
							Vol[0]=*(volumeDomain+indexP-1);
							diffCoeff_w=*(diffusionCoeffdomain+indexP-1);
							r[0]=1/(*(Domain+indexP-1));
						}
						else
						{
							Vol[0] =*(slabVolume[1]+(NJ-1)*k+j);
							diffCoeff_w =*(slabDiffusionCoeff[1]+(NJ-1)*k+j);
							r[0] =1/(*(slabRho[1]+(NJ-1)*k+j));
						}
						
						if(i!=(NI-2))
						{
							Vol[1]=*(volumeDomain+indexP+1);
							diffCoeff_e=*(diffusionCoeffdomain+indexP+1);
							r[1]=1/(*(Domain+indexP+1));
						}
						else
						{
							Vol[1]= *(slabVolume[0]+(NJ-1)*k+j);
							diffCoeff_e=*(slabDiffusionCoeff[0]+(NJ-1)*k+j);
							r[1]= 1/(*(slabRho[0]+(NJ-1)*k+j));
						}
	
						if(j!=0)
						{
							Vol[2]=*(volumeDomain+indexP-J_Increment);
							diffCoeff_s=*(diffusionCoeffdomain+indexP-J_Increment);
							r[2]=1/(*(Domain+indexP-J_Increment));
						}
						else
						{
							Vol[2]=*(slabVolume[3]+i);
							diffCoeff_s=*(slabDiffusionCoeff[3]+(NI-1)*k+i);
							r[2]=1/(*(slabRho[3]+i));
						}
	
						if(j!=(NJ-2))
						{
							Vol[3]=*(volumeDomain+indexP+J_Increment);
							diffCoeff_n=*(diffusionCoeffdomain+indexP+J_Increment);
							r[3]=1/(*(Domain+indexP+J_Increment));
						}
						else
						{
							Vol[3]=*(slabVolume[2]+i);
							diffCoeff_n=*(slabDiffusionCoeff[2]+(NI-1)*k+i);
							r[3]=1/(*(slabRho[2]+i));
						}
						if(k!=0)
						{
							Vol[4]=*(volumeDomain+indexP-K_Increment);
							diffCoeff_b=*(diffusionCoeffdomain+indexP-K_Increment);
							r[4]=1/(*(Domain+indexP-K_Increment));
						}
						else
						{
							Vol[4]=*(slabVolume[5]+(NI-1)*j+i);
							diffCoeff_b=*(slabDiffusionCoeff[5]+(NI-1)*j+i);
							r[4]=1/(*(slabRho[5]+(NI-1)*j+i));
						}
						
						if(k!=(NK-2))
						{
							Vol[5] = *(volumeDomain+indexP+K_Increment);
							diffCoeff_t=*(diffusionCoeffdomain+indexP+K_Increment);
							r[5] = 1/(*(Domain+indexP+K_Increment));
						}
						else
						{
							Vol[5] = *(slabVolume[4]+(NI-1)*j+i);
							diffCoeff_t=*(slabDiffusionCoeff[4]+(NI-1)*j+i);
							r[5] = 1/(*(slabRho[4]+(NI-1)*j+i));
						}
						if(!(i==0 && j==0))
						{
							if(i==0)
							{
			
								r[7]=1/(*(slabRho[1]+(NJ-1)*k+(j-1)));
							}
							else
							{
								if(j==0)
								{
				
									r[7]=1/(*(slabRho[3]+(NI-1)*k+(i-1)));
								}
								else
								{
									r[7]=1/(*(Domain+indexP-J_Increment-1));
				
								}
							}
						}
						else
							r[7]=0.0;
						if(!(i==0 && j==(NJ-2)))
						{
							if(i==0)
							{
			
								r[6]=1/(*(slabRho[1]+(NJ-1)*k+(j+1)));
							}
							else
							{
								if(j==(NJ-2))
								{
				
									r[6]=1/(*(slabRho[2]+(NI-1)*k+(i-1)));
								}
								else
								{
									r[6]=1/(*(Domain+indexP+J_Increment-1));
				
								}
							}
						}
						else
							r[6]=0.0;
						if(!(i==0 && k==0))
						{
							if(i==0)
							{
			
								r[9]=1/(*(slabRho[1]+(NJ-1)*(k-1)+j));
							}
							else
							{
								if(k==0)
								{
				
									r[9]=1/(*(slabRho[5]+(NI-1)*j+(i-1)));
								}
								else
								{
									r[9]=1/(*(Domain+indexP-K_Increment-1));
				
								}
							}
		
						}
						else
							r[9]=0.0;
						if(!(i==0 && k==(NK-2)))
						{
							if(i==0)
							{
			
								r[8]=1/(*(slabRho[1]+(NJ-1)*(k+1)+j));
							}
							else
							{
								if(k==(NK-2))
								{
				
									r[8]=1/(*(slabRho[4]+(NI-1)*j+(i-1)));
								}
								else
								{
									r[8]=1/(*(Domain+indexP+K_Increment-1));
				
								}
							}
		
						}
						else
							r[8]=0.0;
						if(!(i==(NI-2) && j==0))
						{
							if(i==(NI-2))
							{
			
								r[11]=1/(*(slabRho[0]+(NJ-1)*k+(j-1)));
							}
							else
							{
								if(j==0)
								{
				
									r[11]=1/(*(slabRho[3]+(NI-1)*k+(i+1)));
								}
								else
								{
									r[11]=1/(*(Domain+indexP-J_Increment+1));
				
								}
							}
						}
						else
							r[11]=0.0;
						if(!(i==(NI-2) && j==(NJ-2)))
						{
							if(i==(NI-2))
							{
			
								r[10]=1/(*(slabRho[0]+(NJ-1)*k+(j+1)));
							}
							else
							{
								if(j==(NJ-2))
								{
				
									r[10]=1/(*(slabRho[2]+(NI-1)*k+(i+1)));
								}
								else
								{
									r[10]=1/(*(Domain+indexP+J_Increment+1));
				
								}
							}
						}
						else
							r[10]=0.0;
						if(!(i==(NI-2) && k==0))
						{
							if(i==(NI-2))
							{
			
								r[13]=1/(*(slabRho[0]+(NJ-1)*(k-1)+j));
							}
							else
							{
								if(k==0)
								{
				
									r[13]=1/(*(slabRho[5]+(NI-1)*j+(i+1)));
								}
								else
								{
									r[13]=1/(*(Domain+indexP-K_Increment+1));
				
								}
							}
						}
						else
							r[13]=0.0;
						if(!(i==(NI-2) && k==(NK-2)))
						{
							if(i==(NI-2))
							{
			
								r[12]= 1/(*(slabRho[0]+(NJ-1)*(k+1)+j));
							}
							else
							{
								if(k==(NK-2))
								{
				
									r[12]=1/(*(slabRho[4]+(NI-1)*j+(i+1)));
								}
								else
								{
									r[12]=1/(*(Domain+indexP+K_Increment+1));
				
								}
							}
						}
						else
							r[12]=0.0;	
						if(!(k==0 && j==0))
						{
							if(k==0)
							{
			
								r[17]=1/(*(slabRho[5]+(NI-1)*(j-1)+i));
							}
							else
							{
								if(j==0)
								{
				
									r[17]=1/(*(slabRho[3]+(NI-1)*(k-1)+i));
								}
								else
								{
									r[17]=1/(*(Domain+indexP-J_Increment-K_Increment));
				
								}
							}
		
						}
						else
							r[17]=0.0;
						if(!(k==0 && j==(NJ-2)))
						{
							if(k==0 )
							{
			
								r[16]=1/(*(slabRho[5]+(NI-1)*(j+1)+i));
							}
							else
							{
								if(j==(NJ-2))
								{
				
									r[16]=1/(*(slabRho[2]+(NI-1)*(k-1)+i));
								}
								else
								{
									r[16]=1/(*(Domain+indexP+J_Increment-K_Increment));
				
								}
							}
		
						}
						else
							r[16]=0.0;
						if(!(k==(NK-2) && j==0))
						{
							if(k==(NK-2))
							{
			
								r[15]=1/(*(slabRho[4]+(NI-1)*(j-1)+i));
							}
							else
							{
								if(j==0)
								{
				
									r[15]=1/(*(slabRho[3]+(NI-1)*(k+1)+i));
								}
								else
								{
									r[15]=1/(*(Domain+indexP-J_Increment+K_Increment));
				
								}
							}
						}
						else
							r[15]=0.0;
						if(!(k==(NK-2) && j==(NJ-2)))
						{
							if(k==(NK-2))
							{
			
								r[14]=1/(*(slabRho[4]+(NI-1)*(j+1)+i));
							}
							else
							{
								if(j==(NJ-2))
								{
				
									r[14]=1/(*(slabRho[2]+(NI-1)*(k+1)+i));
								}
								else
								{
									r[14]=1/(*(Domain+indexP+J_Increment+K_Increment));
				
								}
							}
						}
						else
							r[14]=0.0;
						diffCoeff_p=*(diffusionCoeffdomain+indexP);
						vol_p=*(volumeDomain+indexP);
						r_p	 = 1/(*(Domain+indexP));
						
						diffCoeff_face[0] = (diffCoeff_w*vol_p + diffCoeff_p*Vol[0])/(vol_p+Vol[0]);
						diffCoeff_face[1] = (diffCoeff_e*vol_p + diffCoeff_p*Vol[1])/(vol_p+Vol[1]);
						diffCoeff_face[2] = (diffCoeff_s*vol_p + diffCoeff_p*Vol[2])/(vol_p+Vol[2]);
						diffCoeff_face[3] = (diffCoeff_n*vol_p + diffCoeff_p*Vol[3])/(vol_p+Vol[3]);
						diffCoeff_face[4] = (diffCoeff_b*vol_p + diffCoeff_p*Vol[4])/(vol_p+Vol[4]);
						diffCoeff_face[5] = (diffCoeff_t*vol_p + diffCoeff_p*Vol[5])/(vol_p+Vol[5]);
						
						*(diffdia  +indexP) = 0.0;
						*(CoeffN[0]+indexP) = 0.0;
						*(CoeffN[1]+indexP) = 0.0;
						*(CoeffN[2]+indexP) = 0.0;
						*(CoeffN[3]+indexP) = 0.0;
						*(CoeffN[4]+indexP) = 0.0;
						*(CoeffN[5]+indexP) = 0.0;
						
						for (face=0;face<6;face++)
						{
							*(diffdia  +indexP) = *(diffdia  +indexP) + ((*(Diag[face]+indexP)) * diffCoeff_face[face])*r_p;
							*(CoeffN[0]+indexP) = *(CoeffN[0]+indexP) + ((*(D[face]+indexP))    * diffCoeff_face[face])*r[0];
							*(CoeffN[1]+indexP) = *(CoeffN[1]+indexP) + ((*(D[face+6]+indexP))    * diffCoeff_face[face])*r[1];
							*(CoeffN[2]+indexP) = *(CoeffN[2]+indexP) + ((*(D[face+12]+indexP))    * diffCoeff_face[face])*r[2];
							*(CoeffN[3]+indexP) = *(CoeffN[3]+indexP) + ((*(D[face+18]+indexP))    * diffCoeff_face[face])*r[3];
							*(CoeffN[4]+indexP) = *(CoeffN[4]+indexP) + ((*(D[face+24]+indexP))    * diffCoeff_face[face])*r[4];
							*(CoeffN[5]+indexP) = *(CoeffN[5]+indexP) + ((*(D[face+30]+indexP))    * diffCoeff_face[face])*r[5];
						}
	
						*(CoeffN[6]+indexP) = ((*(DN[0]+indexP))    * diffCoeff_face[0] + (*(DN[1]+indexP))    * diffCoeff_face[3])*r[6];
						*(CoeffN[7]+indexP) = ((*(DN[2]+indexP))    * diffCoeff_face[0] + (*(DN[3]+indexP))    * diffCoeff_face[2])*r[7];
						
						*(CoeffN[8]+indexP) = ((*(DN[4]+indexP))    * diffCoeff_face[0] + (*(DN[5]+indexP))    * diffCoeff_face[5])*r[8];
						*(CoeffN[9]+indexP) = ((*(DN[6]+indexP))    * diffCoeff_face[0] + (*(DN[7]+indexP))    * diffCoeff_face[4])*r[9];
						
						*(CoeffN[10]+indexP) =((*(DN[8]+indexP))    * diffCoeff_face[1] + (*(DN[9]+indexP))    * diffCoeff_face[3])*r[10];
						*(CoeffN[11]+indexP) =((*(DN[10]+indexP))   * diffCoeff_face[1] +(*(DN[11]+indexP))    * diffCoeff_face[2])*r[11];
						
						*(CoeffN[12]+indexP) = ((*(DN[12]+indexP))    * diffCoeff_face[1] + (*(DN[13]+indexP))    * diffCoeff_face[5])*r[12];
						*(CoeffN[13]+indexP) = ((*(DN[14]+indexP))    * diffCoeff_face[1] + (*(DN[15]+indexP))    * diffCoeff_face[4])*r[13];
						
						*(CoeffN[14]+indexP) = ((*(DN[16]+indexP))    * diffCoeff_face[5] + (*(DN[17]+indexP))    * diffCoeff_face[3])*r[14];
						*(CoeffN[15]+indexP) = ((*(DN[18]+indexP))    * diffCoeff_face[5] + (*(DN[19]+indexP))    * diffCoeff_face[2])*r[15];
						
						*(CoeffN[16]+indexP) = ((*(DN[20]+indexP))    * diffCoeff_face[4] + (*(DN[21]+indexP))    * diffCoeff_face[3])*r[16];
						*(CoeffN[17]+indexP) = ((*(DN[22]+indexP))    * diffCoeff_face[4] + (*(DN[23]+indexP))    * diffCoeff_face[2])*r[17];
						
						indexP++;
					}
				}
			}
		}	
	}
	return(0);
}

