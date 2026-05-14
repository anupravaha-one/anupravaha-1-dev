#include "APKCalcShearStress.h"

/**
This function Calculates Wall noraml shear stress which has three components along x,y and z directions.
*/

int CalcShearStress(Base *base)
{
	FILE *fp;
	int i,j,k,dir,nZone,index[3];
	int NI,NJ,NK,numZbc,zBc,indexP;
	Zone *zone;
	double *surfarea,surfAreaComp1,surfAreaComp2,surfAreaComp3,Area,*vols,*beta;
	double unitVectorComp1,unitVectorComp2,unitVectorComp3;
	double *slabVelocityU[6],*slabVelocityV[6],*slabVelocityW[6],*slabVolume[6];
	double FdiffU[6],FdiffV[6],FdiffW[6];
	double Chi,viscosity,shearStress1,shearStress2,shearStress3;
	int flagVisc;
	
	fp = fopen("SOLUTIONS/shearStress.dat","w");
	for(nZone=0;nZone<base->Number_of_zones;nZone++)
	{
		zone=(base->CgnsZone+nZone);
		numZbc=zone->Number_of_BC_Patches;
		
		NI = *(zone->Size+0);
		NJ = *(zone->Size+1);
		NK = *(zone->Size+2);
		
		MakeSlabPtrArray(zone,slabVelocityU,1);//Called for making slabs of u
		MakeSlabPtrArray(zone,slabVelocityV,2);//Called for making slabs of v 
		MakeSlabPtrArray(zone,slabVelocityW,3);//Called for making slabs of w 
		MakeSlabPtrArray(zone,slabVolume,5);//Called for making slabs of Volumes
		
		flagVisc=zone->Flowsol->flagViscoeff;
		
		for(dir=1;dir<=6;dir++)
		{
			for(zBc=0;zBc<numZbc;zBc++)
			{
				if((zone->ZoneBc+zBc)->Direction==dir)
				{
					fprintf(fp,"ZoneBc dir->%d\n",dir);
					fprintf(fp,"~~~~~~~~~Boundary Name: %s ~~~~~~~~~~\n",((zone->ZoneBc+zBc)->BCName));
					
					switch(dir)
					{
						case 1://East
							i=NI-2;
							vols= zone->ZGeom->Volume;
							beta= zone->ZGeom->Beta;
							surfarea= zone->ZGeom->SurfaceArea;
							
							for (k=0;k<NK-1;k++) 
							for (j=0;j<NJ-1;j++)
							{
								indexP = (NJ-1)*(NI-1)*k+(NI-1)*j+i;
								surfAreaComp1 = *(surfarea+indexP*6*3+3*1+0); //Note SurfArea is outward normal to the surface.
								surfAreaComp2 = *(surfarea+indexP*6*3+3*1+1);
								surfAreaComp3 = *(surfarea+indexP*6*3+3*1+2);
								Area = sqrt(surfAreaComp1*surfAreaComp1 + surfAreaComp2*surfAreaComp2 + surfAreaComp3*surfAreaComp3);
								
								unitVectorComp1 = -1.0*surfAreaComp1/Area; // Negative sign since we want inward normal
								unitVectorComp2 = -1.0*surfAreaComp2/Area;
								unitVectorComp3 = -1.0*surfAreaComp3/Area;
								
								index[0]=k,index[1]=j,index[2]=i;
								Diffusion(index,zone->Flowsol->u,vols,zone->Size,slabVelocityU,slabVolume,beta,FdiffU);
								Diffusion(index,zone->Flowsol->v,vols,zone->Size,slabVelocityV,slabVolume,beta,FdiffV);
								Diffusion(index,zone->Flowsol->w,vols,zone->Size,slabVelocityW,slabVolume,beta,FdiffW);
								// Note FdiffU, FdiffV, FdiffW are outward fluxes, but we need inward fluxes.
								// IMPORTANT---> For all the negative signs see the reason written above
								
								Chi = ((-1.0*FdiffU[1]/Area)*unitVectorComp1) + ((-1.0*FdiffV[1]/Area)*unitVectorComp2) + ((-1.0*FdiffW[1]/Area)*unitVectorComp3);
								viscosity = *(zone->Flowsol->Viscoeff+indexP*flagVisc);
								
								shearStress1 = viscosity*((-1.0*FdiffU[1]/Area) - Chi*unitVectorComp1);
								shearStress2 = viscosity*((-1.0*FdiffV[1]/Area) - Chi*unitVectorComp2);
								shearStress3 = viscosity*((-1.0*FdiffW[1]/Area) - Chi*unitVectorComp3);
								fprintf(fp,"%d %d %d %12.8lf %12.8lf %12.8lf %12.8lf\n",k,j,i,shearStress1,shearStress2,shearStress3,Chi);
							}
						break;
						
						case 2://West
							i=0;
							vols= zone->ZGeom->Volume;
							beta= zone->ZGeom->Beta;
							surfarea= zone->ZGeom->SurfaceArea;
							
							for (k=0;k<NK-1;k++) 
							for (j=0;j<NJ-1;j++)
							{
								indexP = (NJ-1)*(NI-1)*k+(NI-1)*j+i;
								surfAreaComp1 = *(surfarea+indexP*6*3+3*0+0); //Note SurfArea is outward normal to the surface.
								surfAreaComp2 = *(surfarea+indexP*6*3+3*0+1);
								surfAreaComp3 = *(surfarea+indexP*6*3+3*0+2);
								Area = sqrt(surfAreaComp1*surfAreaComp1 + surfAreaComp2*surfAreaComp2 + surfAreaComp3*surfAreaComp3);
								
								unitVectorComp1 = -1.0*surfAreaComp1/Area; // Negative sign sice we want inward normal
								unitVectorComp2 = -1.0*surfAreaComp2/Area;
								unitVectorComp3 = -1.0*surfAreaComp3/Area;
								
								index[0]=k,index[1]=j,index[2]=i;
								Diffusion(index,zone->Flowsol->u,vols,zone->Size,slabVelocityU,slabVolume,beta,FdiffU);
								Diffusion(index,zone->Flowsol->v,vols,zone->Size,slabVelocityV,slabVolume,beta,FdiffV);
								Diffusion(index,zone->Flowsol->w,vols,zone->Size,slabVelocityW,slabVolume,beta,FdiffW);
								// Note FdiffU, FdiffV, FdiffW are outward fluxes, but we need inward fluxes.
								// IMPORTANT---> For all the negative signs see the reason written above
								
								Chi = ((-1.0*FdiffU[0]/Area)*unitVectorComp1) + ((-1.0*FdiffV[0]/Area)*unitVectorComp2) + ((-1.0*FdiffW[0]/Area)*unitVectorComp3);
								viscosity = *(zone->Flowsol->Viscoeff+indexP*flagVisc);
								
								shearStress1 = viscosity*((-1.0*FdiffU[0]/Area) - Chi*unitVectorComp1);
								shearStress2 = viscosity*((-1.0*FdiffV[0]/Area) - Chi*unitVectorComp2);
								shearStress3 = viscosity*((-1.0*FdiffW[0]/Area) - Chi*unitVectorComp3);
								fprintf(fp,"%d %d %d %12.8lf %12.8lf %12.8lf %12.8lf\n",k,j,i,shearStress1,shearStress2,shearStress3,Chi);
							}
						break;
						
						case 3://North
							j=NJ-2;
							vols= zone->ZGeom->Volume;
							beta= zone->ZGeom->Beta;
							surfarea= zone->ZGeom->SurfaceArea;
							
							for(k=0;k<NK-1;k++) 
							for(i=0;i<NI-1;i++)
							{
								indexP = (NJ-1)*(NI-1)*k+(NI-1)*j+i;
								surfAreaComp1 = *(surfarea+indexP*6*3+3*3+0); //Note SurfArea is outward normal to the surface.
								surfAreaComp2 = *(surfarea+indexP*6*3+3*3+1);
								surfAreaComp3 = *(surfarea+indexP*6*3+3*3+2);
								Area = sqrt(surfAreaComp1*surfAreaComp1 + surfAreaComp2*surfAreaComp2 + surfAreaComp3*surfAreaComp3);
								
								unitVectorComp1 = -1.0*surfAreaComp1/Area; // Negative sign sice we want inward normal
								unitVectorComp2 = -1.0*surfAreaComp2/Area;
								unitVectorComp3 = -1.0*surfAreaComp3/Area;
								
								index[0]=k,index[1]=j,index[2]=i;
								Diffusion(index,zone->Flowsol->u,vols,zone->Size,slabVelocityU,slabVolume,beta,FdiffU);
								Diffusion(index,zone->Flowsol->v,vols,zone->Size,slabVelocityV,slabVolume,beta,FdiffV);
								Diffusion(index,zone->Flowsol->w,vols,zone->Size,slabVelocityW,slabVolume,beta,FdiffW);
								// Note FdiffU, FdiffV, FdiffW are outward fluxes, but we need inward fluxes.
								// IMPORTANT---> For all the negative signs see the reason written above
								
								Chi = ((-1.0*FdiffU[3]/Area)*unitVectorComp1) + ((-1.0*FdiffV[3]/Area)*unitVectorComp2) + ((-1.0*FdiffW[3]/Area)*unitVectorComp3);
								viscosity = *(zone->Flowsol->Viscoeff+indexP*flagVisc);
								
								shearStress1 = viscosity*((-1.0*FdiffU[3]/Area) - Chi*unitVectorComp1);
								shearStress2 = viscosity*((-1.0*FdiffV[3]/Area) - Chi*unitVectorComp2);
								shearStress3 = viscosity*((-1.0*FdiffW[3]/Area) - Chi*unitVectorComp3);
								fprintf(fp,"%d %d %d %12.8lf %12.8lf %12.8lf %12.8lf\n",k,j,i,shearStress1,shearStress2,shearStress3,Chi);
							}
						break;
						
						case 4://South
							j=0;
							vols= zone->ZGeom->Volume;
							beta= zone->ZGeom->Beta;
							surfarea= zone->ZGeom->SurfaceArea;
							
							for(k=0;k<NK-1;k++) 
							for(i=0;i<NI-1;i++)
							{
								indexP = (NJ-1)*(NI-1)*k+(NI-1)*j+i;
								surfAreaComp1 = *(surfarea+indexP*6*3+3*2+0); //Note SurfArea is outward normal to the surface.
								surfAreaComp2 = *(surfarea+indexP*6*3+3*2+1);
								surfAreaComp3 = *(surfarea+indexP*6*3+3*2+2);
								Area = sqrt(surfAreaComp1*surfAreaComp1 + surfAreaComp2*surfAreaComp2 + surfAreaComp3*surfAreaComp3);
								
								unitVectorComp1 = -1.0*surfAreaComp1/Area; // Negative sign sice we want inward normal
								unitVectorComp2 = -1.0*surfAreaComp2/Area;
								unitVectorComp3 = -1.0*surfAreaComp3/Area;
								
								index[0]=k,index[1]=j,index[2]=i;
								Diffusion(index,zone->Flowsol->u,vols,zone->Size,slabVelocityU,slabVolume,beta,FdiffU);
								Diffusion(index,zone->Flowsol->v,vols,zone->Size,slabVelocityV,slabVolume,beta,FdiffV);
								Diffusion(index,zone->Flowsol->w,vols,zone->Size,slabVelocityW,slabVolume,beta,FdiffW);
								// Note FdiffU, FdiffV, FdiffW are outward fluxes, but we need inward fluxes.
								// IMPORTANT---> For all the negative signs see the reason written above
								
								Chi = ((-1.0*FdiffU[2]/Area)*unitVectorComp1) + ((-1.0*FdiffV[2]/Area)*unitVectorComp2) + ((-1.0*FdiffW[2]/Area)*unitVectorComp3);
								viscosity = *(zone->Flowsol->Viscoeff+indexP*flagVisc);
								
								shearStress1 = viscosity*((-1.0*FdiffU[2]/Area) - Chi*unitVectorComp1);
								shearStress2 = viscosity*((-1.0*FdiffV[2]/Area) - Chi*unitVectorComp2);
								shearStress3 = viscosity*((-1.0*FdiffW[2]/Area) - Chi*unitVectorComp3);
								fprintf(fp,"%d %d %d %12.8lf %12.8lf %12.8lf %12.8lf\n",k,j,i,shearStress1,shearStress2,shearStress3,Chi);
							}
						break;
						
						case 5://Top
							k=NK-2;
							vols= zone->ZGeom->Volume;
							beta= zone->ZGeom->Beta;
							surfarea= zone->ZGeom->SurfaceArea;
							
							for(j=0;j<NJ-1;j++) 
							for(i=0;i<NI-1;i++)
							{
								indexP = (NJ-1)*(NI-1)*k+(NI-1)*j+i;
								surfAreaComp1 = *(surfarea+indexP*6*3+3*5+0); //Note SurfArea is outward normal to the surface.
								surfAreaComp2 = *(surfarea+indexP*6*3+3*5+1);
								surfAreaComp3 = *(surfarea+indexP*6*3+3*5+2);
								Area = sqrt(surfAreaComp1*surfAreaComp1 + surfAreaComp2*surfAreaComp2 + surfAreaComp3*surfAreaComp3);
								
								unitVectorComp1 = -1.0*surfAreaComp1/Area; // Negative sign sice we want inward normal
								unitVectorComp2 = -1.0*surfAreaComp2/Area;
								unitVectorComp3 = -1.0*surfAreaComp3/Area;
								
								index[0]=k,index[1]=j,index[2]=i;
								Diffusion(index,zone->Flowsol->u,vols,zone->Size,slabVelocityU,slabVolume,beta,FdiffU);
								Diffusion(index,zone->Flowsol->v,vols,zone->Size,slabVelocityV,slabVolume,beta,FdiffV);
								Diffusion(index,zone->Flowsol->w,vols,zone->Size,slabVelocityW,slabVolume,beta,FdiffW);
								// Note FdiffU, FdiffV, FdiffW are outward fluxes, but we need inward fluxes.
								// IMPORTANT---> For all the negative signs see the reason written above
								
								Chi = ((-1.0*FdiffU[5]/Area)*unitVectorComp1) + ((-1.0*FdiffV[5]/Area)*unitVectorComp2) + ((-1.0*FdiffW[5]/Area)*unitVectorComp3);
								viscosity = *(zone->Flowsol->Viscoeff+indexP*flagVisc);
								
								shearStress1 = viscosity*((-1.0*FdiffU[5]/Area) - Chi*unitVectorComp1);
								shearStress2 = viscosity*((-1.0*FdiffV[5]/Area) - Chi*unitVectorComp2);
								shearStress3 = viscosity*((-1.0*FdiffW[5]/Area) - Chi*unitVectorComp3);
								fprintf(fp,"%d %d %d %12.8lf %12.8lf %12.8lf %12.8lf\n",k,j,i,shearStress1,shearStress2,shearStress3,Chi);
							}
						break;
						
						case 6://Bottom
							k=0;
							vols= zone->ZGeom->Volume;
							beta= zone->ZGeom->Beta;
							surfarea= zone->ZGeom->SurfaceArea;
							
							for(j=0;j<NJ-1;j++) 
							for(i=0;i<NI-1;i++)
							{
								indexP = (NJ-1)*(NI-1)*k+(NI-1)*j+i;
								surfAreaComp1 = *(surfarea+indexP*6*3+3*4+0); //Note SurfArea is outward normal to the surface.
								surfAreaComp2 = *(surfarea+indexP*6*3+3*4+1);
								surfAreaComp3 = *(surfarea+indexP*6*3+3*4+2);
								Area = sqrt(surfAreaComp1*surfAreaComp1 + surfAreaComp2*surfAreaComp2 + surfAreaComp3*surfAreaComp3);
								
								unitVectorComp1 = -1.0*surfAreaComp1/Area; // Negative sign sice we want inward normal
								unitVectorComp2 = -1.0*surfAreaComp2/Area;
								unitVectorComp3 = -1.0*surfAreaComp3/Area;
								
								index[0]=k,index[1]=j,index[2]=i;
								Diffusion(index,zone->Flowsol->u,vols,zone->Size,slabVelocityU,slabVolume,beta,FdiffU);
								Diffusion(index,zone->Flowsol->v,vols,zone->Size,slabVelocityV,slabVolume,beta,FdiffV);
								Diffusion(index,zone->Flowsol->w,vols,zone->Size,slabVelocityW,slabVolume,beta,FdiffW);
								// Note FdiffU, FdiffV, FdiffW are outward fluxes, but we need inward fluxes.
								// IMPORTANT---> For all the negative signs see the reason written above
								
								Chi = ((-1.0*FdiffU[4]/Area)*unitVectorComp1) + ((-1.0*FdiffV[4]/Area)*unitVectorComp2) + ((-1.0*FdiffW[4]/Area)*unitVectorComp3);
								viscosity = *(zone->Flowsol->Viscoeff+indexP*flagVisc);
								
								shearStress1 = viscosity*((-1.0*FdiffU[4]/Area) - Chi*unitVectorComp1);
								shearStress2 = viscosity*((-1.0*FdiffV[4]/Area) - Chi*unitVectorComp2);
								shearStress3 = viscosity*((-1.0*FdiffW[4]/Area) - Chi*unitVectorComp3);
								fprintf(fp,"%d %d %d %12.8lf %12.8lf %12.8lf %12.8lf\n",k,j,i,shearStress1,shearStress2,shearStress3,Chi);
							}
						break;
						
						default:
							printf("Direction case not found and exiting from CalcShearStress.c\n");
							exit(1);
						break;
					}
				}//End of if zoneBC->dir statement
			}// End of for loop over zoneBC
		}// End of for loop over direction.
	}// End of for loop over zones
	
	return(0);
}



