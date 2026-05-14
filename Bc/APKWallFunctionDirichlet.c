#include "APKWallFunctionDirichlet.h"
/**
The function WallFunctionDirichlet() returns a nonzero value depending on the variable which is calculated according to
the wall function treatment. 
For details refer to Wall function treatment for the k-epsilon model, Page: 81, Section 6.5, Project Report No. 1, Dec 2005.
 */
double WallFunctionDirichlet(Zone *z,ZoneBoundary *zBc,int vartag,int slabcellindex, double Bv,int flagMPFV)
{
	double variable = 0.0;
	int domaincellindex,NI,NJ,NK,i,j,k,Face,flagrho,flagviscoeff,flagSpecificHeat,flagDiffcoeff;
	double uPtx,uPty,uPtz,mod_uPt,deltaPlus,nx,ny,nz,uPx,uPy,uPz,*surfareaIndex,deltaP;
	double alpha,uTau;
	double Pr,Prt,a1,a2,beta,Tplus,y1plus,y2plus,Tw,TP;
	double rho,visc;//added for testing purpose-->chandrasekhar
	variable = Bv;
	if(z->Flowsol->flagTurb == 1 && z->flagSolidFluid == 0)
	{
		if(flagMPFV==0)	
		{
			if(zBc->Boundary_Tag == 22)
			{
				NI = *(z->Size+0);
				NJ = *(z->Size+1);
				NK = *(z->Size+2);
	
				switch(zBc->Direction)
				{
					case 1: //IMax face,East face
					{
						k = (slabcellindex/(NJ-1));
						j = (slabcellindex % (NJ-1));
						i = (NI-2);
						Face = 1;
					}
					break;
					case 2: //Imin face, West face
					{
						k = (slabcellindex/(NJ-1));
						j = (slabcellindex % (NJ-1));
						i = 0;
						Face = 0;
					}
					break;
					case 3: //JMax face, North face
					{
						k = (slabcellindex/(NI-1));
						i = (slabcellindex % (NI-1));
						j = (NJ-2);
						Face = 3;
					}
					break;
					case 4: //Jmin face,South face
					{
						k = (slabcellindex/(NI-1));
						i = (slabcellindex % (NI-1));
						j = 0;
						Face = 2;
					}
					break;
					case 5: //KMax face, Top face
					{
						j = (slabcellindex/(NI-1));
						i = (slabcellindex % (NI-1));
						k = (NK-2);
						Face = 5;
					}
					break;
					case 6: //Kmin face, Bottom face
					{
						j = (slabcellindex/(NI-1));
						i = (slabcellindex % (NI-1));
						k = 0;
						Face = 4;
					}
					break;
				}
	
				domaincellindex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
	
				if(!(vartag == 601 || vartag == 6))
				{
					surfareaIndex = (z->ZGeom)->SurfaceArea+domaincellindex*6*3+3*Face;
					nx = *(surfareaIndex+0)/sqrt(*(surfareaIndex+0) * *(surfareaIndex+0) + *(surfareaIndex+1) * *(surfareaIndex+1) + *(surfareaIndex+2) * *(surfareaIndex+2));
					ny = *(surfareaIndex+1)/sqrt(*(surfareaIndex+0) * *(surfareaIndex+0) + *(surfareaIndex+1) * *(surfareaIndex+1) + *(surfareaIndex+2) * *(surfareaIndex+2));
					nz = *(surfareaIndex+2)/sqrt(*(surfareaIndex+0) * *(surfareaIndex+0) + *(surfareaIndex+1) * *(surfareaIndex+1) + *(surfareaIndex+2) * *(surfareaIndex+2));
	
					uPx = *(z->Flowsol->u+domaincellindex);//printf("The value of the u velocity at the boundary cell is %lf\n",uPx);
					uPy = *(z->Flowsol->v+domaincellindex);//printf("The value of the v velocity at the boundary cell is %lf\n",uPy);
					if(z->Flowsol->w!=NULL)
						uPz = *(z->Flowsol->w+domaincellindex);
					else
						uPz=0.0;	
					//printf("The value of the u velocity at the boundary cell is %lf\n",uPz);
					if(vartag%100 == 0)
					{
						//printf("***** Implementing WallFunctionTreatment to MassVelocities vartag=%d  Bv = %lf*****\n", vartag,Bv);
						uPx = *(z->Flowsol->umass+domaincellindex);
						uPy = *(z->Flowsol->vmass+domaincellindex);
						if(z->Flowsol->w!=NULL)
							uPz = *(z->Flowsol->wmass+domaincellindex);
						else
							uPz=0.0;
					}
	
					uPtx = uPx - (uPx*nx + uPy*ny + uPz*nz)*nx;
					uPty = uPy - (uPx*nx + uPy*ny + uPz*nz)*ny;
					if(z->Flowsol->w!=NULL)
						uPtz = uPz - (uPx*nx + uPy*ny + uPz*nz)*nz;
					else
						uPtz=0.0;
	
					mod_uPt = sqrt(uPtx*uPtx+uPty*uPty+uPtz*uPtz);
					
					if(mod_uPt == 0)
					{
						mod_uPt = 1.0e-10;
					}
				}
	
				flagrho = z->Flowsol->flagRho;//printf("The flag for Rho is %d\n",flagrho);
				flagviscoeff = z->Flowsol->flagViscoeff;//printf("The flag for viscosity is %d",flagviscoeff);
				flagSpecificHeat=z->Flowsol->flagSpecificHeat;
				flagDiffcoeff=z->Flowsol->flagDiffcoeff; 
	
				uTau = *(zBc->uTau+slabcellindex);//printf("The value of friction velocity is %lf\n",uTau);
				deltaP = *(zBc->deltaP+slabcellindex);//printf("The value of delta p is %lf\n",deltaP);
	
				deltaPlus = uTau * deltaP * *(z->Flowsol->Rho+flagrho*domaincellindex)/ *(z->Flowsol->Viscoeff+flagviscoeff*domaincellindex);
				//printf("The value of deltaPlus is %lf\n",deltaPlus);
				if(deltaPlus < 11.067)	//<---- Scallable Wall Function Treatment; Ref: Heat Transfer Predictions using advanced turbulence models, [CFX], Mentor
				{
					//printf("I am entering the deltaPlus correction");
					deltaPlus =  11.067;
					deltaP = deltaPlus * *(z->Flowsol->Viscoeff+flagviscoeff*domaincellindex)/(*(z->Flowsol->Rho+flagrho*domaincellindex) * uTau);
				}
				
				switch(vartag)
				{
					case 100:
					case 1:
					{
						alpha = (uTau/(mod_uPt * z->Flowsol->Tcoeff->Chi)) * (log(9.0*deltaPlus)-1.0);
						variable = Bv + fabs(alpha)*uPtx;
					}
					break;
	
					case 200:
					case 2:
					{
						alpha = (uTau/(mod_uPt * z->Flowsol->Tcoeff->Chi)) * (log(9.0*deltaPlus)-1.0);
						variable = Bv + fabs(alpha)*uPty;
					}
					break;
	
					case 300:
					case 3:
					{
						alpha = (uTau/(mod_uPt * z->Flowsol->Tcoeff->Chi)) * (log(9.0*deltaPlus)-1.0);
						variable = Bv + fabs(alpha)*uPtz;
					}
					break;
					case 6://wall function treatment for temperature
					{
						Tw=Bv;
						Prt=1.0;
						Pr=*(z->Flowsol->Viscoeff+flagviscoeff*domaincellindex)*(*(z->Flowsol->SpecificHeat+flagSpecificHeat*domaincellindex))/ (*(z->Flowsol->Diffcoeff+flagDiffcoeff*domaincellindex));
						y1plus=10.0/(pow(Pr,(1.0/3.0)));
						a1=Prt/(pow(10.0,3.0));
						a2=15.0*pow(Pr,(2.0/3.0));
						y2plus=sqrt(z->Flowsol->Tcoeff->Chi/a1);
						beta=a2-(Prt/(2.0*z->Flowsol->Tcoeff->Chi))*(1.0+log(z->Flowsol->Tcoeff->Chi/a1));
	
						//Why to calculate 'deltaPlus' again ... Is it necessary? No! Scallable WFT should be 
						//applied to temperature also as and when it is applied to momentum equations, for consistency between flow solution and temperature solution <---Narsimha.
						//deltaPlus = uTau * *(zBc->deltaP+slabcellindex) * *(z->Flowsol->Rho+flagrho*domaincellindex)/(*(z->Flowsol->Viscoeff+flagviscoeff*domaincellindex));
						if(deltaPlus < y1plus)
						{
							Tplus=Pr*deltaPlus;
						}
						else if(deltaPlus < y2plus && deltaPlus >= y1plus)
						{
							Tplus=a2-Prt/(2.0*a1*pow(deltaPlus,2.0));
						}
						else
						{
							Tplus=(Prt/(z->Flowsol->Tcoeff->Chi))*log(deltaPlus)+beta;
						}
						
						TP = *(z->PrevFlowsol->T+domaincellindex);
						variable = TP + (Prt/((z->Flowsol->Tcoeff->Chi)*Tplus))*(Tw-TP);
						
	// 					printf("TP=%e Tw=%e \n", TP,Tw);
					}
					break;
	
					case 601:
					{
						variable = 2.0 * pow(uTau,3)/(deltaP * z->Flowsol->Tcoeff->Chi);
					}
					break;
	
					default:
						variable = Bv;
						break;
				}
			}
			else
				variable = Bv;
		}	
		else if(flagMPFV==1)	
		{	
			if(zBc->Boundary_Tag == 22)
			{
				NI = *(z->Size+0);
				NJ = *(z->Size+1);
				NK = *(z->Size+2);
	
				switch(zBc->Direction)
				{
					case 1: //IMax face,East face
					{
						k = (slabcellindex/(NJ-1));
						j = (slabcellindex % (NJ-1));
						i = (NI-2);
						Face = 1;
					}
					break;
					case 2: //Imin face, West face
					{
						k = (slabcellindex/(NJ-1));
						j = (slabcellindex % (NJ-1));
						i = 0;
						Face = 0;
					}
					break;
					case 3: //JMax face, North face
					{
						k = (slabcellindex/(NI-1));
						i = (slabcellindex % (NI-1));
						j = (NJ-2);
						Face = 3;
					}
					break;
					case 4: //Jmin face,South face
					{
						k = (slabcellindex/(NI-1));
						i = (slabcellindex % (NI-1));
						j = 0;
						Face = 2;
					}
					break;
					case 5: //KMax face, Top face
					{
						j = (slabcellindex/(NI-1));
						i = (slabcellindex % (NI-1));
						k = (NK-2);
						Face = 5;
					}
					break;
					case 6: //Kmin face, Bottom face
					{
						j = (slabcellindex/(NI-1));
						i = (slabcellindex % (NI-1));
						k = 0;
						Face = 4;
					}
					break;
				}
	
				domaincellindex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
				flagrho = z->Flowsol->flagRho;
				flagviscoeff = z->Flowsol->flagViscoeff;
				rho=(*(z->Flowsol->Rho+flagrho*domaincellindex));
				visc=(*(z->Flowsol->Viscoeff+flagviscoeff*domaincellindex));
				if(!(vartag == 601))
				{
					surfareaIndex = (z->ZGeom)->SurfaceArea+domaincellindex*6*3+3*Face;
					nx = *(surfareaIndex+0)/sqrt(*(surfareaIndex+0) * *(surfareaIndex+0) + *(surfareaIndex+1) * *(surfareaIndex+1) + *(surfareaIndex+2) * *(surfareaIndex+2));
					ny = *(surfareaIndex+1)/sqrt(*(surfareaIndex+0) * *(surfareaIndex+0) + *(surfareaIndex+1) * *(surfareaIndex+1) + *(surfareaIndex+2) * *(surfareaIndex+2));
					nz = *(surfareaIndex+2)/sqrt(*(surfareaIndex+0) * *(surfareaIndex+0) + *(surfareaIndex+1) * *(surfareaIndex+1) + *(surfareaIndex+2) * *(surfareaIndex+2));
	
					uPx = *(z->Flowsol->Momentumx+domaincellindex);
					uPy = *(z->Flowsol->Momentumy+domaincellindex);
					if(z->Flowsol->Momentumz!=NULL)
						uPz = *(z->Flowsol->Momentumz+domaincellindex);
					else
						uPz=0.0;
// 					if(slabcellindex==0)
// 					{
// 						printf("\n\n\nThe vartag of the variable is %d\n",vartag);
// 						printf("The value of the u momentum at the boundary cell is %lf\n",uPx);
// 						printf("The value of the v momentum at the boundary cell is %lf\n",uPy);
// 						printf("The value of the w momentum at the boundary cell is %lf\n",uPz);
// 					}
					
					if(vartag%100 == 0)
					{
						//printf("***** Implementing WallFunctionTreatment to MassVelocities vartag=%d  Bv = %lf*****\n", vartag,Bv);
						uPx = *(z->Flowsol->Momentumxstar+domaincellindex);
						uPy = *(z->Flowsol->Momentumystar+domaincellindex);
						if(z->Flowsol->Momentumzstar!=NULL)
							uPz = *(z->Flowsol->Momentumzstar+domaincellindex);
						else
							uPz=0.0;
					}
					
					uPx=uPx*(1.0/rho);
					uPy=uPy*(1.0/rho);
					uPz=uPz*(1.0/rho);
					
					uPtx = uPx - (uPx*nx + uPy*ny + uPz*nz)*nx;
					uPty = uPy - (uPx*nx + uPy*ny + uPz*nz)*ny;
					if(z->Flowsol->Momentumz!=NULL)
						uPtz = uPz - (uPx*nx + uPy*ny + uPz*nz)*nz;
					else
						uPtz=0.0;
	
					mod_uPt = sqrt(uPtx*uPtx+uPty*uPty+uPtz*uPtz);
					
					if(mod_uPt == 0)
					{
						mod_uPt = 1.0e-10;
					}
				}
// 				if(slabcellindex==0)
// 				{
// 					printf("The var tag is %d\n",vartag);
// 					printf("The boundary name is %s\n",zBc->BCName);
// 					printf("The value of uPtx is %lf\n",uPtx);
// 					printf("The value of uPty is %lf\n",uPty);
// 					printf("The value of uPtz is %lf\n",uPtz);
// 				}
				
				uTau = *(zBc->uTau+slabcellindex);
				deltaP = *(zBc->deltaP+slabcellindex);
	
				deltaPlus = uTau * deltaP * *(z->Flowsol->Rho+flagrho*domaincellindex)/ *(z->Flowsol->Viscoeff+flagviscoeff*domaincellindex);
				
				if(deltaPlus < 11.067)	//<---- Scallable Wall Function Treatment; Ref: Heat Transfer Predictions using advanced turbulence models, [CFX], Mentor
				{
					
					deltaPlus =  11.067;
					deltaP = deltaPlus * *(z->Flowsol->Viscoeff+flagviscoeff*domaincellindex)/(*(z->Flowsol->Rho+flagrho*domaincellindex) * uTau);
				}
// 				if(slabcellindex==0)
// 				{
// 					printf("The flag for Rho is %d\n",flagrho);
// 					printf("The value of Rho is %lf\n",rho);
// 					printf("The flag for viscosity is %d\n",flagviscoeff);
// 					printf("The value of viscosity is %lf\n",visc);
// 					printf("The value of friction velocity is %lf\n",uTau);
// 					printf("The value of delta p is %lf\n",deltaP);
// 					printf("The value of deltaPlus is %lf\n",deltaPlus);
// 					printf("The value of flow velocity at the boundary cell is %lf\n",uPx);
// 					printf("The value of flow velocity at the boundary cell is %lf\n",uPy);
// 					printf("The value of flow velocity at the boundary cell is %lf\n",uPz);
// 				}
				switch(vartag)
				{
					case 3100:
					case 31:
					case 41:
					//case 44:
					{
						alpha = (uTau/(mod_uPt * z->Flowsol->Tcoeff->Chi)) * (log(9.0*deltaPlus)-1.0);
						variable = Bv + fabs(alpha)*uPtx;
					}
					break;
					case 44:
					{
						alpha=(uTau/(mod_uPt * z->Flowsol->Tcoeff->Chi)) * (log(9.0*deltaPlus)-1.0);//0.0;
						variable = Bv + fabs(alpha)*uPtx;
					}
					break;
	
					case 3200:
					case 32:
					case 42:
					//case 45:
					{
						alpha = (uTau/(mod_uPt * z->Flowsol->Tcoeff->Chi)) * (log(9.0*deltaPlus)-1.0);
						variable = Bv + fabs(alpha)*uPty;
					}
					break;
					case 45:
					{
						alpha=(uTau/(mod_uPt * z->Flowsol->Tcoeff->Chi)) * (log(9.0*deltaPlus)-1.0);//0.0;
						variable = Bv + fabs(alpha)*uPty;
					}
					break;
	
					case 3300:
					case 33:
					case 43:
					//case 46:
					{
						alpha = (uTau/(mod_uPt * z->Flowsol->Tcoeff->Chi)) * (log(9.0*deltaPlus)-1.0);
						variable = Bv + fabs(alpha)*uPtz;
					}
					break;
					case 46:
					{
						alpha=(uTau/(mod_uPt * z->Flowsol->Tcoeff->Chi)) * (log(9.0*deltaPlus)-1.0);//0.0;
						variable = Bv + fabs(alpha)*uPtz;
					}
					break;
					
	
					case 601:
					{
						variable = 2.0 * pow(uTau,3)/(deltaP * z->Flowsol->Tcoeff->Chi);
					}
					break;
	
					default:
						variable = Bv;
						break;
				}
				
			}
			else
				variable = Bv;
			
// 			if(slabcellindex==0)
// 			{
// 				printf("The value of %d is %lf\n",vartag,variable);
			// 				
			// 				
// 			}
		}
	}
	
	//printf("%s Dir=%d vartag=%d Boundary value from WFD=%e\n", z->zonename,zBc->Direction,vartag, variable);
	return(variable);
}
