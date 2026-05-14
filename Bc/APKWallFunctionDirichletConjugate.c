#include "APKWallFunctionDirichletConjugate.h"
/**
The function WallFunctionDirichletConjugate() returns a nonzero value depending on the variable which is calculated 
according to the wall function treatment.
For details refer to Wall function treatment for the k-epsilon model, Page: 81, Section 6.5, Project Report No. 1, Dec 2005.
 */
double WallFunctionDirichletConjugate(Zone *z,Interface *IfacePtr,int Face,int vartag,int slabcellindex,int domaincellindex,double Bv)
{
	double variable = 0.0;
	//int NI,NJ,NK,i,j,k;
	int flagrho,flagviscoeff,flagSpecificHeat,flagDiffcoeff;
	double uPtx,uPty,uPtz,mod_uPt,deltaPlus,nx,ny,nz,uPx,uPy,uPz,*surfareaIndex,deltaP;
	double alpha,uTau;
	double Pr,Prt,a1,a2,beta,Tplus,y1plus,y2plus,Tw,TP;

	if(!(vartag == 600 || vartag == 601 || vartag == 6))
	{
		surfareaIndex = (z->ZGeom)->SurfaceArea+domaincellindex*6*3+3*Face;
		nx = *(surfareaIndex+0)/sqrt(*(surfareaIndex+0) * *(surfareaIndex+0) + *(surfareaIndex+1) * *(surfareaIndex+1) + *(surfareaIndex+2) * *(surfareaIndex+2));
		ny = *(surfareaIndex+1)/sqrt(*(surfareaIndex+0) * *(surfareaIndex+0) + *(surfareaIndex+1) * *(surfareaIndex+1) + *(surfareaIndex+2) * *(surfareaIndex+2));
		nz = *(surfareaIndex+2)/sqrt(*(surfareaIndex+0) * *(surfareaIndex+0) + *(surfareaIndex+1) * *(surfareaIndex+1) + *(surfareaIndex+2) * *(surfareaIndex+2));
	
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
		else
		{
			uPx = *(z->Flowsol->u+domaincellindex);
			uPy = *(z->Flowsol->v+domaincellindex);
			if(z->Flowsol->w!=NULL)
				uPz = *(z->Flowsol->w+domaincellindex);
			else
				uPz = 0.0;
		}
	
		uPtx = uPx - (uPx*nx + uPy*ny + uPz*nz)*nx;
		uPty = uPy - (uPx*nx + uPy*ny + uPz*nz)*ny;
		if(z->Flowsol->w!=NULL)
			uPtz = uPz - (uPx*nx + uPy*ny + uPz*nz)*nz;
		else
			uPtz = 0.0;
	
		mod_uPt = sqrt(uPtx*uPtx+uPty*uPty+uPtz*uPtz);
		if(mod_uPt == 0)
		{
			mod_uPt = 1.0e-10;
		}
	}

	flagrho = z->Flowsol->flagRho;
	flagviscoeff = z->Flowsol->flagViscoeff;
	flagSpecificHeat=z->Flowsol->flagSpecificHeat;
	flagDiffcoeff=z->Flowsol->flagDiffcoeff;

	uTau = *(IfacePtr->uTau+slabcellindex);
	deltaP = *(IfacePtr->deltaP+slabcellindex);

	deltaPlus = uTau * deltaP * *(z->Flowsol->Rho+flagrho*domaincellindex)/ *(z->Flowsol->Viscoeff+flagviscoeff*domaincellindex);
	if(deltaPlus < 11.067)	//<---- Scallable Wall Function Treatment; Ref: Heat Transfer Predictions using advanced turbulence models, [CFX], Mentor
	{
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
			//Prandtl number should be calculated based on the interface values of viscosity, specific heat and thermal conductivity.
			Pr=*(z->Flowsol->Viscoeff+flagviscoeff*domaincellindex)*(*(z->Flowsol->SpecificHeat+flagSpecificHeat*domaincellindex))/ (*(z->Flowsol->Diffcoeff+flagDiffcoeff*domaincellindex));
			y1plus=10.0/(pow(Pr,(1.0/3.0)));
			a1=Prt/(pow(10.0,3.0));
			a2=15.0*pow(Pr,(2.0/3.0));
			y2plus=sqrt(z->Flowsol->Tcoeff->Chi/a1);
			beta=a2-(Prt/(2.0*z->Flowsol->Tcoeff->Chi))*(1.0+log(z->Flowsol->Tcoeff->Chi/a1));

			//Why to calculate 'deltaPlus' again ... Is it necessary? No! Scallable WFT should be 
			//applied to temperature also as and when it is applied to momentum equations, for consistency between flow solution and temperature solution <---Narsimha.
			//deltaPlus = uTau * *(IfacePtr->deltaP+slabcellindex) * *(z->Flowsol->Rho+flagrho*domaincellindex)/(*(z->Flowsol->Viscoeff+flagviscoeff*domaincellindex));
			if(deltaPlus < y1plus)
			{
				Tplus = Pr * deltaPlus;
			}
			else 
			{
				if(deltaPlus < y2plus)
					Tplus = a2 - Prt/(2.0*a1*deltaPlus*deltaPlus);
				else
					Tplus = (Prt/(z->Flowsol->Tcoeff->Chi)) * log(deltaPlus) + beta;
			}
			//Why to use "PrevFlowsol->T" instead of "Flowsol->T" ?
			TP = *(z->PrevFlowsol->T+domaincellindex);
			variable= TP+(Prt/((z->Flowsol->Tcoeff->Chi)*Tplus))*(Tw-TP);
		}
		break;

		case 600:	//Does this case need to be written? No! check <----- Narsimha.
		{
			variable = *(z->Flowsol->Scal[0]+domaincellindex);
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
	
	return(variable);
}
