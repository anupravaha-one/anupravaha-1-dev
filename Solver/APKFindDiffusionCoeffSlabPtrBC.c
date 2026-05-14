#include "APKFindDiffusionCoeffSlabPtrBC.h"
/**
The function FindDiffusionCoeffSlabPtrBC() takes the variable tag and pointer corresponding to a zone boundary
as input and returns the corresponding Diffusion Coefficient pointer in that boundary slab.
*/
double *FindDiffusionCoeffSlabPtrBC(Zone *zone,ZoneBoundary *facePtr,int vartag)
{
	double *slab,dummy=1.0;
	switch(vartag)
	{
		case 1:
		case 2:
		case 3:
		case 100:
		case 200:
		case 300:
			if(zone->Flowsol->flagTurb==0)
				slab=facePtr->ZSlab->Viscoeff;
			else
				slab=facePtr->ZSlab->EffViscoeff;
		break;
		
		case 4:	
			if(zone->Flowsol->Momentumx!=NULL)
				slab=facePtr->ZSlab->Rho;
			else
			slab=&dummy;
		break;
		/*case 5:
			slab=facePtr->ZSlab->Volume;
		break;*/
		case 6:
			if(zone->Flowsol->flagTurb==0)
				slab=facePtr->ZSlab->Diffcoeff;
			else
				slab=facePtr->ZSlab->EffTurbConductivity;
		break;
		case 7:
                	if(zone->Flowsol->flagTurbComb==1)
			slab=facePtr->ZSlab->EffScalcoeff;
                        else
			slab=facePtr->ZSlab->Diffcoeff;
		break;
		
		case 8:
			slab=(double *)malloc(sizeof(double));
			slab=&dummy;
		break;
		
		case 9:
			slab=facePtr->ZSlab->ChargeDiffcoeff;
		break;
		
                case 14:

			slab=facePtr->ZSlab->EffScalcoeff;
		break;

                case 15:

			slab=facePtr->ZSlab->EffScalcoeff;
		break;

		case 600:
			slab=facePtr->ZSlab->Scalcoeff[0];
		break;
		
		case 601:
			slab=facePtr->ZSlab->Scalcoeff[1];
		break;
		
		case 602:
			slab=facePtr->ZSlab->Scalcoeff[2];
		break;
			
		case 603:
			slab=facePtr->ZSlab->Scalcoeff[3];
		break;
		
		case 604:
			slab=facePtr->ZSlab->Scalcoeff[4];
		break;
		
		case 605:
			slab=facePtr->ZSlab->Scalcoeff[5];
		break;
		
		case 606:
			slab=facePtr->ZSlab->Scalcoeff[6];
		break;
		
		case 607:
			slab=facePtr->ZSlab->Scalcoeff[7];
		break;
		
		case 608:
			slab=facePtr->ZSlab->Scalcoeff[8];
		break;
		
		case 609:
			slab=facePtr->ZSlab->Scalcoeff[9];
		break;
		
		//Added by sushil for multiphase flow
		case 31:
		case 32:
		case 33:
		case 3100:
		case 3200:
		case 3300:
			if(zone->Flowsol->flagTurb==0)
				slab=facePtr->ZSlab->Viscoeff;
			else
				slab=facePtr->ZSlab->EffViscoeff;
		break;

		case 41:
		case 42:
		case 43:
			if(zone->Flowsol->flagTurb==0)
				slab=facePtr->ZSlab->ViscoeffDummy[0];
			else
				slab=facePtr->ZSlab->DummyeffViscoeffofphase[0];
		break;
		case 44:
		case 45:
		case 46:
			if(zone->Flowsol->flagTurb==0)
				slab=facePtr->ZSlab->ViscoeffDummy[1];
			else
				slab=facePtr->ZSlab->DummyeffViscoeffofphase[1];
		break;
		case 47:
		case 48:
		case 49:
			if(zone->Flowsol->flagTurb==0)
				slab=facePtr->ZSlab->ViscoeffDummy[2];
			else
				slab=facePtr->ZSlab->DummyeffViscoeffofphase[2];
		break;
		case 50:
		case 51:
		case 52:
			if(zone->Flowsol->flagTurb==0)
				slab=facePtr->ZSlab->ViscoeffDummy[3];
			else
				slab=facePtr->ZSlab->DummyeffViscoeffofphase[3];
		break;
		
		case 53:
		case 54:
		case 55:
			if(zone->Flowsol->flagTurb==0)
				slab=facePtr->ZSlab->ViscoeffDummy[4];
			else
				slab=facePtr->ZSlab->DummyeffViscoeffofphase[4];
		break;
		case 34:
			slab=facePtr->ZSlab->DiffVolfrac[0];
		break;	
		case 35:
			slab=facePtr->ZSlab->DiffVolfrac[1];
		break;	
		case 36:
			slab=facePtr->ZSlab->DiffVolfrac[2];
		break;	
		case 37:
			slab=facePtr->ZSlab->DiffVolfrac[3];
		break;	
		case 38:
			slab=facePtr->ZSlab->DiffVolfrac[4];
		break;	
			
		break;
		//end of multiphase flow
		case 68:
			slab=facePtr->ZSlab->DropDiffcoeff;
		break;
		case 62:
		case 63:
		case 64:
		case 6200:
		case 6300:
		case 6400:
			slab=facePtr->ZSlab->DummyViscoeff_disperse;
		break;
		case 65:
			slab=facePtr->ZSlab->DiffVolfrac_disperse;
		break;
		default:
			printf("Proper variable tag is not supplied to FindDiffusionCoeffSlabPtrBC\n");
				exit(1);
		break;
		
	}
	return(slab);
}
