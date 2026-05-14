#include "APKTURBCalcmoduPt.h"
/**
This function  is called in CalcUtau when k-Omega model is used.  This calculates |uPt| (the tangential
velocity) at centre of the first cell adjacent to a boundary.
 */
double CalcmoduPt(Zone *z,int domainCellIndex,int Face)
{
	double uPtx,uPty,uPtz=0.0,nx,ny,nz,uPx,uPy,uPz=0.0,*surfareaIndex,moduPt;

	surfareaIndex = (z->ZGeom)->SurfaceArea+domainCellIndex*6*3+3*Face;
	uPx = *(z->Flowsol->u+domainCellIndex);
	uPy = *(z->Flowsol->v+domainCellIndex);
	
	if(z->Flowsol->w!=NULL)
	uPz = *(z->Flowsol->w+domainCellIndex);
	
	nx = *(surfareaIndex+0)/sqrt(*(surfareaIndex+0) * *(surfareaIndex+0) + *(surfareaIndex+1) * *(surfareaIndex+1) + *(surfareaIndex+2) * *(surfareaIndex+2));
	ny = *(surfareaIndex+1)/sqrt(*(surfareaIndex+0) * *(surfareaIndex+0) + *(surfareaIndex+1) * *(surfareaIndex+1) + *(surfareaIndex+2) * *(surfareaIndex+2));
	nz = *(surfareaIndex+2)/sqrt(*(surfareaIndex+0) * *(surfareaIndex+0) + *(surfareaIndex+1) * *(surfareaIndex+1) + *(surfareaIndex+2) * *(surfareaIndex+2));
	uPtx = uPx - (uPx*nx + uPy*ny + uPz*nz)*nx;
	uPty = uPy - (uPx*nx + uPy*ny + uPz*nz)*ny;
	
	if(z->Flowsol->w!=NULL)
	uPtz = uPz - (uPx*nx + uPy*ny + uPz*nz)*nz;
	

	moduPt = sqrt(uPtx*uPtx+uPty*uPty+uPtz*uPtz);
	
	return(moduPt);

}


