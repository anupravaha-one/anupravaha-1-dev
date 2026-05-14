#ifndef __APKMain_h
#define __APKMain_h

#include<stdio.h>
#include<stdlib.h>
#include<malloc.h>
#include<string.h>

#include "Common/APKData_Structures.h"
#include "Common/APKprintfluxes.h"
#include "Common/APKprintslabs.h"
#include "Common/APKprintsolver.h"
#include "Common/APKprintu4bc.h"
#include "Common/APKfile_refresh.h"
#include "Common/APKinputFunc.h"
#include "Common/APKprintmassflux.h"
#include "Common/APKCallMakeSlabPtrArray.h"

#include "Common/APKfind_zone_face_info.h"
#include "Common/APKCalcIFACE_data.h"

#include "Geometry/APKCalculateGeometry.h"
#include "Geometry/APKUpdateInterfaceStruct.h"

#include "Bc/APKBc_input.h"
#include "Bc/APKTimeDepend.h"


#include "Io/cgnslib.h"
#include "Io/APKReadCgnsFile.h"
#include "Io/APKReadCgnsFileRestart.h"
#include "Io/APKReadInitialData.h"
#include "Io/APKWriteCgnsFile.h"
#include "Io/APKWriteCgnsFile_vertex.h"
#include "Io/APKRead_Solution.h"

#include "Initialize/APKinitializeSolidLiquid.h"
#include "Initialize/Combustion/APKInitializeNonSolveCombVariables.h"

#include "Solver/APKEvaluateVarCoefficients.h"
#include "Solver/APKRefreshMassVelocities.h"
#include "Solver/APKsolverCoupled.h"
#include "Solver/APKsolverSegregated.h"
#include "Solver/APKsolverSemiCoupled.h"
#include "Solver/APKsolverCoupled2d.h"
#include "Solver/APKsolverSemiCoupled2d.h"
#include "Solver/APKRefreshSurfaceVolumes.h"
#include "Solver/APKRefreshVariableCoefficients.h"
#include "Solver/APKRefreshSurfaceCells.h"
#include "Solver/APKmassflux.h"/*To be removed later, not in use*/
#include "Solver/APKmassfluxstar.h"/*To be removed later, not in use*/
#include "Solver/APKvolumeflux.h"
#include "Solver/APKvolumefluxstar.h"
#include "Solver/APKPressureForce.h"

#include "Solver/APKCalcNewVolumeFlux.h"
#include "Solver/APKCalcNewVolumefluxRestart.h"
#include "Solver/ESP/APKSolverElectric.h"
#include "Solver/APKfillprevdensity.h"
#include "Solver/APKGridParameters1.h"
#include "Solver/APKGridParameters2.h"
#include "Solver/Radiation/APKRadiationSolverFiniteVolume3D.h"
#include "Solver/Radiation/APKRadiationSolverFiniteVolume2D.h"
#include "Solver/Turbulence/APKTURBCalcMuEff.h"
#include "Solver/Turbulence/APKTURBCalcTurbconductivity.h"
#include "Solver/Turbulence/APKTURBCalcTurbViscosity.h"
 
#include "Solver/ParticleTracking/APKFindMinMaxCoordinatesOfZones.h"
#include "Solver/ParticleTracking/APKInputParticle.h"
#include "Solver/ParticleTracking/APKParticleTracking.h"
#include "Solver/LMMHD/APKMagneticflux.h"

int main(int argc,char *argv[]);

#endif

