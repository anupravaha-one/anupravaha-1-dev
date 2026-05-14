#ifndef __APKsolverSemiCoupled_h
#define __APKsolverSemiCoupled_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<malloc.h>

#include"APKTimeDepend.h"

#include"APKData_Structures.h"
#include"APKprintsolver.h"
#include"APKprintslabs.h"

#include"APKWriteCgnsFile.h"

#include"APKCalcVorticity.h"
#include"APKLimitVariable.h"
#include"APKSolveTurbVariables.h"
#include"APKTURBfreeGamma.h"
#include"APKTURBCalcGamma.h"
#include"APKTURBCalcTurbViscosity.h"
#include"APKTURBCalcUtau.h"
#include"APKTURBCalcC1epsi.h"
#include"APKTURBCalcEpsilonWallTerm.h"
#include"APKTURBCalcKWallTerm.h"
#include"APKTURBCalcMuEff.h"
#include"APKTURBCalcTurbulentProduction.h"
#include"APKTURBCalcVariableOneIterK.h"
#include"APKTURBCalcVariableOneIterEps.h"
#include"APKTURBCalcVariableOneIterOmega.h"
#include"APKTURBCalcVariableOneIterOmegaSST.h"
#include"APKTURBCrossDiffusion.h"
#include"APKTURBDistance2Wall.h"
#include"APKTURBfreeCrossDiff.h"
#include"APKTURBCalcTurbconductivity.h"

#include"APKCDcoeffMemAlloc.h"
#include"APKDiffusionCoeffCalculation.h"
#include"APKDiffusionCoeffCalculationForVariables.h"
#include"APKCDcoeffCalculation.h"

#include"APKCalcNusseltNo.h"
#include"APKCalcVariableOneIterConjugate.h"
#include"APKRefreshSurfaceCells.h"
#include"APKShearStress.h"
#include"APKsolverCoupled.h"
#include "APKsolverSegregated.h"
#include"APKPressureIter.h"
#include"APKCalcVariableOneIter.h"
#include"APKFillNewsolToPrevsol.h"
#include"APKmassVelocityOneIter.h"
#include"APKMassFluxConv.h"
#include"APKVarMassFluxConv.h"
#include"APKCalcNewVolumeFlux.h"
#include"APKvolumefluxstar.h"
#include"APKSaveTimeStepData.h"
#include"APKsavevolumeflux.h"
#include"APKFindVolumeFluxDifference.h"
#include"APKEvaluateVarCoefficients.h"
#include"APKRefreshVariableCoefficients.h"
#include"APKEvaluateEnergysource.h"
#include"APKCalcVariable.h"
#include"APKComputeCombustionVariables.h"
#include"APKSolveCombustionVariables.h"

#include"APKPotentialMHDpoissonRHS.h"
#include"APKPotentialMHDIter.h"
#include"APKCalcCurrentDensity.h"
#include"APKCalcUcrossB.h"

#include"APKPressureHeadIter.h"
#include"APKUpdateWaterContent.h"
#include"APKCalDispersionCoeff.h"	
#include"APKCalVelocity.h"
#include"APKCalGradScalerXY.h"
#include"APKvolumeflux.h"

#include"APKCalcConjugate.h"
#include"APKCalcSteady_vel.h"
#include"APKRadiationSolverWithFlow3D.h"
#include "APKprint_monitors.h"
#include "APKpmonitors.h"
#include"APKPressure_Grad_Calculation.h"
#include"APKCalcSource_opti.h"

#include"APKCalcVariable.h"
#include"APKmassVelocity.h"
#include"APKCallCalcCoefficientsForVariables.h"
#include"APKConv_Coeff_Calculation.h"

#include "APKRefreshMassVelocities.h"
#include"APK_File_Open.h"
#include"APKWrite_Solution.h"

#include"APKConj_PressureIter.h"
#include"APKConj_CalcConjugate.h"
#include"APKModiPressureIter.h"

#include"APKDrag.h"
#include"APKvolumefluxDisperse.h"
#include"APKFindMassFluxDisperseDifference.h"
#include"APKCalcVarDropDia.h"
#include"APKEvaporate.h"
#include"APKCalcEvapoConst.h"
#include"APKCalcFaceValue.h"
#include"APKCalcNewVolumeFluxDisperse.h"
#include"APKDisperseVolumeFluxStar.h"
#include"APKFindDisperseVolumeFluxDifference.h"


#include"APKUpdateVarMPFVDiffusion.h"
#include"APKDragMPFV.h"
#include"APKMPFV.h"
#include"APKvolumefluxMPFV.h"

#include"APKinitialize_variable.h"

#include"APKsolverLMMHD.h"
#include"APKLMMHDmass_Source_Calculation.h"
#include"APKLMMHD_Source_Calculation.h"
#include "APKInducedLMMHD_Source_Calculation.h"

#include "APKprintu4bc.h"
#include "APKLiftMPFV.h"
#include "APKUpdateVarDisperseDiffusion.h"


int solverSemiCoupled(Root CgnsRoot);

#endif

