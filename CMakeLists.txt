CMAKE_MINIMUM_REQUIRED(VERSION 3.10)
PROJECT(APK)

# The sources are pre-C23 C (e.g. old-style "int f();" declarations), which GCC 15+
# rejects under its default gnu23 standard. Pin the dialect explicitly.
SET(CMAKE_C_STANDARD 11)
SET(CMAKE_C_EXTENSIONS ON)
INCLUDE_REGULAR_EXPRESSION("^(lex|APK|png|j).*$")

# Don't build anything unless the version of CMake is high enough.
# The matching ELSE/ENDIF should be the last lines in the file.
IF("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" GREATER 1.5)

#-----------------------------------------------------------------------------
# APK version number.  An even minor number corresponds to releases.
SET(APK_MAJOR_VERSION 1)
SET(APK_MINOR_VERSION 0)
SET(APK_BUILD_VERSION 0)
SET(APK_VERSION
    "${APK_MAJOR_VERSION}.${APK_MINOR_VERSION}.${APK_BUILD_VERSION}")

#-----------------------------------------------------------------------------
# Output directories.
SET(LIBRARY_OUTPUT_PATH ${APK_BINARY_DIR}/bin CACHE INTERNAL "Single output directory for building all libraries.")
SET(EXECUTABLE_OUTPUT_PATH ${APK_BINARY_DIR}/bin CACHE INTERNAL "Single output directory for building all executables.")
SET(APK_LIBRARY_DIR ${LIBRARY_OUTPUT_PATH}/${CMAKE_CFG_INTDIR})
SET(APK_EXECUTABLE_DIR ${EXECUTABLE_OUTPUT_PATH}/${CMAKE_CFG_INTDIR})
SET(CXX_TEST_PATH ${EXECUTABLE_OUTPUT_PATH})

if(COMMAND cmake_policy)
      cmake_policy(SET CMP0015 NEW)
endif(COMMAND cmake_policy)

# The Cygwin maintainers have requested that we not define WIN32. Thus developers
# must use the X11 packages to build cygwin Anupravaha.
IF(CYGWIN)
  unset(WIN32)
ENDIF()

# Put "." in the rpath of all APK shared libraries.
LINK_DIRECTORIES(.)


# APK requires special compiler flags on some platforms.
IF(CMAKE_COMPILER_IS_GNUCXX)
  # A GCC compiler.  Quiet warning about strstream deprecation.
  SET(APK_REQUIRED_CXX_FLAGS "${APK_REQUIRED_CXX_FLAGS} -Wno-deprecated")
  IF(WIN32)
    # The platform is gcc on cygwin.
    SET(APK_REQUIRED_CXX_FLAGS "${APK_REQUIRED_CXX_FLAGS} -mwin32")
    SET(APK_REQUIRED_C_FLAGS "${APK_REQUIRED_C_FLAGS} -mwin32")
  ENDIF(WIN32)
  IF(MINGW)
    LINK_LIBRARIES(-lgdi32)
  ENDIF(MINGW)
  IF(CMAKE_SYSTEM MATCHES "SunOS.*")
    # Disable warnings that occur in X11 headers.
    IF(DART_ROOT AND BUILD_TESTING)
      SET(APK_REQUIRED_CXX_FLAGS "${APK_REQUIRED_CXX_FLAGS} -Wno-unknown-pragmas")
      SET(APK_REQUIRED_C_FLAGS "${APK_REQUIRED_C_FLAGS} -Wno-unknown-pragmas")
    ENDIF(DART_ROOT AND BUILD_TESTING)
  ENDIF(CMAKE_SYSTEM MATCHES "SunOS.*")
ELSE(CMAKE_COMPILER_IS_GNUCXX)
  IF(CMAKE_ANSI_CFLAGS)
    SET(APK_REQUIRED_C_FLAGS "${APK_REQUIRED_C_FLAGS} ${CMAKE_ANSI_CFLAGS}")
  ENDIF(CMAKE_ANSI_CFLAGS)
  IF(CMAKE_SYSTEM MATCHES "OSF1-V.*")
     SET(APK_REQUIRED_CXX_FLAGS
         "${APK_REQUIRED_CXX_FLAGS} -timplicit_local -no_implicit_include")
  ENDIF(CMAKE_SYSTEM MATCHES "OSF1-V.*")
  #silence duplicate symbol warnings on AIX
  IF(CMAKE_SYSTEM MATCHES "AIX.*")
    SET(APK_REQUIRED_EXE_LINKER_FLAGS "${APK_REQUIRED_EXE_LINKER_FLAGS} -bhalt:5")
  ENDIF(CMAKE_SYSTEM MATCHES "AIX.*")
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

IF(APPLE)
  SET(APK_REQUIRED_C_FLAGS "${APK_REQUIRED_C_FLAGS} -no-cpp-precomp")
  SET(APK_REQUIRED_CXX_FLAGS "${APK_REQUIRED_CXX_FLAGS} -no-cpp-precomp")
  IF(APK_USE_CARBON)
    SET(APK_REQUIRED_CXX_FLAGS "${APK_REQUIRED_CXX_FLAGS} -fpascal-strings")
  ENDIF(APK_USE_CARBON)
  IF(BUILD_SHARED_LIBS)
    SET(APK_APPLE_SHARED_FLAGS_NEEDED 1)
  ENDIF(BUILD_SHARED_LIBS)
ENDIF(APPLE)

IF(APK_USE_X)
  IF(CMAKE_X_CFLAGS)
    SET(APK_REQUIRED_C_FLAGS "${APK_REQUIRED_C_FLAGS} ${CMAKE_X_CFLAGS}")
    SET(APK_REQUIRED_CXX_FLAGS "${APK_REQUIRED_CXX_FLAGS} ${CMAKE_X_CFLAGS}")
  ENDIF(CMAKE_X_CFLAGS)
ENDIF(APK_USE_X)

#-----------------------------------------------------------------------------
# Configure the default APK_DATA_ROOT for the location of APKData.
FIND_PATH(APK_DATA_ROOT APKData.readme ${APK_SOURCE_DIR}/../APKData $ENV{APK_DATA_ROOT})

#-----------------------------------------------------------------------------
# Configure files with settings for use by the build.
CONFIGURE_FILE(
${APK_SOURCE_DIR}/APKConfigure.h.in
${APK_BINARY_DIR}/APKConfigure.h)

#-----------------------------------------------------------------------------
# The entire APK tree should use the same include path.

# This should be the only INCLUDE_DIRECTORIES command in the entire
# tree, except for the CMake, Utilities, and Examples directories.  We
# need to do this in one place to make sure the order is correct.
INCLUDE_DIRECTORIES(
  ${APK_INCLUDE_DIRS_BUILD_TREE}
  ${APK_INCLUDE_DIRS_SOURCE_TREE}
  ${APK_INCLUDE_DIRS_BUILD_TREE_CXX}
  ${APK_INCLUDE_DIRS_SYSTEM}

  ${APK_SOURCE_DIR}/Bc
  ${APK_SOURCE_DIR}/Bc/Disperse
  ${APK_SOURCE_DIR}/Bc/ESP
  ${APK_SOURCE_DIR}/Bc/Groundwater
  ${APK_SOURCE_DIR}/Bc/MPFV
  ${APK_SOURCE_DIR}/Bc/Radiation

  ${APK_SOURCE_DIR}/Cgns
  ${APK_SOURCE_DIR}/Cgns/adf
  ${APK_SOURCE_DIR}/Cgns/adfh
  ${APK_SOURCE_DIR}/Cgns/f2c
  ${APK_SOURCE_DIR}/Cgns/tests
  ${APK_SOURCE_DIR}/Cgns/tools

  ${APK_SOURCE_DIR}/Common
  ${APK_SOURCE_DIR}/Common/ESP
  ${APK_SOURCE_DIR}/Common/LMMHD
  ${APK_SOURCE_DIR}/Common/Radiation

  ${APK_SOURCE_DIR}/Geometry
  ${APK_SOURCE_DIR}/Geometry/Radiation

  ${APK_SOURCE_DIR}/Initialize
  ${APK_SOURCE_DIR}/Initialize/Combustion
  ${APK_SOURCE_DIR}/Initialize/Disperse
  ${APK_SOURCE_DIR}/Initialize/ESP
  ${APK_SOURCE_DIR}/Initialize/Groundwater
  ${APK_SOURCE_DIR}/Initialize/LMMHD
  ${APK_SOURCE_DIR}/Initialize/MPFV
  ${APK_SOURCE_DIR}/Initialize/Radiation

  ${APK_SOURCE_DIR}/Io
  ${APK_SOURCE_DIR}/Io/Disperse
  ${APK_SOURCE_DIR}/Io/MPFV

  ${APK_SOURCE_DIR}/Solver
  ${APK_SOURCE_DIR}/Solver/CgMethod
  ${APK_SOURCE_DIR}/Solver/Combustion
  ${APK_SOURCE_DIR}/Solver/Disperse
  ${APK_SOURCE_DIR}/Solver/ESP
  ${APK_SOURCE_DIR}/Solver/Groundwater
  ${APK_SOURCE_DIR}/Solver/LMMHD
  ${APK_SOURCE_DIR}/Solver/MPFV
  ${APK_SOURCE_DIR}/Solver/Optimial_Functions
  ${APK_SOURCE_DIR}/Solver/ParticleTracking
  ${APK_SOURCE_DIR}/Solver/Radiation
  ${APK_SOURCE_DIR}/Solver/Solidification
  ${APK_SOURCE_DIR}/Solver/Turbulence

  ${APK_SOURCE_DIR}/User
  ${APK_SOURCE_DIR}/User/Radiation
)

#-----------------------------------------------------------------------------
# Create the Executable
SET( EXE_SRCS
APKMain.c
)


SET( Bc_SRCS_00
Bc/APKBc_input.c
Bc/APKboundary_d_t_s.c
Bc/APKboundary_d_t_sTurb.c
Bc/APKCall_user_func.c
Bc/APKcheck_bc_type.c
Bc/APKcheck_function.c
Bc/APKInputBoundarycondition.c
Bc/APKPoint_to_Coord.c
Bc/APKprint_tagname.c
Bc/APKTimeDepend.c
Bc/APKupdateBCslabsVarCoeff.c
Bc/APKupdate_boundary.c
Bc/APKUser_Enter.c
Bc/APKWallFunctionDirichlet.c
Bc/APKWallFunctionDirichletConjugate.c
Bc/APKfind_name.c

Bc/ESP/APKWireCurrentDensity.c
Bc/ESP/APKRefreshWireCharge.c
Bc/ESP/APKGetBoundaryFaces.c
Bc/ESP/APKUpdate_Boundary_Charge.c
Bc/ESP/APKImplementChargeWireBC.c
Bc/ESP/APKWriteWireCharge.c

Bc/MPFV/APKupdate_boundaryDMPFV.c
Bc/MPFV/APKupdate_boundaryTurbDMPFV.c
Bc/Disperse/APKupdate_boundaryDisperse.c

Bc/Groundwater/APKupdateBCslabsVar_Coeff_PorousMedia.c
Bc/Groundwater/APKUpdateBCslabsKunsat.c

Bc/Radiation/APKInputRadiationBoundarycondition.c
Bc/Radiation/APKpoint_coordRadiation.c
Bc/Radiation/APKboundary_d_t_sRadiation.c
Bc/Radiation/APKboundary_d_t_sRadsymmetry.c
Bc/Radiation/APKboundary_d_t_sRadinsulation.c
Bc/Radiation/APKboundary_d_t_sRobin.c
Bc/Radiation/APKboundary_d_t_sInsulation.c
)


SET(Cgns_SRCS_00
Cgns/adf_cond.c
Cgns/adf_ftoc.c
Cgns/cg_ftoc.c
Cgns/cg_malloc.c
Cgns/cgns_error.c
Cgns/cgns_internals.c
Cgns/cgnslib.c
Cgns/adf/ADF_fortran_2_c.c
Cgns/adf/ADF_interface.c
Cgns/adf/ADF_internals.c
)

SET( Common_SRCS_00
Common/APKCalcDiffusionCoefCommon.c
Common/APKCalculateTimeStep.c
Common/APKCallMakeSlabPtrArray.c
Common/APKCountCells.c
Common/APKDiffusionCommon.c
Common/APKDuplicateSlabTemperatures.c
Common/APKfile_refresh.c
Common/APKFillerBoundary.c
Common/APKFillerConjugateBoundary.c
Common/APKFillerConjugateInterface.c
Common/APKFillerConjugateInterfaceTurb.c
Common/APKFillerIfaceConjugatePressure.c
Common/APKFillerIfaceConjugateTemperatureTurb.c
Common/APKFillerIfaceConjugateVariable.c
Common/APKFillerInterface.c
Common/APKFillerInterfaceVariables.c
Common/APKFindCellCentreArray.c
Common/APKFindConjuIfaceNeighbor.c
Common/APKFind_Direction.c
Common/APKFindInterfaceBCPtr.c
Common/APKFindInterfaceBCPtrMod.c
Common/APKFindmaXConv.c
Common/APKFindmaXDiff.c
Common/APKFindmiN.c
Common/APKFindSlabPtrBC.c
Common/APKFindSlabPtrIface.c
Common/APKfindZone2Variable.c
Common/APKfindZone2Volumes.c
Common/APKHomogeneousDirichletIC.c
Common/APKHomogeneousNeumannIC.c
Common/APKinputFunc.c
Common/APKMakeDomainArray.c
Common/APKMakeSlabPtrArray.c
Common/APKModifiedDiffCoeff.c
Common/APKNonOrthogonalHomogeneousNeumannIC.c
Common/APKNonOrthogonalNeumannIC.c
Common/APKprintfluxes.c
Common/APKprintmassflux.c
Common/APKprintslabs.c
Common/APKprintsolver.c
Common/APKprintu4bc.c
Common/APKTerminate.c
Common/APKTransformCellIndicesCommon.c
Common/APKfind_nv.c

Common/APKMFillerInterfaceVariables.c
Common/APKMFillerInterfaceVariables2d.c
Common/APKMfindZone2Variable.c
Common/APKMfindZone2Variable2d.c
Common/APKpmonitors.c
Common/APKprint_monitors.c
Common/APKfind_zone_face_info.c
Common/APKMakeSlabPtrArray_com.c
Common/APKCalcIFACE_data.c

Common/ESP/APKinputFuncElectric.c

Common/LMMHD/APKFillerIfaceConjugatePotentialMHD.c

Common/Radiation/APKFillerInterfaceRadiationVariables.c
Common/Radiation/APKCalGradRad.c
)

SET( Geometry_SRCS_00
Geometry/APKAdjoiningZonePointerXYZ.c
Geometry/APKCalculateGeometry.c
Geometry/APKCalSurf.c
Geometry/APKCheckOutward.c
Geometry/APKDerive_Transformation_Matrix.c
Geometry/APKface_center.c
Geometry/APKFindBeta.c
Geometry/APKFindCellCentre.c
Geometry/APKFindCellTransformRule.c
Geometry/APKFindFaceCenterBC.c
Geometry/APKFindZonePtr.c
Geometry/APKGetSurfacePoints.c
Geometry/APKTransformCellIndices1.c
Geometry/APKUpdateInterfaceStruct.c

Geometry/Radiation/APKCalculateDmF.c
Geometry/Radiation/APKCalculateNormalVector.c
)

SET( Initialize_SRCS_00
Initialize/APKinitialize.c
Initialize/APKinitialize_Effviscosity.c
Initialize/APKinitialize_Effviscosityofphase.c
Initialize/APKinitialize_Dummyeffviscosityofphase.c
Initialize/APKinitializeslab.c
Initialize/APKinitslabmag.c
Initialize/APKinitializeSlabVarCoeff.c
Initialize/APKinitializeSolidLiquid.c
Initialize/APKinitializeTurbVariables.c
Initialize/APKinitializeVarCoefficient.c
Initialize/APKinitialize_variable.c
Initialize/APKinitialize_T.c
Initialize/APKinitialize_enthalpy.c
Initialize/APKinitialize_TwoPhase.c
Initialize/APKinitialize_monitors.c
Initialize/Combustion/APKinitializeCombParameters.c
Initialize/Combustion/APKinitializeCombVariables.c
Initialize/Combustion/APKInitializeNonSolveCombVariables.c

Initialize/ESP/APKinitializeElectric.c

Initialize/LMMHD/APKinitializeLMMHD.c

Initialize/MPFV/APKinitializeMPFV.c
Initialize/MPFV/APKinitializeVarMPFV.c

Initialize/Disperse/APKinitializeDisperse.c
Initialize/Disperse/APKinitializeVarDisperse.c

Initialize/Groundwater/APKinitializeGroundwater.c
Initialize/Groundwater/APKinitializeDispersionCoeff.c
Initialize/Groundwater/APKinitializeSlabDispersionCoeff.c
Initialize/Groundwater/APKupdateBCDeltaInverse.c

Initialize/Radiation/APKinitializeRadiation.c
Initialize/Radiation/APKinitializeRadiationSlab.c
Initialize/Radiation/APKinitializeRadiationVarCoefficient.c
Initialize/Radiation/APKinitializeRadiationSlabVarCoeff.c
Initialize/Radiation/APKinitializeRadiation_variable.c

)

SET( Io_SRCS_00
Io/APKArrayDuplicator.c
Io/APKCGNSErrorCheck.c
Io/APKReadCgnsFile.c
Io/APKReadCgnsFileRestart.c
Io/APKReadInitialData.c
Io/APKWriteCgnsFile.c
Io/APKWriteCgnsFile_vertex.c

Io/APKVariable_AT_Vertex.c
Io/APKRead_Solution.c
Io/APKWrite_Solution.c
#Io/APKWrite_gnuplot_style_solution.c

Io/MPFV/APKWriteCgnsFileMPFV.c
Io/Disperse/APKWriteCgnsFileDisperse.c
)

SET( Solver_SRCS_00
Solver/APKCalcCoefficients.c
Solver/APKCalcCoefficients2d.c
Solver/APKCalcCoefficientsCorner.c
Solver/APKCalcCoefficientsCorner2d.c
Solver/APKCalcCoefficientsForVariables.c
Solver/APKCalcCoefficientsCornerForVariables.c
Solver/APKCalcCoefficients2dForVariables.c
Solver/APKCalcCoefficientsCorner2dForVariables.c
Solver/APKCalcCoefficientsForVariablesVarDiffCoeff.c
Solver/APKCalcCoefficientsCornerForVariablesVarDiffCoeff.c
Solver/APKCalcCoefficients2dForVariablesVarDiffCoeff.c
Solver/APKCalcCoefficientsCorner2dForVariablesVarDiffCoeff.c
Solver/APKCalcDiffusionCoef.c
Solver/APKCalcExtraTermsOfVarDensity.c
Solver/APKCalcExtraTermsOfVarDensity2d.c
Solver/APKCalcExtraTermsOfVarDensityMassflux.c
Solver/APKCalcExtraTermsOfVarDensityMassflux2d.c
Solver/APKCalcGradVar.c
Solver/APKCalcNewMassFlux2d.c
Solver/APKCalcNewMassFlux.c
Solver/APKCalcNewMassFlux_restart.c
Solver/APKCalcNewVolumeFlux.c
Solver/APKCalcNewVelocity.c
Solver/APKCalcNewVolumefluxRestart.c
Solver/APKCalcNusseltNo.c
Solver/APKCalcPressureGrad.c
Solver/APKCalcSource.c
Solver/APKCalcShearStress.c
Solver/APKCalculateErrorPressure.c
Solver/APKCalcvarCoeffAndPrevCoeff.c
Solver/APKCalcvarCoeffAndPrevCoeffConjugate.c
Solver/APKCalcVariable.c
Solver/APKCalcVariableConjugate.c
Solver/APKCalcVariableOneIter.c
Solver/APKCalcVariableOneIterConjugate.c
Solver/APKCalcVorticity.c
Solver/APKCallCalcCoefficients.c
Solver/APKCallCalcCoefficientsForVariables.c
Solver/APKCDcoeffCalculation.c
Solver/APKCDcoeffMemAlloc.c
Solver/APKCNterm.c
Solver/APKCoarseGridCoefficients.c
Solver/APKCoarseGridCoefficientsDiffusion.c

Solver/APKconvfluxC2D.c
Solver/APKconvfluxC3D.c
Solver/APKConvfluxExt2D.c
Solver/APKConvfluxExt3D.c
Solver/APKConvfluxExtCorner2D.c
Solver/APKConvfluxExtCorner3D.c
Solver/APKConvfluxExtCother2D.c
Solver/APKConvfluxExtCother3D.c
Solver/APKConvfluxExtQ2D.c
Solver/APKConvfluxExtQ3D.c
Solver/APKConvfluxIntC2D.c
Solver/APKConvfluxIntC3D.c
Solver/APKConvfluxIntQ2D.c
Solver/APKConvfluxIntQ3D.c
Solver/APKconvfluxQ2D.c
Solver/APKconvfluxQ3D.c

Solver/APKDiffusionCoeffCalculation.c
Solver/APKDiffusionCoeffCalculationForVariables.c
Solver/APKDiffusion2d.c
Solver/APKDiffusion2dWithConstDiffCoeff.c
Solver/APKDiffusion2dWithVaryingDiffCoeff.c
Solver/APKDiffusion3dWithConstDiffCoeff.c
Solver/APKDiffusion3dWithVaryingDiffCoeff.c
Solver/APKDiffusion.c
Solver/APKEvaluateEnergysource.c
Solver/APKEvaluateVarCoefficients.c
Solver/APKExtractMassVelocities.c
Solver/APKExtractVariableCoefficients.c
Solver/APKExtractVariables.c
Solver/APKExtractVolume.c
Solver/APKFillNewsolToPrevsol.c
Solver/APKfillnewtoprev.c
Solver/APKfind_convcri.c
Solver/APKfind_GSRelaxation.c
Solver/APKFindDiffusionCoeffSlabPtrBC.c
Solver/APKFindDiffusionCoeffSlabPtrIface.c
Solver/APKFindgradient.c
Solver/APKFindMassFluxDifference.c
Solver/APKfindresidual.c
Solver/APKfindresidualCD.c
Solver/APKfindresidualCD2d.c
Solver/APKfindresidualQVK.c
Solver/APKfindresidualQVK2d.c
Solver/APKFindVolumeFluxDifference.c
Solver/APKFineGridCoefficients1.c
Solver/APKFineGridCoefficients2.c
Solver/APKFineGridDefQVKCoefficients1.c
Solver/APKFineGridDefQVKCoefficients2.c
Solver/APKFineGridDefCDCoefficients1.c
Solver/APKFineGridDefCDCoefficients2.c
Solver/APKGridParameters1.c
Solver/APKGridParameters2.c
Solver/APKMakeDiffusionCoeffDomain.c
Solver/APKMakeDiffusionCoeffFlag.c
Solver/APKMakeDiffusionCoeffSlabPtrArray.c
Solver/APKmassflux2d.c
Solver/APKmassflux.c
Solver/APKMassFluxConv.c
Solver/APKmassfluxstar2d.c
Solver/APKmassfluxstar.c
Solver/APKmassVelocity.c
Solver/APKmassVelocityOneIter.c
Solver/APKMatrixMult.c
Solver/APKMatrixMultCorner.c
Solver/APKPressureForce.c
Solver/APKPressureIterArtificial.c
Solver/APKPressureIter.c
Solver/APKRefreshMassVelocities.c
Solver/APKRefreshSurfaceCells.c
Solver/APKRefreshSurfaceVolumes.c
Solver/APKRefreshVariableCoefficients.c
Solver/APKsavemassflux.c
Solver/APKSaveTimeStepData.c
Solver/APKsavevolumeflux.c
Solver/APKShearStress.c
Solver/APKSolutions2d.c
Solver/APKSolutions.c
Solver/APKsolverCoupled2d.c
Solver/APKsolverCoupled.c
Solver/APKsolverSegregated.c
Solver/APKsolverSemiCoupled2d.c
Solver/APKsolverSemiCoupled.c
Solver/APKTDMA.c

Solver/APKCorrection.c
Solver/APKDummyVariableRestrict.c
Solver/APKDummyVariableRestrict2d.c
Solver/APKDummyVariableProlongate.c
Solver/APKDummyVariableProlongate2d.c
Solver/APKFillVariable2Dummy.c
Solver/APKFillVariable2Dummy2d.c
Solver/APKInterZoneVariableTransfer.c
Solver/APKInterZoneVariableTransfer2d.c
Solver/APKMultigrid.c
Solver/APKMultigrid2d.c
Solver/APKMultigridExtractVariables.c
Solver/APKPressureIterMultigrid.c
Solver/APKProlongation.c
Solver/APKProlongation2d.c
Solver/APKResidFind0.c
Solver/APKResidFind2d0.c
Solver/APKResidFind.c
Solver/APKResidFind2d.c
Solver/APKRestriction.c
Solver/APKSmoothing.c
Solver/APKSmoothing2d.c
Solver/APKCalcConjugate.c
Solver/APKCalcSteady_vel.c
Solver/APKConv_Coeff_Calculation.c
Solver/APKConv_Diff_Assembly.c

Solver/Optimial_Functions/APKDiffusion2dWithConstDiffCoeff_opti.c
Solver/Optimial_Functions/APKDiffusion3dWithConstDiffCoeff_opti.c
Solver/Optimial_Functions/APKPressure_Grad_Calculation.c
Solver/Optimial_Functions/APKCalcPressureGrad_Common.c
Solver/Optimial_Functions/APKCalcSource_opti.c
Solver/Optimial_Functions/APKFineGridCoefficients1_solid.c
Solver/Optimial_Functions/APKCalcCoefficients_opti.c
Solver/Optimial_Functions/APKCalcCoefficientsCorner_opti.c
Solver/Optimial_Functions/APKCalcCoefficients2d_opti.c
Solver/Optimial_Functions/APKCalcCoefficientsCorner2d_opti.c
Solver/Optimial_Functions/APKCallCalcCoefficients_opti.c
Solver/Optimial_Functions/APKFind_Mushy_Source.c

Solver/ESP/APKSolverElectric.c
Solver/ESP/APKSolverElectricSemiCoupled2d.c
Solver/ESP/APKSolverElectricSemiCoupled.c
Solver/ESP/APKCalcVariablePotential.c
Solver/ESP/APKCalcVariableCharge.c
Solver/ESP/APKElectricFieldXYZ.c
Solver/ESP/APKElectricFluxXYZ.c
Solver/ESP/APKPotentialIter.c
Solver/ESP/APKUpdatePotential2d.c
Solver/ESP/APKCoeffcalculation_Potential.c
Solver/ESP/APKCoeffcalculation_ESP.c
Solver/ESP/APKFree_Coefficients_ESP.c
Solver/ESP/APKESP_Source_Calculation.c

Solver/LMMHD/APKFindgradientPhiMHD.c
Solver/LMMHD/APKPotentialMHDIter.c
Solver/LMMHD/APKUpdatePotentialMHD.c
Solver/LMMHD/APKUpdatePotentialMHD2d.c
Solver/LMMHD/APKPotentialMHDpoissonRHS.c
Solver/LMMHD/APKCalcCurrentDensity.c
Solver/LMMHD/APKCalcUcrossB.c
Solver/LMMHD/APKsolverLMMHD.c
Solver/LMMHD/APKInducedMagneticSource.c
Solver/LMMHD/APKInducedMagneticSource2d.c
Solver/LMMHD/APKUpdateinducedB.c
Solver/LMMHD/APKUpdateinducedB2d.c
Solver/LMMHD/APKCalcNewMagneticflux.c
Solver/LMMHD/APKMagneticflux.c
Solver/LMMHD/APKCalcPsiSource.c
Solver/LMMHD/APKsavemagneticflux.c
Solver/LMMHD/APKInduced_Source_Calculation.c
Solver/LMMHD/APKLMMHD_Source_Calculation.c
Solver/LMMHD/APKLMMHDmass_Source_Calculation.c
Solver/LMMHD/APKInducedLMMHD_Source_Calculation.c
Solver/LMMHD/APKInducedMomentumSource2d.c
Solver/LMMHD/APKInducedMomentumSource.c
Solver/LMMHD/APKCalcKhiGrad_Common.c



Solver/MPFV/APKDiffusion3dWithConstDiffCoeff_optiMPFV.c
Solver/MPFV/APKDiffusion2dWithConstDiffCoeff_optiMPFV.c
Solver/MPFV/APKCalcCoefficientsCorner2dForVariablesVarDiffCoeffMPFV.c
Solver/MPFV/APKCalcCoefficients2dForVariablesVarDiffCoeffMPFV.c
Solver/MPFV/APKCalcCoefficientsCornerForVariablesVarDiffCoeffMPFV.c
Solver/MPFV/APKCalcCoefficientsForVariablesVarDiffCoeffMPFV.c
Solver/MPFV/APKUpdateVarMPFVDiffusion.c
Solver/MPFV/APKDragMPFV.c
Solver/MPFV/APKLiftMPFV.c
Solver/MPFV/APKMPFV.c
Solver/MPFV/APKvolumefluxMPFV.c

Solver/Groundwater/APKCalDispersionCoeff.c
Solver/Groundwater/APKCalGradScalerXY.c
Solver/Groundwater/APKCalVelocity.c
Solver/Groundwater/APKdatum2d.c
Solver/Groundwater/APKdatumUnsat2d.c
Solver/Groundwater/APKDisperisonflux2d.c
Solver/Groundwater/APKPressureHeadIter.c
Solver/Groundwater/APKUpdateWaterContent.c
Solver/Groundwater/APKUpdatePressureHead2d.c
Solver/Groundwater/APKUpdatePressureHead3D.c
Solver/Groundwater/APKExtract_Var_Coeff_PorousMedia.c
Solver/Groundwater/APKUpdateScaler.c

Solver/Turbulence/APKCalcNuDiffCoef.c
Solver/Turbulence/APKLimitVariable.c
Solver/Turbulence/APKSolveTurbVariables.c
Solver/Turbulence/APKSpalartAllmarasSource.c
Solver/Turbulence/APKTURBCalcC1epsi.c
Solver/Turbulence/APKTURBCalcDoubleDerivative2d.c
Solver/Turbulence/APKTURBCalcDoubleDerivative3d.c
Solver/Turbulence/APKTURBCalcEpsilonWallTerm.c
Solver/Turbulence/APKTURBCalcGamma.c
Solver/Turbulence/APKTURBCalcGradSquare.c
Solver/Turbulence/APKTURBCalcKWallTerm.c
Solver/Turbulence/APKTURBCalcmoduPt.c
Solver/Turbulence/APKTURBCalcMuEff.c
Solver/Turbulence/APKTURBCalcTurbconductivity.c
Solver/Turbulence/APKTURBCalcTurbulentProduction.c
Solver/Turbulence/APKTURBCalcTurbViscosity.c
Solver/Turbulence/APKTURBCalcUtau.c
Solver/Turbulence/APKTURBCalcVariableOneIterEps.c
Solver/Turbulence/APKTURBCalcVariableOneIterK.c
Solver/Turbulence/APKTURBCalcVariableOneIterOmega.c
Solver/Turbulence/APKTURBCalcVariableOneIterOmegaSST.c
Solver/Turbulence/APKTURBChangeIndices.c
Solver/Turbulence/APKTURBCrossDiffusion.c
Solver/Turbulence/APKTURBDistance2Wall.c
Solver/Turbulence/APKTURBfreeCrossDiff.c
Solver/Turbulence/APKTURBfreeGamma.c
Solver/Turbulence/APKTURBgrad2D.c
Solver/Turbulence/APKTURBgrad3D.c
Solver/Turbulence/APKTURBkturb.c
Solver/Turbulence/APKTURBKgradient2D.c
Solver/Turbulence/APKTURBKgradient3D.c
Solver/Turbulence/APKTURBminmax.c
Solver/Turbulence/APKTURBOmega.c
Solver/Turbulence/APKTURBOmegaSST.c
Solver/Turbulence/APKTURBUpdateOldZSlabsKEpsilon.c
Solver/Turbulence/APKTURBKSource.c

Solver/APKUpdateBoundarySlabs.c
Solver/APKUpdatePressure2d.c
Solver/APKUpdatePressure2dVarDensity.c
Solver/APKUpdatePressure2dVarDensityMassFlux.c
Solver/APKUpdatePressureVarDensityMassFlux.c
Solver/APKUpdatePressureArtificial2d.c
Solver/APKUpdatePressureArtificial.c
Solver/APKUpdatePressure.c
Solver/APKUpdatePressureVarDensity.c
Solver/APKUpdateVar.c
Solver/APKUpdateVarCD.c
Solver/APKUpdateVarCD2d.c
Solver/APKUpdateVarQVK.c
Solver/APKUpdateVarQVK2d.c
Solver/APKUpdateVarConjugate.c
Solver/APKUpdateVarConjugate2d.c
Solver/APKUpdateOldZSlabs.c
Solver/APKVarMassFluxConv.c
Solver/APKvelocitymassBc.c
Solver/APKvolumeflux.c
Solver/APKvolumefluxstar.c
Solver/Solidification/APKupdate_dele_iteration.c
Solver/Solidification/APKprint_enthalpy.c
Solver/Solidification/APKSolid_fraction.c
Solver/Combustion/APKSolveCombustionVariables.c
Solver/Combustion/APKComputeCombustionVariables.c
Solver/APKfillprevdensity.c

Solver/Radiation/APKRadiationSolverFiniteVolume3D.c
Solver/Radiation/APKRadiationSolverFiniteVolume2D.c
Solver/Radiation/APKFillNewRadsolToPrevRadsol.c
Solver/Radiation/APKCalculateIntensity.c
Solver/Radiation/APKSolution3DConstRefIndex.c
Solver/Radiation/APKSolution2DConstRefIndex.c
Solver/Radiation/APKRadiationSolutionInt3D.c
Solver/Radiation/APKIntExtremeFace3D.c
Solver/Radiation/APKIntExtremeCell3D.c
Solver/Radiation/APKRadiationSolutionInt2D.c
Solver/Radiation/APKIntExtremeFace2D.c
Solver/Radiation/APKIntExtremeCell2D.c
Solver/Radiation/APKFillOutgoingZoneBCIntensity.c
Solver/Radiation/APKRadiationSolverWithFlow3D.c
Solver/Radiation/APKRadiationSolverWithFlow2D.c

Solver/APKCalcVolumeFlux.c
Solver/Disperse/APKvolumefluxDisperse.c
Solver/Disperse/APKDrag.c
Solver/Disperse/APKCalcVarDropDia.c
Solver/Disperse/APKEvaporate.c
Solver/Disperse/APKCalcEvapoConst.c
Solver/Disperse/APKFindMassFluxDisperseDifference.c
Solver/Disperse/APKCheckConservation.c
Solver/Disperse/APKUpdateVarDisperseDiffusion.c
Solver/Disperse/APKCalcGradVar_Disperse.c
Solver/Disperse/APKFindgradient_Disperse.c
Solver/Disperse/APKCalcGradVar_Volfrac.c
Solver/Disperse/APKCalcFaceValue.c
Solver/Disperse/APKCalcVolfracFaceValue.c
Solver/Disperse/APKCalcRhoFaceValue.c
Solver/Disperse/APKCalcNewVolumeFluxDisperse.c
Solver/Disperse/APKDisperseVolumeFluxStar.c
Solver/Disperse/APKFindDisperseVolumeFluxDifference.c

Solver/CgMethod/APKConj_PressureIter.c
Solver/CgMethod/APKCalcInitPressureResidue2d.c
Solver/CgMethod/APKCalcLHS.c
Solver/CgMethod/APKConj_CalcConjugate.c
Solver/CgMethod/APKCalcInitConjVarResidue2d.c
Solver/CgMethod/APKCalcInitConjVarResidue.c

Solver/CgMethod/APKCalcInitPressureResidue.c
Solver/CgMethod/APKCalcLHS3d.c
Solver/CgMethod/APKModiPressureIter.c
Solver/CgMethod/APKModiCalcConjugate.c
Solver/CgMethod/APKCalcCorrection.c
Solver/CgMethod/APKCalcCorrection3d.c
Solver/CgMethod/APKCalcInitPressureResidue2dVarDensityMassFlux.c
Solver/CgMethod/APKCalcLHSVarDensity.c

Solver/ParticleTracking/APKFindMinMaxCoordinatesOfZones.c
Solver/ParticleTracking/APKParticleTracking.c
Solver/ParticleTracking/APKGlobalParticleTracking.c
Solver/ParticleTracking/APKQuickParticleTracking.c
Solver/ParticleTracking/APKInterpolateVelocity.c
Solver/ParticleTracking/APKInterpolateEfield.c
Solver/ParticleTracking/APKCalcNewParticleVelocity.c
Solver/ParticleTracking/APKInputParticle.c

)

SET( User_SRCS_00
User/APKBC_userDefinedVarCoefficients.c
User/APKuser_defined_funcs.c
User/APKuser_defined_funcs2.c
User/APKuserDefinedVarCoefficients.c
User/Radiation/APKChoiceScatteringfunction.c
)

SET( User_SRCS_02
User/APKBC_userDefinedVarCoefficients.c
User/APKuser_defined_funcs2.c
User/APKuserDefinedVarCoefficients.c
User/Radiation/APKChoiceScatteringfunction.c
)

SET(CMAKE_BUILD_TYPE DEBUG)
# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Debug' as none was specified.")
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

IF (WIN32)
  SET( EXE_SRCS			${EXE_SRCS}		)
  SET( Bc_SRCS_00		${Bc_SRCS_00}		)
  SET( Cgns_SRCS_00		${Cgns_SRCS_00}		)
  SET( Common_SRCS_00		${Common_SRCS_00}	)
  SET( Geometry_SRCS_00		${Geometry_SRCS_00}	)
  SET( Initialize_SRCS_00	${Initialize_SRCS_00}	)
  SET( Io_SRCS_00		${Io_SRCS_00}		)
  SET( Solver_SRCS_00		${Solver_SRCS_00}	)
  SET( User_SRCS_00		${User_SRCS_00}		)
ENDIF (WIN32)

if(COMMAND cmake_policy)
      cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

IF (WIN32)
	ADD_EXECUTABLE(Anupravaha_Win32 ${EXE_SRCS} ${Bc_SRCS_00} ${Cgns_SRCS_00} ${Common_SRCS_00} ${Geometry_SRCS_00} ${Initialize_SRCS_00} ${Io_SRCS_00} ${Solver_SRCS_00} ${User_SRCS_00} )
ENDIF (WIN32)

#-----------------------------------------------------------------------------
# Choose static or shared libraries.
OPTION(BUILD_SHARED_LIBS "Build APK with shared libraries." ON)
SET(APK_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})

#-----------------------------------------------------------------------------
# Dispatch the build into the proper subdirectories.
#SUBDIRS(Cgns Solver Bc Geometry User Initialize Io Common)
#IF( NOT (WIN32))
# ADD_SUBDIRECTORY(Bc)
# ADD_SUBDIRECTORY(Cgns)
# ADD_SUBDIRECTORY(Common)
# ADD_SUBDIRECTORY(Geometry)
# ADD_SUBDIRECTORY(Initialize)
# ADD_SUBDIRECTORY(Io)
# ADD_SUBDIRECTORY(Solver)
# ADD_SUBDIRECTORY(User)
#ENDIF(NOT (WIN32))

#IF( NOT WIN32)
# ADD_EXECUTABLE(Anupravaha_1 ${EXE_SRCS})
#ENDIF( NOT WIN32)

# This object contains everything except user defined functions
# This library along with APKuser_defined_funcs.c APKuser_defined_funcs.h APKData_Structures.h files can be distributed to the colleges.
# The user can then create the executable using "gcc  APKuser_defined_funcs.c libAPK.so -lm"
# Make sure that you set your LD_LIBRARY_PATH environment variable point to the directory that contains libAPK.so to run the executable
ADD_LIBRARY(APK ${EXE_SRCS} ${Bc_SRCS_00} ${Cgns_SRCS_00} ${Common_SRCS_00} ${Geometry_SRCS_00} ${Initialize_SRCS_00} ${Io_SRCS_00} ${Solver_SRCS_00} ${User_SRCS_02} )

IF (APPLE)
	#TARGET_LINK_LIBRARIES(Anupravaha_1 "-framework AppKit" APKCommon APKIo APKInitialize APKUser APKGeometry APKBc APKSolver APKCgns)
	TARGET_LINK_LIBRARIES(Anupravaha_1 "-framework AppKit" APKSolver APKBc APKGeometry APKUser APKInitialize APKIo APKCommon APKCgns)
ENDIF (APPLE)

IF (UNIX)
	#TARGET_LINK_LIBRARIES(Anupravaha_1 ${CMAKE_THREAD_LIBS} ${CMAKE_DL_LIBS} -lm  APKCommon APKIo APKInitialize APKUser APKGeometry APKBc APKSolver APKCgns)
	#TARGET_LINK_LIBRARIES(Anupravaha_1 ${CMAKE_THREAD_LIBS} ${CMAKE_DL_LIBS} -lm  APKSolver APKBc APKGeometry APKUser APKInitialize APKIo APKCommon APKCgns)
	ADD_EXECUTABLE(Anupravaha_1 ${EXE_SRCS} ${Bc_SRCS_00} ${Cgns_SRCS_00} ${Common_SRCS_00} ${Geometry_SRCS_00} ${Initialize_SRCS_00} ${Io_SRCS_00} ${Solver_SRCS_00} ${User_SRCS_00} )
	TARGET_LINK_LIBRARIES(Anupravaha_1 ${CMAKE_THREAD_LIBS} ${CMAKE_DL_LIBS} -lm  )
ENDIF (UNIX)


#-----------------------------------------------------------------------------
# Add compiler flags APK needs to work on this platform.  This must be
# done AFTER the call to CMAKE_EXPORT_BUILD_SETTINGS.

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${APK_REQUIRED_C_FLAGS}")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${APK_REQUIRED_CXX_FLAGS}")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${APK_REQUIRED_EXE_LINKER_FLAGS}")

# If the version of CMake was too old, complain and build nothing.
# These should be the last lines in this file.
ELSE("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" GREATER 1.5)
  MESSAGE(SEND_ERROR
          "This version of CMake is too old to build APK.  "
          "Please upgrade to CMake 1.6.")
ENDIF("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" GREATER 1.5)
