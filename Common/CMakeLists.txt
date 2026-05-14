SET( Common_SRCS
APKCalcDiffusionCoefCommon.c
APKCalculateTimeStep.c
APKCallMakeSlabPtrArray.c
APKCountCells.c
APKDiffusionCommon.c
APKDuplicateSlabTemperatures.c
APKfile_refresh.c
APKFillerBoundary.c
APKFillerConjugateBoundary.c
APKFillerConjugateInterface.c
APKFillerConjugateInterfaceTurb.c
APKFillerIfaceConjugatePressure.c
APKFillerIfaceConjugateTemperatureTurb.c
APKFillerIfaceConjugateVariable.c
APKFillerInterface.c
APKFillerInterfaceVariables.c
APKFindCellCentreArray.c
APKFindConjuIfaceNeighbor.c
APKFind_Direction.c
APKFindInterfaceBCPtr.c
APKFindInterfaceBCPtrMod.c
APKFindmaXConv.c
APKFindmaXDiff.c
APKFindmiN.c
APKFindSlabPtrBC.c
APKFindSlabPtrIface.c
APKfindZone2Variable.c
APKfindZone2Volumes.c
APKHomogeneousDirichletIC.c
APKHomogeneousNeumannIC.c
APKinputFunc.c
APKMakeDomainArray.c
APKMakeSlabPtrArray.c
APKModifiedDiffCoeff.c
APKNonOrthogonalHomogeneousNeumannIC.c
APKNonOrthogonalNeumannIC.c
APKprintfluxes.c
APKprintmassflux.c
APKprintslabs.c
APKprintsolver.c
APKprintu4bc.c
APKTerminate.c
APKTransformCellIndicesCommon.c
APKfind_nv.c

APKMFillerInterfaceVariables.c
APKMFillerInterfaceVariables2d.c
APKMfindZone2Variable.c
APKMfindZone2Variable2d.c
APKpmonitors.c
APKprint_monitors.c
APKfind_zone_face_info.c
APKMakeSlabPtrArray_com.c
APKCalcIFACE_data.c

ESP/APKinputFuncElectric.c

LMMHD/APKFillerIfaceConjugatePotentialMHD.c

Radiation/APKFillerInterfaceRadiationVariables.c
Radiation/APKCalGradRad.c
)



IF (WIN32)
  SET(Common_SRCS ${Common_SRCS})
ENDIF (WIN32)


ADD_LIBRARY(APKCommon ${Common_SRCS})

IF(APK_APPLE_SHARED_FLAGS_NEEDED)
  SET_TARGET_PROPERTIES(APKCommon PROPERTIES LINK_FLAGS "-init _oft_initCom")
ENDIF(APK_APPLE_SHARED_FLAGS_NEEDED)
IF (APPLE)
  TARGET_LINK_LIBRARIES(APKCommon "-framework AppKit" APKCommon)
ENDIF (APPLE)

IF (UNIX)
  TARGET_LINK_LIBRARIES(APKCommon ${CMAKE_THREAD_LIBS} ${CMAKE_DL_LIBS} -lm APKCommon)
ENDIF (UNIX)

