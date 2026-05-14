SET( Initialize_SRCS
APKinitialize.c
APKinitialize_Effviscosity.c
APKinitialize_Effviscosityofphase.c
APKinitialize_Dummyeffviscosityofphase.c
APKinitializeslab.c
APKinitslabmag.c
APKinitializeSlabVarCoeff.c
APKinitializeSolidLiquid.c
APKinitializeTurbVariables.c
APKinitializeVarCoefficient.c
APKinitialize_variable.c
APKinitialize_T.c
APKinitialize_enthalpy.c
APKinitialize_TwoPhase.c
APKinitialize_monitors.c
Combustion/APKinitializeCombParameters.c
Combustion/APKinitializeCombVariables.c
Combustion/APKInitializeNonSolveCombVariables.c

ESP/APKinitializeElectric.c

LMMHD/APKinitializeLMMHD.c

MPFV/APKinitializeMPFV.c
MPFV/APKinitializeVarMPFV.c

Disperse/APKinitializeDisperse.c
Disperse/APKinitializeVarDisperse.c

Groundwater/APKinitializeGroundwater.c
Groundwater/APKinitializeDispersionCoeff.c
Groundwater/APKinitializeSlabDispersionCoeff.c
Groundwater/APKupdateBCDeltaInverse.c

Radiation/APKinitializeRadiation.c
Radiation/APKinitializeRadiationSlab.c
Radiation/APKinitializeRadiationVarCoefficient.c
Radiation/APKinitializeRadiationSlabVarCoeff.c
Radiation/APKinitializeRadiation_variable.c

)
INCLUDE_DIRECTORIES(
  ${APK_SOURCE_DIR}/Initialize
  ${APK_SOURCE_DIR}/Initialize/Combustion
  ${APK_SOURCE_DIR}/Initialize/ESP
  ${APK_SOURCE_DIR}/Initialize/LMMHD
  ${APK_SOURCE_DIR}/Initialize/Groundwater
  ${APK_SOURCE_DIR}/Initialize/MPFV
  ${APK_SOURCE_DIR}/Initialize/Radiation
)


IF (WIN32)
  SET(Initialize_SRCS ${Initialize_SRCS})
ENDIF (WIN32)


ADD_LIBRARY(APKInitialize ${Initialize_SRCS})

IF(APK_APPLE_SHARED_FLAGS_NEEDED)
  SET_TARGET_PROPERTIES(APKInitialize PROPERTIES LINK_FLAGS "-init _oft_initCom")
ENDIF(APK_APPLE_SHARED_FLAGS_NEEDED)
IF (APPLE)
  TARGET_LINK_LIBRARIES(APKInitialize "-framework AppKit" APKInitialize)
ENDIF (APPLE)

IF (UNIX)
  TARGET_LINK_LIBRARIES(APKInitialize ${CMAKE_THREAD_LIBS} ${CMAKE_DL_LIBS} -lm APKInitialize)
ENDIF (UNIX)

