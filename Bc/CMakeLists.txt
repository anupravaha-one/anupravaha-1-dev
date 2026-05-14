SET( Bc_SRCS
APKBc_input.c
APKboundary_d_t_s.c
APKboundary_d_t_sTurb.c
APKCall_user_func.c
APKcheck_bc_type.c
APKcheck_function.c
APKInputBoundarycondition.c
APKPoint_to_Coord.c
APKprint_tagname.c
APKTimeDepend.c
APKupdateBCslabsVarCoeff.c
APKupdate_boundary.c
APKUser_Enter.c
APKWallFunctionDirichlet.c
APKWallFunctionDirichletConjugate.c
APKfind_name.c

ESP/APKWireCurrentDensity.c
ESP/APKRefreshWireCharge.c
ESP/APKGetBoundaryFaces.c
ESP/APKUpdate_Boundary_Charge.c
ESP/APKImplementChargeWireBC.c
ESP/APKWriteWireCharge.c


MPFV/APKupdate_boundaryDMPFV.c
MPFV/APKupdate_boundaryTurbDMPFV.c
Disperse/APKupdate_boundaryDisperse.c

Groundwater/APKupdateBCslabsVar_Coeff_PorousMedia.c
Groundwater/APKUpdateBCslabsKunsat.c

Radiation/APKInputRadiationBoundarycondition.c
Radiation/APKpoint_coordRadiation.c
Radiation/APKboundary_d_t_sRadiation.c
Radiation/APKboundary_d_t_sRadsymmetry.c
Radiation/APKboundary_d_t_sRadinsulation.c
Radiation/APKboundary_d_t_sRobin.c
Radiation/APKboundary_d_t_sInsulation.c
)


IF (WIN32)
  SET(Bc_SRCS ${Bc_SRCS})
ENDIF (WIN32)


ADD_LIBRARY(APKBc ${Bc_SRCS})

IF(APK_APPLE_SHARED_FLAGS_NEEDED)
  SET_TARGET_PROPERTIES(APKBc PROPERTIES LINK_FLAGS "-init _oft_initCom")
ENDIF(APK_APPLE_SHARED_FLAGS_NEEDED)
IF (APPLE)
  TARGET_LINK_LIBRARIES(APKBc "-framework AppKit" APKCommon)
ENDIF (APPLE)

IF (UNIX)
  TARGET_LINK_LIBRARIES(APKBc ${CMAKE_THREAD_LIBS} ${CMAKE_DL_LIBS} -lm APKCommon)
ENDIF (UNIX)

