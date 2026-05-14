SET( Geometry_SRCS
APKAdjoiningZonePointerXYZ.c
APKCalculateGeometry.c
APKCalSurf.c
APKCheckOutward.c
APKDerive_Transformation_Matrix.c
APKface_center.c
APKFindBeta.c
APKFindCellCentre.c
APKFindCellTransformRule.c
APKFindFaceCenterBC.c
APKFindZonePtr.c
APKGetSurfacePoints.c
APKTransformCellIndices1.c
APKUpdateInterfaceStruct.c

Radiation/APKCalculateDmF.c
Radiation/APKCalculateNormalVector.c
)



IF (WIN32)
  SET(Geometry_SRCS ${Geometry_SRCS})
ENDIF (WIN32)


ADD_LIBRARY(APKGeometry ${Geometry_SRCS})

IF(APK_APPLE_SHARED_FLAGS_NEEDED)
  SET_TARGET_PROPERTIES(APKGeometry PROPERTIES LINK_FLAGS "-init _oft_initCom")
ENDIF(APK_APPLE_SHARED_FLAGS_NEEDED)
IF (APPLE)
  TARGET_LINK_LIBRARIES(APKGeometry "-framework AppKit" APKGeometry)
ENDIF (APPLE)

IF (UNIX)
  TARGET_LINK_LIBRARIES(APKGeometry ${CMAKE_THREAD_LIBS} ${CMAKE_DL_LIBS} -lm APKGeometry)
ENDIF (UNIX)

