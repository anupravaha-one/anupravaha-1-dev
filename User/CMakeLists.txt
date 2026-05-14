SET( User_SRCS
APKBC_userDefinedVarCoefficients.c
APKuser_defined_funcs.c
APKuser_defined_funcs2.c
APKuserDefinedVarCoefficients.c
Radiation/APKChoiceScatteringfunction.c
)



IF (WIN32)
  SET(User_SRCS ${User_SRCS})
ENDIF (WIN32)


ADD_LIBRARY(APKUser ${User_SRCS})

IF(APK_APPLE_SHARED_FLAGS_NEEDED)
  SET_TARGET_PROPERTIES(APKUser PROPERTIES LINK_FLAGS "-init _oft_initCom")
ENDIF(APK_APPLE_SHARED_FLAGS_NEEDED)
IF (APPLE)
  TARGET_LINK_LIBRARIES(APKUser "-framework AppKit" APKUser)
ENDIF (APPLE)

IF (UNIX)
  TARGET_LINK_LIBRARIES(APKUser ${CMAKE_THREAD_LIBS} ${CMAKE_DL_LIBS} -lm APKUser)
ENDIF (UNIX)

