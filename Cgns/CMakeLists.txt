SET(Cgns_SRCS
adf_cond.c
adf_ftoc.c
cg_ftoc.c
cg_malloc.c
cgns_error.c
cgns_internals.c
cgnslib.c
adf/ADF_fortran_2_c.c
adf/ADF_interface.c
adf/ADF_internals.c
)

INCLUDE_DIRECTORIES(
  ${APK_SOURCE_DIR}/Cgns
  ${APK_SOURCE_DIR}/Cgns/adf
  ${APK_SOURCE_DIR}/Cgns/f2c
)

IF (WIN32)
  SET(Cgns_SRCS ${Cgns_SRCS})
ENDIF (WIN32)


ADD_LIBRARY(APKCgns ${Cgns_SRCS})

IF(APK_APPLE_SHARED_FLAGS_NEEDED)
  SET_TARGET_PROPERTIES(APKCgns PROPERTIES LINK_FLAGS "-init _oft_initCom")
ENDIF(APK_APPLE_SHARED_FLAGS_NEEDED)
IF (APPLE)
  TARGET_LINK_LIBRARIES(APKCgns "-framework AppKit" APKCgns)
ENDIF (APPLE)

IF (UNIX)
  TARGET_LINK_LIBRARIES(APKCgns ${CMAKE_THREAD_LIBS} ${CMAKE_DL_LIBS} -lm APKCgns)
ENDIF (UNIX)

