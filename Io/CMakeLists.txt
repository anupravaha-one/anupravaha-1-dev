SET( Io_SRCS
APKArrayDuplicator.c
APKCGNSErrorCheck.c
APKReadCgnsFile.c
APKReadCgnsFileRestart.c
APKReadInitialData.c
APKWriteCgnsFile.c
APKWriteCgnsFile_vertex.c

APKVariable_AT_Vertex.c
APKRead_Solution.c
APKWrite_Solution.c
#APKWrite_gnuplot_style_solution.c

MPFV/APKWriteCgnsFileMPFV.c
Disperse/APKWriteCgnsFileDisperse.c
)

INCLUDE_DIRECTORIES(
  ${APK_SOURCE_DIR}/Io
  ${APK_SOURCE_DIR}/Io/MPFV
  ${APK_SOURCE_DIR}/Io/Disperse
)


IF (WIN32)
  SET(Io_SRCS ${Io_SRCS})
ENDIF (WIN32)


ADD_LIBRARY(APKIo ${Io_SRCS})

IF(APK_APPLE_SHARED_FLAGS_NEEDED)
  SET_TARGET_PROPERTIES(APKIo PROPERTIES LINK_FLAGS "-init _oft_initCom")
ENDIF(APK_APPLE_SHARED_FLAGS_NEEDED)
IF (APPLE)
  TARGET_LINK_LIBRARIES(APKIo "-framework AppKit" APKIo)
ENDIF (APPLE)

IF (UNIX)
  TARGET_LINK_LIBRARIES(APKIo ${CMAKE_THREAD_LIBS} ${CMAKE_DL_LIBS} -lm APKIo)
ENDIF (UNIX)

