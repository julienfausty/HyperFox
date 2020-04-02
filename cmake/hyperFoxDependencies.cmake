SUBDIRLIST(hyperFoxSrcNames "${CMAKE_SOURCE_DIR}/src")

foreach(name ${hyperFoxSrcNames})
  if(NOT "${name}" STREQUAL "globals")
    set(hyperFoxLibs ${hyperFoxLibs} ${name})
  endif()
    set(hyperFoxIncludes ${hyperFoxIncludes} "${CMAKE_SOURCE_DIR}/src/${name}")
endforeach()

set(hyperFoxIncludes ${hyperFoxIncludes} ${INCLUDES_OUTPUT_PATH})

SUBDIRLIST(hyperFoxTestNames "${CMAKE_SOURCE_DIR}/tests/unittests")

foreach(name ${hyperFoxTestNames})
  file(GLOB
    tempTestSrcFiles
    "${CMAKE_SOURCE_DIR}/tests/unittests/${name}/*.cpp")
  set(testSrcFiles ${testSrcFiles} ${tempTestSrcFiles})
endforeach()
