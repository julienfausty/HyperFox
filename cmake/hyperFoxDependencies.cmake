SUBDIRLIST(hyperFoxSrcNames "${CMAKE_SOURCE_DIR}/src")

foreach(name ${hyperFoxSrcNames})
  if(NOT "${name}" STREQUAL "globals")
    set(hyperFoxLibs ${hyperFoxLibs} ${name})
  endif()
    set(hyperFoxIncludes ${hyperFoxIncludes} "${CMAKE_SOURCE_DIR}/src/${name}")
endforeach()

SUBDIRLIST(hyperFoxTestNames "${CMAKE_SOURCE_DIR}/tests/unittests")

foreach(name ${hyperFoxTestNames})
  file(GLOB
    tempTestSrcFiles
    "${CMAKE_SOURCE_DIR}/tests/unittests/${name}/*.cpp")
  set(testSrcFiles ${testSrcFiles} ${tempTestSrcFiles})
endforeach()
