SUBDIRLIST(hyperFoxSrcNames "${CMAKE_SOURCE_DIR}/src")

foreach(name ${hyperFoxSrcNames})
  set(hyperFoxLibs ${hyperFoxLibs} ${name})
  set(hyperFoxIncludes ${hyperFoxIncludes} "${CMAKE_SOURCE_DIR}/src/${name}")
endforeach()

SUBDIRLIST(hyperFoxTestNames "${CMAKE_SOURCE_DIR}/tests")

foreach(name ${hyperFoxTestNames})
  file(GLOB
    tempTestSrcFiles
    "${CMAKE_SOURCE_DIR}/tests/${name}/*.cpp")
  set(testSrcFiles ${testSrcFiles} ${tempTestSrcFiles})
endforeach()
