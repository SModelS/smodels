set(resummino_SRCS
  debug.cc
  hxs.cc
  hxs_dlnm2.cc
  hxs_dpt2.cc
  kinematics.cc
  maths.cc
  params.cc
  pdf.cc
  resummation.cc
  resummino.cc
  dipoles.cc
  integration_method.cc
  tensors.cc
)



include_directories(gaugino_squark/generated)
file(GLOB_RECURSE gaugino_squark_SRC "gaugino_squark/*.cc")
file(GLOB_RECURSE gaugino_gluino_SRC "gaugino_gluino/*.cc")
file(GLOB_RECURSE gauginos_SRC "gauginos/*.cc")
file(GLOB_RECURSE leptons_SRC "leptons/*.cc")
file(GLOB_RECURSE sleptons_SRC "sleptons/*.cc")
file(GLOB_RECURSE matrix_SRC "matrix/*.cc")

include_directories(".")
include_directories(${install_dir}/external/include)

add_library( libresummino STATIC ${resummino_SRCS} ${gaugino_squark_SRC} ${gaugino_gluino_SRC} ${gauginos_SRC} ${leptons_SRC} ${sleptons_SRC} ${matrix_SRC})

if(BUILD_LHAPDF)
	add_dependencies(libresummino lhapdf)
endif()
if(BUILD_LOOPTOOLS)
	add_dependencies(libresummino lt)
endif()

#if(CI)
#	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage -fprofile-arcs -ftest-coverage -O0 -g")
#endif()


## changed the looptools configure script that the folder is always called lib
target_link_libraries(libresummino ${GSL_LIBRARIES}
	${LHAPDF_LIBRARIES} 
	${LOOPTOOLS_LIBRARIES} 
	 gfortran)

add_executable(resummino  main.cc)
target_link_libraries(resummino libresummino)

install(TARGETS resummino RUNTIME DESTINATION bin)
install(TARGETS libresummino ARCHIVE DESTINATION lib)
