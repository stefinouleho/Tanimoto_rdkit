cmake_minimum_required( VERSION 2.8.9 )

project(RDKitSV)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall" )


set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#if(OPENVDB_CORE_SHARED)
#   set(BUILD_SHARED_LIBS ON)
#   set(Boost_USE_STATIC_LIBS OFF)
#endif()

  find_library(BOOST_IOSTREAMS_LIBRARY
    NAMES
      boost_iostreams${BOOST_LIB_SUFFIX}
    PATHS
      ${BOOST_SEARCH_DIRS}
      /usr/lib
      /usr/local/lib
      /opt/local/lib
      /sw/lib
  )
#message(${BOOST_IOSTREAMS_LIBRARY})
  find_library(BOOST_FILESYSTEM_LIBRARY
    NAMES
      boost_filesystem${BOOST_LIB_SUFFIX}
    PATHS
      ${BOOST_SEARCH_DIRS}
      /usr/lib
      /usr/local/lib
      /opt/local/lib
      /sw/lib
  )

  find_library(BOOST_SYSTEM_LIBRARY
    NAMES
      boost_system${BOOST_LIB_SUFFIX}
    PATHS
      ${BOOST_SEARCH_DIRS}
      /usr/lib
      /usr/local/lib
      /opt/local/lib
      /sw/lib
  )
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED OFF)
set(Boost_USE_STATIC_RUNTIME ON)
#find_package( Boost ${MINIMUM_BOOST_VERSION} REQUIRED COMPONENTS iostreams filesystem system)


find_library(CAIRO_LIB_DEBUG cairod)
find_library(CAIRO_LIB_RELEASE cairo)
find_library(CAIRO_GOBJECT_LIB_DEBUG cairo-gobjectd)
find_library(CAIRO_GOBJECT_LIB_RELEASE cairo-gobject)

find_path(CAIRO_INCLUDE_DIR cairo/cairo.h)

find_package( MyRDKit COMPONENTS ChemReactions FileParsers SmilesParse Depictor RDGeometryLib
  RDGeneral SubstructMatch Subgraphs MolDraw2D GraphMol DistGeometry DistGeomHelpers
  MolAlign Optimizer ForceField ForceFieldHelpers Alignment Fingerprints DataStructs
  ForceField  MolTransforms EigenSolvers )

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR})

find_package (Threads)
set(RDKit_THREAD_LIBS Threads::Threads)

set( LIBS ${RDKIT_LIBRARIES} ${BOOST_IOSTREAMS_LIBRARY} ${BOOST_FILESYSTEM_LIBRARY} ${BOOST_SYSTEM_LIBRARY} ${RDKit_THREAD_LIBS} z  )

include_directories(${RDKIT_INCLUDE_DIR})


add_executable( empreintes empreintes.cpp )
target_link_libraries( empreintes ${LIBS} )
