cmake_minimum_required(VERSION 3.16)

if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_CURRENT_BINARY_DIR)
  message(FATAL_ERROR "Use a separate build directory, for example cmake -S . -B cmake-build")
endif()

project(ISFFT LANGUAGES Fortran)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/install" CACHE PATH "Installation prefix" FORCE)
endif()

find_package(MPI REQUIRED COMPONENTS Fortran)
if(NOT MPI_Fortran_HAVE_F90_MODULE)
  message(FATAL_ERROR "ISFFT requires the MPI Fortran module mpi")
endif()

option(ISFFT_BUILD_EXAMPLE "Build the CFD toy solver" OFF)

set(ISFFT_MODULE_DIR "${CMAKE_CURRENT_BINARY_DIR}/include")
file(MAKE_DIRECTORY "${ISFFT_MODULE_DIR}")

add_library(isfft STATIC
  src/ISFFT_Communicator.f90
  src/ISFFT_Precision.f90
  src/ISFFT_Problem_State.f90
  src/ISFFT_Discrete_Operators.f90
  src/ISFFT_Workspace.f90
  src/ISFFT_Plans.f90
  src/ISFFT_MPI_Part.f90
  src/ISFFT_Halos_Exchange.f90
  src/ISFFT_FFT_kernel.f90
  src/ISFFT_DST_kernel.f90
  src/ISFFT_DCT_kernel.f90
  src/ISFFT_Modal_Operations.f90
  src/ISFFT_Communication.f90
  src/ISFFT_Axis_Execution.f90
  src/ISFFT_Transform_Execution.f90
  src/ISFFT_Allocation.f90
  src/ISFFT_Precompute.f90
  src/ISFFT_Planning.f90
  src/ISFFT_Initialization.f90
  src/ISFFT_Redistribution.f90
  src/ISFFT_Integrated_Solver.f90
  src/ISFFT.f90
)
set_target_properties(isfft PROPERTIES
  Fortran_MODULE_DIRECTORY "${ISFFT_MODULE_DIR}"
  ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib"
)
target_include_directories(isfft PUBLIC "${ISFFT_MODULE_DIR}")
target_link_libraries(isfft PUBLIC MPI::MPI_Fortran)

if(ISFFT_BUILD_EXAMPLE)
  set(ISFFT_EXAMPLE_DIR "${CMAKE_CURRENT_BINARY_DIR}/example")
  file(MAKE_DIRECTORY "${ISFFT_EXAMPLE_DIR}")
  # menu.f90 includes all external solver sources.
  add_executable(Poisson example/menu.f90)
  set_target_properties(Poisson PROPERTIES
    OUTPUT_NAME Poisson.out
    RUNTIME_OUTPUT_DIRECTORY "${ISFFT_EXAMPLE_DIR}"
    Fortran_MODULE_DIRECTORY "${ISFFT_EXAMPLE_DIR}"
  )
  target_include_directories(Poisson PRIVATE
    "${ISFFT_EXAMPLE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/example"
  )
  target_link_libraries(Poisson PRIVATE isfft)
  configure_file(example/Poisson_solver_UI.in example/Poisson_solver_UI.in COPYONLY)
  configure_file(example/Post_Poisson.m example/Post_Poisson.m COPYONLY)
endif()

install(TARGETS isfft ARCHIVE DESTINATION lib)
install(DIRECTORY "${ISFFT_MODULE_DIR}/" DESTINATION include FILES_MATCHING PATTERN "*.mod")
