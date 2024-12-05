cmake_minimum_required(VERSION 3.19)

project(
  qsi_defect_cluster
)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
find_package(xdiag REQUIRED)

add_executable(main qsi_defect_cluster.cpp geometry.cpp)
target_link_libraries(main PRIVATE xdiag::xdiag)
