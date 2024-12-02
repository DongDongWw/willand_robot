cmake_minimum_required(VERSION 3.0.2)
project(traj_tracking)
find_package(catkin REQUIRED COMPONENTS
  roscpp
  std_msgs
  geometry_msgs
  nav_msgs
)
find_package(Eigen3 REQUIRED)
add_subdirectory(${PROJECT_SOURCE_DIR}/third_party/osqp)
add_subdirectory(${PROJECT_SOURCE_DIR}/third_party/osqp-eigen)

# protobuffer
find_package(Protobuf 3 REQUIRED)
set(PROTO_DATA_PATH "/tmp/ros/proto/traj_tracking")

SET(PROTO_DIR ${PROJECT_SOURCE_DIR}/proto/)
SET(PB_CXX_DIR ${PROJECT_SOURCE_DIR}/generated)
SET(PB_PYTHON_DIR ${PROJECT_SOURCE_DIR}/script/generated)

file(REMOVE_RECURSE ${PROTO_DATA_PATH})

# file(REMOVE_RECURSE ${PB_CXX_DIR})
# file(REMOVE_RECURSE ${PB_PYTHON_DIR})
file(MAKE_DIRECTORY ${PROTO_DATA_PATH})

# file(MAKE_DIRECTORY ${PB_CXX_DIR})
# file(MAKE_DIRECTORY ${PB_PYTHON_DIR})
message(STATUS "Protobuf compiler: ${PROTOBUF_PROTOC_EXECUTABLE}")
message(STATUS "[PROTOBUF_INCLUDE_DIRS]" ${Protobuf_INCLUDE_DIRS})
message("[PROTO_DIR]" ${PROTO_DIR})

LIST(APPEND PROTO_FLAGS -I${PROTO_DIR})
file(GLOB_RECURSE MSG_PROTOS ${PROTO_DIR}/*.proto)

foreach(msg ${MSG_PROTOS})
  get_filename_component(FIL_NAME ${msg} NAME_WE)
  get_filename_component(FILE_PATH_ABS ${msg} ABSOLUTE)
  string(REPLACE "${PROTO_DIR}" "" FILE_PATH_REL ${FILE_PATH_ABS})
  string(REPLACE ".proto" ".pb.cc" FILE_NAME_CC ${FILE_PATH_REL})
  string(REPLACE ".proto" ".pb.h" FILE_NAME_H ${FILE_PATH_REL})
  message("[pb file absolute path] " ${PB_CXX_DIR}/${FILE_NAME_H})
  message("[pb file absolute path] " ${PB_CXX_DIR}/${FILE_NAME_CC})
  add_custom_command(
    OUTPUT
    "${PB_CXX_DIR}/${FILE_NAME_H}"
    "${PB_CXX_DIR}/${FILE_NAME_CC}"
    COMMAND ${PROTOBUF_PROTOC_EXECUTABLE}
    ARGS
    --cpp_out ${PB_CXX_DIR}
    --python_out ${PB_PYTHON_DIR}
    -I ${PROTO_DIR}
    ${msg}
    DEPENDS ${msg}
    COMMENT "Running C++ protocol buffer compiler on ${msg}"
    VERBATIM
  )

  set(PROTO_SRC ${PROTO_SRC} ${PB_CXX_DIR}/${FILE_NAME_CC})
  set(PROTO_INC ${PROTO_INC} ${PB_CXX_DIR}/${FILE_NAME_H})
endforeach()

set_source_files_properties(${PROTO_SRC} ${PROTO_INC} PROPERTIES GENERATED TRUE)

catkin_package(
  INCLUDE_DIRS include
  LIBRARIES traj_tracking
  CATKIN_DEPENDS roscpp std_msgs geometry_msgs
)

# build target
include_directories(
  include
  ${EIGEN3_INCLUDE_DIR}
  ${catkin_INCLUDE_DIRS}
  ${PB_CXX_DIR}
)

add_executable(${PROJECT_NAME}_sim
  src/traj_tracking_sim.cpp
  src/trajectory_tracker.cpp
  ${PROTO_SRC}
)
target_link_libraries(${PROJECT_NAME}_sim
  ${catkin_LIBRARIES}
  OsqpEigen::OsqpEigen
  protobuf::libprotobuf
)
