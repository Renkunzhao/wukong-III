cmake_minimum_required(VERSION 3.5)
project(rkz_wk LANGUAGES CXX)
add_definitions(-DEIGEN_STACK_ALLOCATION_LIMIT=0)                       # 不加这个mpctest报错栈溢出

find_package(raisim CONFIG REQUIRED)
find_package(Eigen3 REQUIRED HINTS ${Eigen3_HINT})
find_package(ROOT 6.26 CONFIG REQUIRED)
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
find_package(Python3 COMPONENTS NumPy)
find_package(PythonLibs REQUIRED)

message("${raisim_FOUND}")
message("${raisim_INCLUDE_DIRS}")
message("${raisim_LIBRARIES}")

message("${EIGEN3_FOUND}")
message("${EIGEN3_INCLUDE_DIRS}")
message("${EIGEN3_LIBRARIES}")

message("${ROOT_FOUND}")
message("${ROOT_INCLUDE_DIRS}")
message("${ROOT_LIBRARIES}")

message("${ROOT_CXX_FLAGS}")

message("${PYTHON3_FOUND}")
message("${PYTHON_INCLUDE_DIRS}")
message("${PYTHON_LIBRARIES}")

include_directories (   ${raisim_INCLUDE_DIRS}
                        ${EIGEN3_INCLUDE_DIRS} 
                        ${ROOT_INCLUDE_DIRS} 
                        ${PYTHON_INCLUDE_DIRS} 
                        /home/rkz/.local/lib/python3.8/site-packages/numpy/core/include
                        ${CMAKE_CURRENT_SOURCE_DIR}/include )

link_libraries( ${raisim_LIBRARIES} 
                pthread 
                ${ROOT_LIBRARIES} 
                ${PYTHON_LIBRARIES} )

#==
# Build Targets
#==
add_library(kinematics src/kinematics.cpp)
add_library(kinematics_foot src/kinematics_foot.cpp)
target_link_libraries(kinematics_foot PUBLIC kinematics)

add_executable(kinematics_test src/kinematics_test.cpp)
target_link_libraries(kinematics_test PUBLIC kinematics kinematics_foot)

add_library(adaptive src/adaptive.cpp)
target_link_libraries(adaptive PUBLIC kinematics)

add_executable(wukong src/wukong.cpp)
target_link_libraries(wukong PUBLIC kinematics adaptive kinematics_foot)

add_executable(wukong_test src/wukong_test.cpp)
target_link_libraries(wukong_test PUBLIC kinematics adaptive)

add_executable(root_tutorials src/root_tutorials.cpp)
add_executable(matplotlib src/matplotlib.cpp)