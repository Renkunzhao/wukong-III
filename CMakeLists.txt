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

add_executable(forward_kinematics src/forward_kinematics.cpp)
target_link_libraries(forward_kinematics PUBLIC kinematics)

add_executable(kinematics_test src/kinematics_test.cpp)
target_link_libraries(kinematics_test PUBLIC kinematics)

add_library(adaptive src/adaptive.cpp)
target_link_libraries(adaptive PUBLIC kinematics)

add_library(adaptive3 src/adaptive3.cpp)
target_link_libraries(adaptive3 PUBLIC kinematics)

add_executable(parameter_identify src/parameter_identify.cpp)
target_link_libraries(parameter_identify PUBLIC kinematics adaptive)

add_executable(wukong src/wukong.cpp)
target_link_libraries(wukong PUBLIC kinematics adaptive)

add_executable(wukong2 src/wukong2.cpp)
target_link_libraries(wukong2 PUBLIC kinematics adaptive)

add_executable(wukong3 src/wukong3.cpp)
target_link_libraries(wukong3 PUBLIC kinematics adaptive3)

add_executable(wukong3plot src/wukong3plot.cpp)
target_link_libraries(wukong3plot PUBLIC kinematics adaptive3)

add_executable(wukong_up src/wukong_up.cpp)
target_link_libraries(wukong_up PUBLIC kinematics adaptive3)

add_executable(wukong_VMC src/wukong_VMC.cpp)
target_link_libraries(wukong_VMC PUBLIC kinematics adaptive kinematics_foot)

add_executable(wukong_test src/wukong_test.cpp)
target_link_libraries(wukong_test PUBLIC kinematics adaptive)

add_executable(root_tutorials src/root_tutorials.cpp)
add_executable(matplotlib src/matplotlib.cpp)