# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/rkz/桌面/simulation

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/rkz/桌面/simulation/build

# Include any dependencies generated for this target.
include CMakeFiles/wukong.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/wukong.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/wukong.dir/flags.make

CMakeFiles/wukong.dir/src/wukong.cpp.o: CMakeFiles/wukong.dir/flags.make
CMakeFiles/wukong.dir/src/wukong.cpp.o: ../src/wukong.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/rkz/桌面/simulation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/wukong.dir/src/wukong.cpp.o"
	/bin/x86_64-linux-gnu-g++-9  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/wukong.dir/src/wukong.cpp.o -c /home/rkz/桌面/simulation/src/wukong.cpp

CMakeFiles/wukong.dir/src/wukong.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/wukong.dir/src/wukong.cpp.i"
	/bin/x86_64-linux-gnu-g++-9 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/rkz/桌面/simulation/src/wukong.cpp > CMakeFiles/wukong.dir/src/wukong.cpp.i

CMakeFiles/wukong.dir/src/wukong.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/wukong.dir/src/wukong.cpp.s"
	/bin/x86_64-linux-gnu-g++-9 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/rkz/桌面/simulation/src/wukong.cpp -o CMakeFiles/wukong.dir/src/wukong.cpp.s

# Object files for target wukong
wukong_OBJECTS = \
"CMakeFiles/wukong.dir/src/wukong.cpp.o"

# External object files for target wukong
wukong_EXTERNAL_OBJECTS =

wukong: CMakeFiles/wukong.dir/src/wukong.cpp.o
wukong: CMakeFiles/wukong.dir/build.make
wukong: /home/rkz/software/root/lib/libCore.so
wukong: /home/rkz/software/root/lib/libImt.so
wukong: /home/rkz/software/root/lib/libRIO.so
wukong: /home/rkz/software/root/lib/libNet.so
wukong: /home/rkz/software/root/lib/libHist.so
wukong: /home/rkz/software/root/lib/libGraf.so
wukong: /home/rkz/software/root/lib/libGraf3d.so
wukong: /home/rkz/software/root/lib/libGpad.so
wukong: /home/rkz/software/root/lib/libROOTDataFrame.so
wukong: /home/rkz/software/root/lib/libTree.so
wukong: /home/rkz/software/root/lib/libTreePlayer.so
wukong: /home/rkz/software/root/lib/libRint.so
wukong: /home/rkz/software/root/lib/libPostscript.so
wukong: /home/rkz/software/root/lib/libMatrix.so
wukong: /home/rkz/software/root/lib/libPhysics.so
wukong: /home/rkz/software/root/lib/libMathCore.so
wukong: /home/rkz/software/root/lib/libThread.so
wukong: /home/rkz/software/root/lib/libMultiProc.so
wukong: /home/rkz/software/root/lib/libROOTVecOps.so
wukong: /usr/lib/x86_64-linux-gnu/libpython3.8.so
wukong: libkinematics.a
wukong: libadaptive.a
wukong: libkinematics_foot.a
wukong: libkinematics.a
wukong: libadaptive.a
wukong: libkinematics_foot.a
wukong: libkinematics.a
wukong: /opt/raisim/lib/libraisim.so
wukong: /opt/raisim/lib/libraisimPng.so
wukong: /opt/raisim/lib/libraisimZ.so
wukong: /opt/raisim/lib/libraisimODE.so
wukong: /opt/raisim/lib/libraisimMine.so
wukong: /home/rkz/software/root/lib/libCore.so
wukong: /home/rkz/software/root/lib/libImt.so
wukong: /home/rkz/software/root/lib/libRIO.so
wukong: /home/rkz/software/root/lib/libNet.so
wukong: /home/rkz/software/root/lib/libHist.so
wukong: /home/rkz/software/root/lib/libGraf.so
wukong: /home/rkz/software/root/lib/libGraf3d.so
wukong: /home/rkz/software/root/lib/libGpad.so
wukong: /home/rkz/software/root/lib/libROOTDataFrame.so
wukong: /home/rkz/software/root/lib/libTree.so
wukong: /home/rkz/software/root/lib/libTreePlayer.so
wukong: /home/rkz/software/root/lib/libRint.so
wukong: /home/rkz/software/root/lib/libPostscript.so
wukong: /home/rkz/software/root/lib/libMatrix.so
wukong: /home/rkz/software/root/lib/libPhysics.so
wukong: /home/rkz/software/root/lib/libMathCore.so
wukong: /home/rkz/software/root/lib/libThread.so
wukong: /home/rkz/software/root/lib/libMultiProc.so
wukong: /home/rkz/software/root/lib/libROOTVecOps.so
wukong: /usr/lib/x86_64-linux-gnu/libpython3.8.so
wukong: CMakeFiles/wukong.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/rkz/桌面/simulation/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable wukong"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/wukong.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/wukong.dir/build: wukong

.PHONY : CMakeFiles/wukong.dir/build

CMakeFiles/wukong.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/wukong.dir/cmake_clean.cmake
.PHONY : CMakeFiles/wukong.dir/clean

CMakeFiles/wukong.dir/depend:
	cd /home/rkz/桌面/simulation/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/rkz/桌面/simulation /home/rkz/桌面/simulation /home/rkz/桌面/simulation/build /home/rkz/桌面/simulation/build /home/rkz/桌面/simulation/build/CMakeFiles/wukong.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/wukong.dir/depend

