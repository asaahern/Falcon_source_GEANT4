# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.14

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
CMAKE_COMMAND = /opt/local/bin/cmake

# The command to remove a file.
RM = /opt/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Applications/GEANT4/geant4_10_05_p01/NAUSICAA/FALCON_Xray_source

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Applications/GEANT4/geant4_10_05_p01/NAUSICAA/FALCON_Xray_source/build

# Include any dependencies generated for this target.
include CMakeFiles/FALCON.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/FALCON.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/FALCON.dir/flags.make

CMakeFiles/FALCON.dir/FALCON.cc.o: CMakeFiles/FALCON.dir/flags.make
CMakeFiles/FALCON.dir/FALCON.cc.o: ../FALCON.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Applications/GEANT4/geant4_10_05_p01/NAUSICAA/FALCON_Xray_source/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/FALCON.dir/FALCON.cc.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/FALCON.dir/FALCON.cc.o -c /Applications/GEANT4/geant4_10_05_p01/NAUSICAA/FALCON_Xray_source/FALCON.cc

CMakeFiles/FALCON.dir/FALCON.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/FALCON.dir/FALCON.cc.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Applications/GEANT4/geant4_10_05_p01/NAUSICAA/FALCON_Xray_source/FALCON.cc > CMakeFiles/FALCON.dir/FALCON.cc.i

CMakeFiles/FALCON.dir/FALCON.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/FALCON.dir/FALCON.cc.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Applications/GEANT4/geant4_10_05_p01/NAUSICAA/FALCON_Xray_source/FALCON.cc -o CMakeFiles/FALCON.dir/FALCON.cc.s

CMakeFiles/FALCON.dir/src/ActionInitialization.cc.o: CMakeFiles/FALCON.dir/flags.make
CMakeFiles/FALCON.dir/src/ActionInitialization.cc.o: ../src/ActionInitialization.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Applications/GEANT4/geant4_10_05_p01/NAUSICAA/FALCON_Xray_source/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/FALCON.dir/src/ActionInitialization.cc.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/FALCON.dir/src/ActionInitialization.cc.o -c /Applications/GEANT4/geant4_10_05_p01/NAUSICAA/FALCON_Xray_source/src/ActionInitialization.cc

CMakeFiles/FALCON.dir/src/ActionInitialization.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/FALCON.dir/src/ActionInitialization.cc.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Applications/GEANT4/geant4_10_05_p01/NAUSICAA/FALCON_Xray_source/src/ActionInitialization.cc > CMakeFiles/FALCON.dir/src/ActionInitialization.cc.i

CMakeFiles/FALCON.dir/src/ActionInitialization.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/FALCON.dir/src/ActionInitialization.cc.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Applications/GEANT4/geant4_10_05_p01/NAUSICAA/FALCON_Xray_source/src/ActionInitialization.cc -o CMakeFiles/FALCON.dir/src/ActionInitialization.cc.s

CMakeFiles/FALCON.dir/src/DetectorConstruction.cc.o: CMakeFiles/FALCON.dir/flags.make
CMakeFiles/FALCON.dir/src/DetectorConstruction.cc.o: ../src/DetectorConstruction.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Applications/GEANT4/geant4_10_05_p01/NAUSICAA/FALCON_Xray_source/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/FALCON.dir/src/DetectorConstruction.cc.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/FALCON.dir/src/DetectorConstruction.cc.o -c /Applications/GEANT4/geant4_10_05_p01/NAUSICAA/FALCON_Xray_source/src/DetectorConstruction.cc

CMakeFiles/FALCON.dir/src/DetectorConstruction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/FALCON.dir/src/DetectorConstruction.cc.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Applications/GEANT4/geant4_10_05_p01/NAUSICAA/FALCON_Xray_source/src/DetectorConstruction.cc > CMakeFiles/FALCON.dir/src/DetectorConstruction.cc.i

CMakeFiles/FALCON.dir/src/DetectorConstruction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/FALCON.dir/src/DetectorConstruction.cc.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Applications/GEANT4/geant4_10_05_p01/NAUSICAA/FALCON_Xray_source/src/DetectorConstruction.cc -o CMakeFiles/FALCON.dir/src/DetectorConstruction.cc.s

CMakeFiles/FALCON.dir/src/PhysicsList.cc.o: CMakeFiles/FALCON.dir/flags.make
CMakeFiles/FALCON.dir/src/PhysicsList.cc.o: ../src/PhysicsList.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Applications/GEANT4/geant4_10_05_p01/NAUSICAA/FALCON_Xray_source/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/FALCON.dir/src/PhysicsList.cc.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/FALCON.dir/src/PhysicsList.cc.o -c /Applications/GEANT4/geant4_10_05_p01/NAUSICAA/FALCON_Xray_source/src/PhysicsList.cc

CMakeFiles/FALCON.dir/src/PhysicsList.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/FALCON.dir/src/PhysicsList.cc.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Applications/GEANT4/geant4_10_05_p01/NAUSICAA/FALCON_Xray_source/src/PhysicsList.cc > CMakeFiles/FALCON.dir/src/PhysicsList.cc.i

CMakeFiles/FALCON.dir/src/PhysicsList.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/FALCON.dir/src/PhysicsList.cc.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Applications/GEANT4/geant4_10_05_p01/NAUSICAA/FALCON_Xray_source/src/PhysicsList.cc -o CMakeFiles/FALCON.dir/src/PhysicsList.cc.s

CMakeFiles/FALCON.dir/src/PhysicsListMessenger.cc.o: CMakeFiles/FALCON.dir/flags.make
CMakeFiles/FALCON.dir/src/PhysicsListMessenger.cc.o: ../src/PhysicsListMessenger.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Applications/GEANT4/geant4_10_05_p01/NAUSICAA/FALCON_Xray_source/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/FALCON.dir/src/PhysicsListMessenger.cc.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/FALCON.dir/src/PhysicsListMessenger.cc.o -c /Applications/GEANT4/geant4_10_05_p01/NAUSICAA/FALCON_Xray_source/src/PhysicsListMessenger.cc

CMakeFiles/FALCON.dir/src/PhysicsListMessenger.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/FALCON.dir/src/PhysicsListMessenger.cc.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Applications/GEANT4/geant4_10_05_p01/NAUSICAA/FALCON_Xray_source/src/PhysicsListMessenger.cc > CMakeFiles/FALCON.dir/src/PhysicsListMessenger.cc.i

CMakeFiles/FALCON.dir/src/PhysicsListMessenger.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/FALCON.dir/src/PhysicsListMessenger.cc.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Applications/GEANT4/geant4_10_05_p01/NAUSICAA/FALCON_Xray_source/src/PhysicsListMessenger.cc -o CMakeFiles/FALCON.dir/src/PhysicsListMessenger.cc.s

CMakeFiles/FALCON.dir/src/PrimaryGeneratorAction.cc.o: CMakeFiles/FALCON.dir/flags.make
CMakeFiles/FALCON.dir/src/PrimaryGeneratorAction.cc.o: ../src/PrimaryGeneratorAction.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Applications/GEANT4/geant4_10_05_p01/NAUSICAA/FALCON_Xray_source/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/FALCON.dir/src/PrimaryGeneratorAction.cc.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/FALCON.dir/src/PrimaryGeneratorAction.cc.o -c /Applications/GEANT4/geant4_10_05_p01/NAUSICAA/FALCON_Xray_source/src/PrimaryGeneratorAction.cc

CMakeFiles/FALCON.dir/src/PrimaryGeneratorAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/FALCON.dir/src/PrimaryGeneratorAction.cc.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Applications/GEANT4/geant4_10_05_p01/NAUSICAA/FALCON_Xray_source/src/PrimaryGeneratorAction.cc > CMakeFiles/FALCON.dir/src/PrimaryGeneratorAction.cc.i

CMakeFiles/FALCON.dir/src/PrimaryGeneratorAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/FALCON.dir/src/PrimaryGeneratorAction.cc.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Applications/GEANT4/geant4_10_05_p01/NAUSICAA/FALCON_Xray_source/src/PrimaryGeneratorAction.cc -o CMakeFiles/FALCON.dir/src/PrimaryGeneratorAction.cc.s

CMakeFiles/FALCON.dir/src/SteppingAction.cc.o: CMakeFiles/FALCON.dir/flags.make
CMakeFiles/FALCON.dir/src/SteppingAction.cc.o: ../src/SteppingAction.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Applications/GEANT4/geant4_10_05_p01/NAUSICAA/FALCON_Xray_source/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/FALCON.dir/src/SteppingAction.cc.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/FALCON.dir/src/SteppingAction.cc.o -c /Applications/GEANT4/geant4_10_05_p01/NAUSICAA/FALCON_Xray_source/src/SteppingAction.cc

CMakeFiles/FALCON.dir/src/SteppingAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/FALCON.dir/src/SteppingAction.cc.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Applications/GEANT4/geant4_10_05_p01/NAUSICAA/FALCON_Xray_source/src/SteppingAction.cc > CMakeFiles/FALCON.dir/src/SteppingAction.cc.i

CMakeFiles/FALCON.dir/src/SteppingAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/FALCON.dir/src/SteppingAction.cc.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Applications/GEANT4/geant4_10_05_p01/NAUSICAA/FALCON_Xray_source/src/SteppingAction.cc -o CMakeFiles/FALCON.dir/src/SteppingAction.cc.s

CMakeFiles/FALCON.dir/src/XenonHit.cc.o: CMakeFiles/FALCON.dir/flags.make
CMakeFiles/FALCON.dir/src/XenonHit.cc.o: ../src/XenonHit.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Applications/GEANT4/geant4_10_05_p01/NAUSICAA/FALCON_Xray_source/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/FALCON.dir/src/XenonHit.cc.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/FALCON.dir/src/XenonHit.cc.o -c /Applications/GEANT4/geant4_10_05_p01/NAUSICAA/FALCON_Xray_source/src/XenonHit.cc

CMakeFiles/FALCON.dir/src/XenonHit.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/FALCON.dir/src/XenonHit.cc.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Applications/GEANT4/geant4_10_05_p01/NAUSICAA/FALCON_Xray_source/src/XenonHit.cc > CMakeFiles/FALCON.dir/src/XenonHit.cc.i

CMakeFiles/FALCON.dir/src/XenonHit.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/FALCON.dir/src/XenonHit.cc.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Applications/GEANT4/geant4_10_05_p01/NAUSICAA/FALCON_Xray_source/src/XenonHit.cc -o CMakeFiles/FALCON.dir/src/XenonHit.cc.s

CMakeFiles/FALCON.dir/src/XenonSD.cc.o: CMakeFiles/FALCON.dir/flags.make
CMakeFiles/FALCON.dir/src/XenonSD.cc.o: ../src/XenonSD.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Applications/GEANT4/geant4_10_05_p01/NAUSICAA/FALCON_Xray_source/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object CMakeFiles/FALCON.dir/src/XenonSD.cc.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/FALCON.dir/src/XenonSD.cc.o -c /Applications/GEANT4/geant4_10_05_p01/NAUSICAA/FALCON_Xray_source/src/XenonSD.cc

CMakeFiles/FALCON.dir/src/XenonSD.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/FALCON.dir/src/XenonSD.cc.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Applications/GEANT4/geant4_10_05_p01/NAUSICAA/FALCON_Xray_source/src/XenonSD.cc > CMakeFiles/FALCON.dir/src/XenonSD.cc.i

CMakeFiles/FALCON.dir/src/XenonSD.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/FALCON.dir/src/XenonSD.cc.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Applications/GEANT4/geant4_10_05_p01/NAUSICAA/FALCON_Xray_source/src/XenonSD.cc -o CMakeFiles/FALCON.dir/src/XenonSD.cc.s

# Object files for target FALCON
FALCON_OBJECTS = \
"CMakeFiles/FALCON.dir/FALCON.cc.o" \
"CMakeFiles/FALCON.dir/src/ActionInitialization.cc.o" \
"CMakeFiles/FALCON.dir/src/DetectorConstruction.cc.o" \
"CMakeFiles/FALCON.dir/src/PhysicsList.cc.o" \
"CMakeFiles/FALCON.dir/src/PhysicsListMessenger.cc.o" \
"CMakeFiles/FALCON.dir/src/PrimaryGeneratorAction.cc.o" \
"CMakeFiles/FALCON.dir/src/SteppingAction.cc.o" \
"CMakeFiles/FALCON.dir/src/XenonHit.cc.o" \
"CMakeFiles/FALCON.dir/src/XenonSD.cc.o"

# External object files for target FALCON
FALCON_EXTERNAL_OBJECTS =

FALCON: CMakeFiles/FALCON.dir/FALCON.cc.o
FALCON: CMakeFiles/FALCON.dir/src/ActionInitialization.cc.o
FALCON: CMakeFiles/FALCON.dir/src/DetectorConstruction.cc.o
FALCON: CMakeFiles/FALCON.dir/src/PhysicsList.cc.o
FALCON: CMakeFiles/FALCON.dir/src/PhysicsListMessenger.cc.o
FALCON: CMakeFiles/FALCON.dir/src/PrimaryGeneratorAction.cc.o
FALCON: CMakeFiles/FALCON.dir/src/SteppingAction.cc.o
FALCON: CMakeFiles/FALCON.dir/src/XenonHit.cc.o
FALCON: CMakeFiles/FALCON.dir/src/XenonSD.cc.o
FALCON: CMakeFiles/FALCON.dir/build.make
FALCON: /Applications/GEANT4/geant4_10_05_p01-install/lib/libG4Tree.dylib
FALCON: /Applications/GEANT4/geant4_10_05_p01-install/lib/libG4GMocren.dylib
FALCON: /Applications/GEANT4/geant4_10_05_p01-install/lib/libG4visHepRep.dylib
FALCON: /Applications/GEANT4/geant4_10_05_p01-install/lib/libG4RayTracer.dylib
FALCON: /Applications/GEANT4/geant4_10_05_p01-install/lib/libG4VRML.dylib
FALCON: /Applications/GEANT4/geant4_10_05_p01-install/lib/libG4OpenGL.dylib
FALCON: /Applications/GEANT4/geant4_10_05_p01-install/lib/libG4gl2ps.dylib
FALCON: /Applications/GEANT4/geant4_10_05_p01-install/lib/libG4interfaces.dylib
FALCON: /Applications/GEANT4/geant4_10_05_p01-install/lib/libG4persistency.dylib
FALCON: /Applications/GEANT4/geant4_10_05_p01-install/lib/libG4error_propagation.dylib
FALCON: /Applications/GEANT4/geant4_10_05_p01-install/lib/libG4readout.dylib
FALCON: /Applications/GEANT4/geant4_10_05_p01-install/lib/libG4physicslists.dylib
FALCON: /Applications/GEANT4/geant4_10_05_p01-install/lib/libG4parmodels.dylib
FALCON: /Applications/GEANT4/geant4_10_05_p01-install/lib/libG4FR.dylib
FALCON: /Applications/GEANT4/geant4_10_05_p01-install/lib/libG4vis_management.dylib
FALCON: /Applications/GEANT4/geant4_10_05_p01-install/lib/libG4modeling.dylib
FALCON: /opt/local/lib/libSM.dylib
FALCON: /opt/local/lib/libICE.dylib
FALCON: /opt/local/lib/libX11.dylib
FALCON: /opt/local/lib/libXext.dylib
FALCON: /opt/X11/lib/libGL.dylib
FALCON: /opt/X11/lib/libGLU.dylib
FALCON: /opt/X11/lib/libXmu.dylib
FALCON: //anaconda3/lib/libQt5OpenGL.5.9.7.dylib
FALCON: //anaconda3/lib/libQt5PrintSupport.5.9.7.dylib
FALCON: //anaconda3/lib/libQt5Widgets.5.9.7.dylib
FALCON: //anaconda3/lib/libQt5Gui.5.9.7.dylib
FALCON: //anaconda3/lib/libQt5Core.5.9.7.dylib
FALCON: /opt/local/lib/libxerces-c.dylib
FALCON: /Applications/GEANT4/geant4_10_05_p01-install/lib/libG4run.dylib
FALCON: /Applications/GEANT4/geant4_10_05_p01-install/lib/libG4event.dylib
FALCON: /Applications/GEANT4/geant4_10_05_p01-install/lib/libG4tracking.dylib
FALCON: /Applications/GEANT4/geant4_10_05_p01-install/lib/libG4processes.dylib
FALCON: /Applications/GEANT4/geant4_10_05_p01-install/lib/libG4analysis.dylib
FALCON: /opt/local/lib/libfreetype.dylib
FALCON: /Applications/GEANT4/geant4_10_05_p01-install/lib/libG4zlib.dylib
FALCON: /opt/local/lib/libexpat.dylib
FALCON: /Applications/GEANT4/geant4_10_05_p01-install/lib/libG4digits_hits.dylib
FALCON: /Applications/GEANT4/geant4_10_05_p01-install/lib/libG4track.dylib
FALCON: /Applications/GEANT4/geant4_10_05_p01-install/lib/libG4particles.dylib
FALCON: /Applications/GEANT4/geant4_10_05_p01-install/lib/libG4geometry.dylib
FALCON: /Applications/GEANT4/geant4_10_05_p01-install/lib/libG4materials.dylib
FALCON: /Applications/GEANT4/geant4_10_05_p01-install/lib/libG4graphics_reps.dylib
FALCON: /Applications/GEANT4/geant4_10_05_p01-install/lib/libG4intercoms.dylib
FALCON: /Applications/GEANT4/geant4_10_05_p01-install/lib/libG4global.dylib
FALCON: /opt/local/lib/libCLHEP-2.4.1.0.dylib
FALCON: CMakeFiles/FALCON.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Applications/GEANT4/geant4_10_05_p01/NAUSICAA/FALCON_Xray_source/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Linking CXX executable FALCON"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/FALCON.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/FALCON.dir/build: FALCON

.PHONY : CMakeFiles/FALCON.dir/build

CMakeFiles/FALCON.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/FALCON.dir/cmake_clean.cmake
.PHONY : CMakeFiles/FALCON.dir/clean

CMakeFiles/FALCON.dir/depend:
	cd /Applications/GEANT4/geant4_10_05_p01/NAUSICAA/FALCON_Xray_source/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Applications/GEANT4/geant4_10_05_p01/NAUSICAA/FALCON_Xray_source /Applications/GEANT4/geant4_10_05_p01/NAUSICAA/FALCON_Xray_source /Applications/GEANT4/geant4_10_05_p01/NAUSICAA/FALCON_Xray_source/build /Applications/GEANT4/geant4_10_05_p01/NAUSICAA/FALCON_Xray_source/build /Applications/GEANT4/geant4_10_05_p01/NAUSICAA/FALCON_Xray_source/build/CMakeFiles/FALCON.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/FALCON.dir/depend

