#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux
CND_DLIB_EXT=so
CND_CONF=Debug
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/Atom.o \
	${OBJECTDIR}/Cube.o \
	${OBJECTDIR}/Mat33.o \
	${OBJECTDIR}/Parser.o \
	${OBJECTDIR}/Usage.o \
	${OBJECTDIR}/Utils.o \
	${OBJECTDIR}/Vec3.o \
	${OBJECTDIR}/cbAtom.o \
	${OBJECTDIR}/kabsch.o \
	${OBJECTDIR}/main.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=`pkg-config --libs gsl`  

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/cubemaps

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/cubemaps: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/cubemaps ${OBJECTFILES} ${LDLIBSOPTIONS}

${OBJECTDIR}/Atom.o: Atom.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g `pkg-config --cflags gsl` -std=c++14  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Atom.o Atom.cpp

${OBJECTDIR}/Cube.o: Cube.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g `pkg-config --cflags gsl` -std=c++14  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Cube.o Cube.cpp

${OBJECTDIR}/Mat33.o: Mat33.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g `pkg-config --cflags gsl` -std=c++14  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Mat33.o Mat33.cpp

${OBJECTDIR}/Parser.o: Parser.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g `pkg-config --cflags gsl` -std=c++14  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Parser.o Parser.cpp

${OBJECTDIR}/Usage.o: Usage.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g `pkg-config --cflags gsl` -std=c++14  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Usage.o Usage.cpp

${OBJECTDIR}/Utils.o: Utils.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g `pkg-config --cflags gsl` -std=c++14  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Utils.o Utils.cpp

${OBJECTDIR}/Vec3.o: Vec3.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g `pkg-config --cflags gsl` -std=c++14  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Vec3.o Vec3.cpp

${OBJECTDIR}/cbAtom.o: cbAtom.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g `pkg-config --cflags gsl` -std=c++14  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/cbAtom.o cbAtom.cpp

${OBJECTDIR}/kabsch.o: kabsch.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g `pkg-config --cflags gsl` -std=c++14  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/kabsch.o kabsch.cpp

${OBJECTDIR}/main.o: main.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g `pkg-config --cflags gsl` -std=c++14  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/main.o main.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
