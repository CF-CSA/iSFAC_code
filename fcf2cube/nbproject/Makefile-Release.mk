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
CC=clang
CCC=clang++
CXX=clang++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=CLang-Linux
CND_DLIB_EXT=so
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/Atom.o \
	${OBJECTDIR}/Cubefile.o \
	${OBJECTDIR}/FCFInfo.o \
	${OBJECTDIR}/FCFfile.o \
	${OBJECTDIR}/FCFitem.o \
	${OBJECTDIR}/HKL.o \
	${OBJECTDIR}/HKLops.o \
	${OBJECTDIR}/Int3x3.o \
	${OBJECTDIR}/Mat33.o \
	${OBJECTDIR}/Parser.o \
	${OBJECTDIR}/Reflex.o \
	${OBJECTDIR}/ResFile.o \
	${OBJECTDIR}/Usage.o \
	${OBJECTDIR}/Utils.o \
	${OBJECTDIR}/Vec3.o \
	${OBJECTDIR}/main.o \
	${OBJECTDIR}/sxfft.o


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
LDLIBSOPTIONS=`pkg-config --libs kissfft-float`  

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/fcf2cube

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/fcf2cube: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/fcf2cube ${OBJECTFILES} ${LDLIBSOPTIONS}

${OBJECTDIR}/Atom.o: Atom.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 `pkg-config --cflags kissfft-float` -std=c++14  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Atom.o Atom.cpp

${OBJECTDIR}/Cubefile.o: Cubefile.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 `pkg-config --cflags kissfft-float` -std=c++14  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Cubefile.o Cubefile.cpp

${OBJECTDIR}/FCFInfo.o: FCFInfo.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 `pkg-config --cflags kissfft-float` -std=c++14  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/FCFInfo.o FCFInfo.cpp

${OBJECTDIR}/FCFfile.o: FCFfile.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 `pkg-config --cflags kissfft-float` -std=c++14  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/FCFfile.o FCFfile.cpp

${OBJECTDIR}/FCFitem.o: FCFitem.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 `pkg-config --cflags kissfft-float` -std=c++14  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/FCFitem.o FCFitem.cpp

${OBJECTDIR}/HKL.o: HKL.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 `pkg-config --cflags kissfft-float` -std=c++14  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/HKL.o HKL.cpp

${OBJECTDIR}/HKLops.o: HKLops.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 `pkg-config --cflags kissfft-float` -std=c++14  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/HKLops.o HKLops.cpp

${OBJECTDIR}/Int3x3.o: Int3x3.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 `pkg-config --cflags kissfft-float` -std=c++14  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Int3x3.o Int3x3.cpp

${OBJECTDIR}/Mat33.o: Mat33.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 `pkg-config --cflags kissfft-float` -std=c++14  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Mat33.o Mat33.cpp

${OBJECTDIR}/Parser.o: Parser.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 `pkg-config --cflags kissfft-float` -std=c++14  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Parser.o Parser.cpp

${OBJECTDIR}/Reflex.o: Reflex.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 `pkg-config --cflags kissfft-float` -std=c++14  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Reflex.o Reflex.cpp

${OBJECTDIR}/ResFile.o: ResFile.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 `pkg-config --cflags kissfft-float` -std=c++14  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/ResFile.o ResFile.cpp

${OBJECTDIR}/Usage.o: Usage.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 `pkg-config --cflags kissfft-float` -std=c++14  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Usage.o Usage.cpp

${OBJECTDIR}/Utils.o: Utils.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 `pkg-config --cflags kissfft-float` -std=c++14  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Utils.o Utils.cpp

${OBJECTDIR}/Vec3.o: Vec3.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 `pkg-config --cflags kissfft-float` -std=c++14  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Vec3.o Vec3.cpp

${OBJECTDIR}/main.o: main.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 `pkg-config --cflags kissfft-float` -std=c++14  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/main.o main.cpp

${OBJECTDIR}/sxfft.o: sxfft.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 `pkg-config --cflags kissfft-float` -std=c++14  -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/sxfft.o sxfft.cpp

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
