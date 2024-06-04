cmake_minimum_required(VERSION 3.14)
# because we need target specific compile and link flags,
# which were introduced in 3.14

MESSAGE(STATUS "====================================================")
MESSAGE(STATUS "============ Configuring GTECTON ===================")
MESSAGE(STATUS "====================================================")

project(GTECTON VERSION 1.0
                DESCRIPTION "GTECTON"
                LANGUAGES C Fortran)

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/sparseModFiles)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/nonSparseModFiles)

set(sparseModDir ${CMAKE_BINARY_DIR}/sparseModFiles)
set(nonSparseModDir ${CMAKE_BINARY_DIR}/nonSparseModFiles)

SET(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS f;F;f77;F77;f90;F90;f95;F95)

# modules are in the path 'cmake'

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmakeModules")

MESSAGE(STATUS "Using CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}")

##########################################################################
# check for dependencies
##########################################################################

if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
  if(CMAKE_C_COMPILER_VERSION VERSION_LESS 7.0.0)
    message(FATAL_ERROR "Needs at least gcc version 7")
  endif()
else()
  message(FATAL ERROR "GTECTON needs the GNU compiler set.")
endif()

# showme needs X11.
find_package(X11 REQUIRED)
message(STATUS "X11_FOUND = ${X11_FOUND}")
message(STATUS "X11_INCLUDE_DIR = ${X11_INCLUDE_DIR}")
message(STATUS "X11_LIBRARIES = ${X11_LIBRARIES}")

# partition needs METIS.
find_package(METIS REQUIRED)
message(STATUS "METIS_FOUND = ${METIS_FOUND}")
message(STATUS "METIS_INCLUDE_DIR = ${METIS_INCLUDE_DIR}")
message(STATUS "METIS_LIBRARIES = ${METIS_LIBRARIES}")


find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

# Check whether we have the correct PETSc
IF(DEFINED ENV{PETSC_DIR})
    set(PETSC_DIR $ENV{PETSC_DIR})
    MESSAGE(STATUS "Found environment variable PETSC_DIR: ${PETSC_DIR}")
ELSE()
    MESSAGE(FATAL_ERROR "environment variable PETSC_DIR no set. Must be set to root or PETSc installation")
ENDIF()

IF(DEFINED ENV{PETSC_ARCH})
    set(PETSC_ARCH $ENV{PETSC_ARCH})
    MESSAGE(STATUS "Found environment variable PETSC_DIR: ${PETSC_ARCH}")
ELSE()
    MESSAGE(FATAL_ERROR "environment variable PETSC_ARCH no set. Must be set to architecture name within root of PETSc installation, such as arch-linux2-c-opt")
ENDIF()


set(PETSC_INCLUDES " ")
set(PETSC_LIBRARIES " ")
set(PETSC_EXECUTABLE_RUNS " ")

find_package(PETSc REQUIRED)

if(${PETSC_FOUND} STREQUAL "TRUE")
    MESSAGE(STATUS "Found PETSc.")
else()
    MESSAGE(FATAL_ERROR "Could not find PETSc. Is environment variable PETSC_DIR set correctly?")
endif()

#if(${PETSC_VERSION} STREQUAL "3.17.2")
#    MESSAGE(STATUS "Confirmed correct PETSc version 3.10.4")
#else()
#    MESSAGE(FATAL_ERROR "Need PETSc version 3.10.4 Found version" ${PETSC_VERSION})
#endif()


########## BLAS ###########

find_package(BLAS)

if(${BLAS_FOUND} STREQUAL "TRUE")
	message(STATUS "Found BLAS libraries ${BLAS_LIBRARIES}")
    message(STATUS "Found BLAS linker flags ${BLAS_LINKER_FLAGS}")
else()
	message(FATAL_ERROR "Did not find BLAS")
endif()

########## LAPACK ###########

find_package(LAPACK)

if(${LAPACK_FOUND} STREQUAL "TRUE")
    message(STATUS "Found LAPACK libraries ${LAPACK_LIBRARIES}")
    message(STATUS "Found LAPACK linker flags ${LAPACK_LINKER_FLAGS}")
else()
    message(FATAL_ERROR "Did not find LAPACK")
endif()

########## MPI ###########

find_package(MPI)

if(${MPI_FOUND} STREQUAL "TRUE")
    message(STATUS "Found MPI libraries ${MPI_Fortran_LIBRARIES}")
    message(STATUS "Found MPI compiler ${MPI_Fortran_COMPILER}")

else()
    message(FATAL_ERROR "Did not find MPI")
endif()


# Commented block below prints every known cmake variable.
# There are hundreds, hence commented out.
# Its usefulness for debugging cannot be understated, so please leave it in.

#message(STATUS "************************************************************")
#get_cmake_property(_variableNames VARIABLES)
#list (SORT _variableNames)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()
#message(STATUS "************************************************************")

##########################################################################
# Set paths to all the individual components
##########################################################################

# 01 meshing
set(dirTriangle        src/01_meshing/2D/triangle)
set(dirTriangleHelp    src/01_meshing/2D/triangle-helpers)
set(dirGmshHelp        src/01_meshing/3D/gmsh-helpers)

# 02 partioning
set(dirPart            src/02_partitioning/partition)

# 03 preprocessing
set(dirPre             src/03_preprocessing)
set(dirPreElmside      src/03_preprocessing/elmside)

# 04 simulating
set(dirSimGeneral      src/04_simulating/general)
set(dirSimSpecificPLN  src/04_simulating/2D/pln)
set(dirSimSpecificOPN  src/04_simulating/2D/opn)
set(dirSimSpecificAXI  src/04_simulating/2D/axi)
set(dirSimSpecificF3D  src/04_simulating/3D/f3d)

# 05 merging
set(dirMerge           src/05_merging)

# 06 postprocessing
set(dirPostGeneral2D   src/06_postprocessing/2D/general)
set(dirPostSpecificPLN src/06_postprocessing/2D/pln)
set(dirPostSpecificOPN src/06_postprocessing/2D/opn)
set(dirPostSpecificAXI src/06_postprocessing/2D/axi)

set(dirPostGeneral3D   src/06_postprocessing/3D)

# lib misc
set(dirLibMisc         src/lib_misc)

# tools
set(dirTools           src/tools)



##########################################################################
# set compiler options
##########################################################################

set(CMAKE_Fortran_COMPILER "${MPI_Fortran_COMPILER}")

# todo, mets include should only be at the relevant c files.
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} \
-DFORTRAN_UNDERSCORE \
-I${METIS_INCLUDE_DIR}")

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} \
-cpp \
-fno-align-commons \
-DPETSC_USE_REAL_DOUBLE \
-std=f2018 \
-ffree-form \
-L${PETSC_DIR}/${PETSC_ARCH}/lib \
-fopenmp \
-I../src \
-I${PETSC_DIR}/include \
-I${PETSC_DIR}/${PETSC_ARCH}/include \
-I../${dirSimGeneral} \
-ffree-line-length-none \
-DLongData \
-DGIA \
-g \
-fcheck=bounds \
-DDIFORC_SET \
-fimplicit-none")

#-DEXPERIMENTALTRANSIENTDISCREEP \
# -DSHEARZONE_TRACTIONS \
# -DEXPERIMENTAL_CYCLIC \
# -DZISFAULTNORMAL \
#-DEXPERIMENTAL_ReducedSlipAlloc \
#-fallow-argument-mismatch")
#-Wall \
#-Wunused \

#-DEXPERIMENTAL_FRICTION \

INCLUDE_DIRECTORIES(BEFORE "${dirSimGeneral}")


##########################################################################
# group the source files, so that they can be used for multiple
# executables
##########################################################################

set(picknpsSource ${dirPre}/picknps.c)


list(APPEND libMiscSource
"${dirLibMisc}/ipeven.f"
"${dirLibMisc}/openf.F"
"${dirLibMisc}/ifopen.F"
"${dirLibMisc}/nextlu.F"
"${dirLibMisc}/prerr.F"
"${dirLibMisc}/lnblk.F"
"${dirLibMisc}/lueof.F"
"${dirLibMisc}/exitp.F"
"${dirLibMisc}/exitpc.c"
"${dirLibMisc}/iflu.F"
"${dirLibMisc}/inmax.F"
"${dirLibMisc}/closef.F"
"${dirLibMisc}/dirname.F"
"${dirLibMisc}/tmpnm.F"
"${dirLibMisc}/flushf.F"
"${dirLibMisc}/flushc.c"
"${dirLibMisc}/float.F"
"${dirLibMisc}/nint.F"
"${dirLibMisc}/fnblnk.F"
"${dirLibMisc}/fbopen.F"
"${dirLibMisc}/how2use.F"
"${dirLibMisc}/odd.F"
"${dirLibMisc}/adjustl.F"
"${dirLibMisc}/tmpdir.F"
"${dirLibMisc}/triostat.F"
"${dirLibMisc}/rm.c"
"${dirLibMisc}/gethome.c"
"${dirLibMisc}/soelim.c"
"${dirLibMisc}/csystem.c"
"${dirLibMisc}/gethome.c"
"${dirLibMisc}/numeric.F"
"${dirLibMisc}/ifbopen.F"
"${dirLibMisc}/atan.F"
"${dirLibMisc}/signalc.c"
"${dirLibMisc}/wGMTfile.F"
"${dirLibMisc}/indx.F"
"${dirLibMisc}/chreal.F"
"${dirLibMisc}/wround.F"
"${dirLibMisc}/AsciF.F"
"${dirLibMisc}/cstat.c"
"${dirLibMisc}/input.F"
"${dirLibMisc}/isbinary.c"
"${dirLibMisc}/varifmt.F"
"${dirLibMisc}/cartopolar.F"
"${dirLibMisc}/getwdir.F"
"${dirLibMisc}/cwd.c"
"${dirLibMisc}/smallSleep.c")


# put all the modules in their own target, because they depend
# on each other and need to be compiled in the correct sequence
list(APPEND moduleGITSource
"${dirSimGeneral}/gitModule.F")
list(APPEND modulePostProcessData
"${dirSimGeneral}/postProcessData.F")
list(APPEND moduleConstantsSource
"${dirSimGeneral}/constants.F")
list(APPEND moduleErrorsSource
"${dirSimGeneral}/errorModule.F")
list(APPEND moduleTopologySource
"${dirSimGeneral}/modelTopology.F")
list(APPEND moduleAlgebraSource
"${dirSimGeneral}/algebra.F")
list(APPEND moduleCTXSource
"${dirSimGeneral}/modelctx.F")
list(APPEND moduleDocumentationSource
"${dirSimGeneral}/documentation.F")
list(APPEND moduleConversionsSource
"${dirSimGeneral}/conversions.F")
list(APPEND moduleFileSource
"${dirSimGeneral}/filemodule.F")
list(APPEND moduleAOSource
"${dirSimGeneral}/AOmodule.F")
list(APPEND moduleDebugSource
"${dirSimGeneral}/debugmodule.F"
"${dirLibMisc}/smallSleep.c")
list(APPEND moduleMeshdataSource
"${dirSimGeneral}/meshdatamodule.F")
list(APPEND moduleFormatsSource
"${dirSimGeneral}/formatsmodule.F")
list(APPEND moduleIOSource
"${dirSimGeneral}/iomodule.F")
list(APPEND moduleLabelSource
"${dirSimGeneral}/labelmodule.F")
list(APPEND moduleMaterialsSource
"${dirSimGeneral}/materialsmodule.F")
list(APPEND moduleModelDefinitionSource
"${dirSimGeneral}/modeldefinition.F")
list(APPEND moduleModelDataSource
"${dirSimGeneral}/modeldatamodule.F")
list(APPEND modulePETScSource
"${dirSimGeneral}/PETScmodule.F")
list(APPEND moduleTimeStepSource
"${dirSimGeneral}/timestepmodule.F")
list(APPEND moduleExtremaSource
"${dirSimGeneral}/extremamodule.F")
list(APPEND moduleFrictionSource
"${dirSimGeneral}/frictionmodule.F")
list(APPEND moduleTestSource
"${dirSimGeneral}/testmodule.F")


list(APPEND libGtectonGeneralSource
"${dirSimGeneral}/tecin.F"
"${dirSimGeneral}/rdthrm.F"
"${dirSimGeneral}/tempdat.F"
"${dirSimGeneral}/visco.F"
"${dirSimGeneral}/elastc.F"
"${dirSimGeneral}/setprp.F"
"${dirSimGeneral}/timsol.F"
"${dirSimGeneral}/tmarch.F"
"${dirSimGeneral}/getfil.F"
"${dirSimGeneral}/steady.F"
"${dirSimGeneral}/addan.F"
"${dirSimGeneral}/tecdat.F"
"${dirSimGeneral}/temp0.F"
"${dirSimGeneral}/setio.F")


list(APPEND libGtectonMainSource
"${dirSimGeneral}/gtecton.F")

list(APPEND f3dlibSource
"${dirSimSpecificF3D}/f3dlib.F")

list(APPEND opnlibSource
"${dirSimSpecificOPN}/opnlib.F")

list(APPEND plnlibSource
"${dirSimSpecificPLN}/plnlib.F")

list(APPEND axilibSource
"${dirSimSpecificAXI}/axilib.F")

list(APPEND plt2DSource
"${dirPostGeneral2D}/feplt2D.F"
"${dirPostGeneral2D}/setpltio.F"
"${dirPostGeneral2D}/pltout.F"
"${dirPostGeneral2D}/pltdef2D.F"
"${dirPostGeneral2D}/parseCommandLineArguments2D.F"
"${dirPostGeneral2D}/igetcmd2D.F"
"${dirPostGeneral2D}/xeclib2D.F"
"${dirPostGeneral2D}/pltinp.F"
"${dirPostGeneral2D}/setcmd2D.F"
"${dirPostGeneral2D}/bufcmd2D.F"
"${dirPostGeneral2D}/show2D.F"
"${dirPostGeneral2D}/bc.F"
"${dirPostGeneral2D}/arrow.F"
"${dirPostGeneral2D}/stress.F"
"${dirPostGeneral2D}/postel.F"
"${dirPostGeneral2D}/strain.F"
"${dirPostGeneral2D}/pltio.F"
"${dirPostGeneral2D}/shvar2D.F"
"${dirPostGeneral2D}/elmsys.F"
"${dirPostGeneral2D}/setmon2D.F"
"${dirPostGeneral2D}/contur.F"
"${dirPostGeneral2D}/timinc2D.F"
"${dirPostGeneral2D}/dt.F"
"${dirPostGeneral2D}/mats.F"
"${dirPostGeneral2D}/varfit2D.F"
"${dirPostGeneral2D}/makeij.F"
"${dirPostGeneral2D}/flux.F"
"${dirPostGeneral2D}/trian.F"
"${dirPostGeneral2D}/displ.F"
"${dirPostGeneral2D}/partitionBoundaries.F")

list(APPEND plt3DSource
"${dirPostGeneral3D}/base64encode.F"
"${dirPostGeneral3D}/vtkmodule.F"
"${dirPostGeneral3D}/feplt3D.F"
"${dirPostGeneral3D}/parseCommandLineArguments3D.F"
"${dirPostGeneral3D}/setcmd3D.F"
"${dirPostGeneral3D}/xeclib3D.F"
"${dirPostGeneral3D}/bufcmd3D.F"
"${dirPostGeneral3D}/show3D.F"
"${dirPostGeneral3D}/setmon3D.F"
"${dirPostGeneral3D}/timinc3D.F"
"${dirPostGeneral3D}/varfit3D.F"
"${dirPostGeneral3D}/feio.F"
"${dirPostGeneral3D}/pltcon.F"
"${dirPostGeneral3D}/pltdef3D.F"
"${dirPostGeneral3D}/igetcmd3D.F"
"${dirPostGeneral3D}/wait.F"
"${dirPostGeneral3D}/iprinfo.F"
"${dirPostGeneral3D}/echo.F"
"${dirPostGeneral3D}/shvar3D.F"
"${dirPostGeneral3D}/ofopen.F"
"${dirPostGeneral3D}/contour.F"
"${dirPostGeneral3D}/ofclose.F"
"${dirPostGeneral3D}/wreffs.F"
"${dirPostGeneral3D}/wrtens.F"
"${dirPostGeneral3D}/wrvect.F"
"${dirPostGeneral3D}/sample.F"
"${dirPostGeneral3D}/calcs.F"
"${dirPostGeneral3D}/ebc.F"
"${dirPostGeneral3D}/wreslv.F"
"${dirPostGeneral3D}/infrot.F"
"${dirPostGeneral3D}/wcalct.F"
"${dirPostGeneral3D}/wrdx.c"
"${dirPostGeneral3D}/wrealt.F"
"${dirPostGeneral3D}/updatefault.F"
"${dirPostGeneral3D}/grav3D.F"
"${dirPostGeneral3D}/error.f"
"${dirPostGeneral3D}/readGmsh.F")

list(APPEND plnpltSource
"${dirPostSpecificPLN}/plncon.F"
"${dirPostSpecificPLN}/plneff.F"
"${dirPostSpecificPLN}/plngrav.F"
"${dirPostSpecificPLN}/plnplot.F"
"${dirSimSpecificPLN}/plnlib.F")

list(APPEND axipltSource
"${dirPostSpecificAXI}/axicon.F"
"${dirPostSpecificAXI}/axieff.F"
"${dirPostSpecificAXI}/axigrav.F"
"${dirPostSpecificAXI}/axiplot.F"
"${dirSimSpecificAXI}/axilib.F")

list(APPEND opnpltSource
"${dirPostSpecificOPN}/opncon.F"
"${dirPostSpecificOPN}/opneff.F"
"${dirPostSpecificOPN}/opngrav.F"
"${dirPostSpecificOPN}/opnplot.F"
"${dirSimSpecificOPN}/opnlib.F")

list(APPEND triangleSource
"${dirTriangle}/triangle.c")

list(APPEND showmeSource
"${dirTriangle}/showme.c")

list(APPEND elmsideSource
"${dirPreElmside}/globals.f"
"${dirPreElmside}/stringHandling.f"
"${dirPreElmside}/numInArray.f"
"${dirPreElmside}/elmside.F"
"${dirPreElmside}/algorithm1_circleCenters.f"
"${dirPreElmside}/algorithm2_classic.f"
"${dirPreElmside}/algorithm3_classic.f"
"${dirPreElmside}/readCommandLineArgs.f"
"${dirPreElmside}/readInput.f"
"${dirPreElmside}/writeOutput.f")

list(APPEND setAreaSource
"${dirTriangleHelp}/setarea.f")

list(APPEND polysortSource
"${dirPre}/polysort.f")

list(APPEND mergefilesSource
"${dirMerge}/mergermodule.F"
"${dirMerge}/mergers.F")

list(APPEND gmsh2tectonSource
"${dirGmshHelp}/GmshToTecton.F")

list(APPEND rmdoublesSource
"${dirGmshHelp}/rmdoubles.F")

list(APPEND getnormalsSource
"${dirGmshHelp}/getnormals.f90")

list(APPEND setpremSource
"${dirGmshHelp}/setprem.f90")

list(APPEND partition2paraviewSource
"${dirPart}/partition2paraview.f")

list(APPEND partitionSource
"${dirPart}/partition.c")

list(APPEND tri2feSource
"${dirTriangleHelp}/tri2fe.f")

list(APPEND domainUnitySource
"${dirTriangleHelp}/domainUnity.F")

list(APPEND PartitionedGmshToTectonSource
"${dirGmshHelp}/PartitionedGmshToTecton.F")

list(APPEND partitionInterfaceCountSource
"${dirGmshHelp}/partitionInterfaceCount.f")

list(APPEND relvelSource
"${dirLibMisc}/relvel.f")

list(APPEND regridSource
"${dirSimGeneral}/regrid.F")



##########################################################################
# Group the module targets for compactness, as they are on multiple
# occasions required all together
##########################################################################

list(APPEND  allSparseModules
     SmoduleGIT
     SmoduleConstants
     SmoduleErrors
     SmodulePostProcessData
     SmoduleTopology
     SmoduleAlgebra
     SmoduleDocumentation
     SmoduleConversions
     SmoduleCTX
     SmoduleFile
     SmoduleAO
     SmoduleDebug
     SmoduleMeshdata
     SmoduleFormats
     SmoduleIO
     SmoduleLabel
     SmoduleMaterials
     SmoduleModelDefinition
     SmoduleModelData
     SmodulePETSc
     SmoduleTimeStep
     SmoduleExtrema
     SmoduleFriction
     SmoduleTest)

list(APPEND  allNonSparseModules
     NSmoduleGIT
     NSmoduleConstants
     NSmoduleErrors
     NSmodulePostProcessData
     NSmoduleTopology
     NSmoduleAlgebra
     NSmoduleDocumentation
     NSmoduleConversions
     NSmoduleCTX
     NSmoduleFile
     NSmoduleAO
     NSmoduleDebug
     NSmoduleMeshdata
     NSmoduleFormats
     NSmoduleIO
     NSmoduleLabel
     NSmoduleMaterials
     NSmoduleModelDefinition
     NSmoduleModelData
     NSmodulePETSc
     NSmoduleTimeStep
     NSmoduleExtrema
     NSmoduleFriction
     NSmoduleTest)


#MESSAGE(STATUS "Extracting GIT hash")
# put GIT hash invariable so that it can be used in setio and appear in GTecton.rc
# before anything is compiled.
#execute_process(
#        COMMAND ./getGitData.sh
#        WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
#        OUTPUT_VARIABLE GIT_HASH
#        )



##########################################################################
# define compile targets
##########################################################################

# Note that Constants, Topology, Algebra, Documentation and Conversions
# do not have SPARSE tags, and in theory could be build only once for both the
# SPARSE and the non-SPARSE archive.

# However, in case such a split will arise in the future, for robustness
# they are built separately in both forms.


############# NON SPARSE version of the modules, for plnplt etc ###################

add_library(NSmoduleGIT STATIC ${moduleGITSource})
# no dependencies

add_library(NSmoduleConstants STATIC ${moduleConstantsSource})
# no dependencies

add_library(NSmoduleErrors STATIC ${moduleErrorsSource})
# no dependencies

add_library(NSmodulePostProcessData STATIC ${modulePostProcessData})
add_dependencies(NSmodulePostProcessData NSmoduleModelDefinition)
# no dependencies

add_library(NSmoduleTopology STATIC ${moduleTopologySource})
# no dependencies

add_library(NSmoduleAlgebra STATIC ${moduleAlgebraSource})
add_dependencies(NSmoduleAlgebra NSmoduleConstants NSmoduleIO)

add_library(NSmoduleDocumentation STATIC ${moduleDocumentationSource})
add_dependencies(NSmoduleDocumentation NSmoduleConstants)

add_library(NSmoduleConversions STATIC ${moduleConversionsSource})
# no dependencies

add_library(NSmoduleCTX STATIC ${moduleCTXSource})
# no dependencies

add_library(NSmoduleFile STATIC ${moduleFileSource})
add_dependencies(NSmoduleFile NSmoduleCTX)

add_library(NSmoduleAO STATIC ${moduleAOSource})
add_dependencies(NSmoduleAO NSmoduleFile NSmoduleDebug NSmoduleTimeStep NSmoduleModelDefinition NSmoduleMeshdata NSmoduleModelData NSmoduleCTX NSmoduleTopology NSmodulePETSc NSmoduleAlgebra NSmoduleConstants)

add_library(NSmoduleDebug STATIC ${moduleDebugSource})
add_dependencies(NSmoduleDebug NSmoduleAlgebra NSmoduleCTX NSmoduleFile NSmodulePETSc)

add_library(NSmoduleMeshdata STATIC ${moduleMeshdataSource})
add_dependencies(NSmoduleMeshdata NSmoduleCTX NSmoduleFile NSmoduleTopology NSmoduleConstants)

add_library(NSmoduleFormats STATIC ${moduleFormatsSource})
# no dependencies

add_library(NSmoduleIO STATIC ${moduleIOSource})
# no dependencies

add_library(NSmoduleLabel STATIC ${moduleLabelSource})
# no dependencies

add_library(NSmoduleMaterials STATIC ${moduleMaterialsSource})
add_dependencies(NSmoduleMaterials NSmoduleMeshdata NSmoduleDebug NSmoduleCTX)

add_library(NSmoduleModelDefinition STATIC ${moduleModelDefinitionSource})
# no dependencies

add_library(NSmoduleModelData STATIC ${moduleModelDataSource})
add_dependencies(NSmoduleModelData NSmoduleMeshdata NSmoduleModelDefinition NSmoduleDebug NSmoduleTopology)
add_dependencies(NSmoduleModelData triangle)

add_library(NSmodulePETSc STATIC ${modulePETScSource})
add_dependencies(NSmodulePETSc NSmoduleIO NSmoduleMeshdata NSmoduleFile NSmoduleCTX NSmoduleModelDefinition)

add_library(NSmoduleTimeStep STATIC ${moduleTimeStepSource})
add_dependencies(NSmoduleTimeStep NSmoduleDebug NSmoduleCTX NSmoduleConstants NSmoduleModelDefinition)
target_link_libraries(NSmoduleTimeStep libMisc)

add_library(NSmoduleExtrema STATIC ${moduleExtremaSource})
# no dependencies

add_library(NSmoduleFriction STATIC ${moduleFrictionSource})
add_dependencies(NSmoduleFriction NSmoduleModelData)
target_link_libraries(NSmoduleFriction ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

add_library(NSmoduleTest STATIC ${moduleTestSource})
add_dependencies(NSmoduleTest NSmoduleFriction)

### make sure the .mod files all end up in the proper directory ${nonSparseModDir},
### in stead of in the main compile directory.
set_target_properties(${allNonSparseModules} PROPERTIES Fortran_MODULE_DIRECTORY ${nonSparseModDir})



############# SPARSE version of the modules, for pln etc ###################

# Not al modules have an explicit sparse equivalent, (constants, labels, etc)
# So the sparse modules also depend on the non-sparse equivalents.


add_library(SmoduleGIT STATIC ${moduleGITSource})
target_compile_options(SmoduleGIT PUBLIC "-DSPARSE")
# no dependencies

add_library(SmoduleConstants STATIC ${moduleConstantsSource})
target_compile_options(SmoduleConstants PUBLIC "-DSPARSE")
# no dependencies

add_library(SmoduleErrors STATIC ${moduleErrorsSource})
target_compile_options(SmoduleErrors PUBLIC "-DSPARSE")
# no dependencies

add_library(SmodulePostProcessData STATIC ${modulePostProcessData})
add_dependencies(SmodulePostProcessData SmoduleModelDefinition)
target_compile_options(SmodulePostProcessData PUBLIC "-DSPARSE")
# no dependencies

add_library(SmoduleTopology STATIC ${moduleTopologySource})
target_compile_options(SmoduleTopology PUBLIC "-DSPARSE")
# no dependencies

add_library(SmoduleAlgebra STATIC ${moduleAlgebraSource})
target_compile_options(SmoduleAlgebra PUBLIC "-DSPARSE")
add_dependencies(SmoduleAlgebra SmoduleConstants SmoduleIO)

add_library(SmoduleDocumentation STATIC ${moduleDocumentationSource})
add_dependencies(SmoduleDocumentation SmoduleConstants)
target_compile_options(SmoduleDocumentation PUBLIC "-DSPARSE")

add_library(SmoduleConversions STATIC ${moduleConversionsSource})
target_compile_options(SmoduleConversions PUBLIC "-DSPARSE")
# no dependencies

add_library(SmoduleCTX STATIC ${moduleCTXSource})
target_compile_options(SmoduleCTX PUBLIC "-DSPARSE")
# no dependencies

add_library(SmoduleFile STATIC ${moduleFileSource})
target_compile_options(SmoduleFile PUBLIC "-DSPARSE")
add_dependencies(SmoduleFile SmoduleCTX)

add_library(SmoduleAO STATIC ${moduleAOSource})
target_compile_options(SmoduleAO PUBLIC "-DSPARSE")
add_dependencies(SmoduleAO SmoduleFile SmoduleDebug SmoduleTimeStep SmoduleModelDefinition SmoduleMeshdata SmoduleModelData SmoduleCTX SmoduleTopology SmodulePETSc SmoduleAlgebra SmoduleConstants)

add_library(SmoduleDebug STATIC ${moduleDebugSource})
target_compile_options(SmoduleDebug PUBLIC "-DSPARSE")
add_dependencies(SmoduleDebug SmoduleAlgebra SmoduleCTX SmoduleFile SmodulePETSc)

add_library(SmoduleMeshdata STATIC ${moduleMeshdataSource})
target_compile_options(SmoduleMeshdata PUBLIC "-DSPARSE")
add_dependencies(SmoduleMeshdata SmoduleCTX SmoduleFile SmoduleTopology SmoduleConstants)

add_library(SmoduleFormats STATIC ${moduleFormatsSource})
target_compile_options(SmoduleFormats PUBLIC "-DSPARSE")
# no dependencies

add_library(SmoduleIO STATIC ${moduleIOSource})
target_compile_options(SmoduleIO PUBLIC "-DSPARSE")
# no dependencies

add_library(SmoduleLabel STATIC ${moduleLabelSource})
target_compile_options(SmoduleLabel PUBLIC "-DSPARSE")
# no dependencies

add_library(SmoduleMaterials STATIC ${moduleMaterialsSource})
target_compile_options(SmoduleMaterials PUBLIC "-DSPARSE")
add_dependencies(SmoduleMaterials SmoduleMeshdata SmoduleDebug SmoduleCTX)

add_library(SmoduleModelDefinition STATIC ${moduleModelDefinitionSource})
target_compile_options(SmoduleModelDefinition PUBLIC "-DSPARSE")
# no dependencies

add_library(SmoduleModelData STATIC ${moduleModelDataSource})
target_compile_options(SmoduleModelData PUBLIC "-DSPARSE")
add_dependencies(SmoduleModelData SmoduleMeshdata SmoduleModelDefinition SmoduleDebug SmoduleTopology)


add_library(SmodulePETSc STATIC ${modulePETScSource})
target_compile_options(SmodulePETSc PUBLIC "-DSPARSE")
add_dependencies(SmodulePETSc SmoduleIO SmoduleMeshdata SmoduleFile SmoduleCTX SmoduleModelDefinition)

add_library(SmoduleTimeStep STATIC ${moduleTimeStepSource})
target_compile_options(SmoduleTimeStep PUBLIC "-DSPARSE")
add_dependencies(SmoduleTimeStep SmoduleDebug SmoduleCTX SmoduleConstants SmoduleModelDefinition)
target_link_libraries(SmoduleTimeStep libMisc)

add_library(SmoduleExtrema STATIC ${moduleExtremaSource})
target_compile_options(SmoduleExtrema PUBLIC "-DSPARSE")
# no dependencies

add_library(SmoduleFriction STATIC ${moduleFrictionSource})
target_compile_options(SmoduleFriction PUBLIC "-DSPARSE")
add_dependencies(SmoduleFriction SmoduleModelData)
target_link_libraries(SmoduleFriction ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

add_library(SmoduleTest STATIC ${moduleTestSource})
target_compile_options(SmoduleTest PUBLIC "-DSPARSE")
add_dependencies(SmoduleTest SmoduleFriction)

set_target_properties(${allSparseModules} PROPERTIES Fortran_MODULE_DIRECTORY ${sparseModDir})


###### other dependencies of main targets

# this sequence is the one used for linking.
# NOT the dependency list above.
# So do not change the sequence of these modules.


add_library(libModulesSparse STATIC ${moduleConstantsSource})
add_dependencies(libModulesSparse ${allSparseModules})
target_link_libraries(libModulesSparse ${allSparseModules})
target_link_libraries(libModulesSparse $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib/libpetsc.a ${METIS_LIBRARIES} -ldl ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
target_link_libraries(libModulesSparse libMisc)
target_compile_options(libModulesSparse PUBLIC "-DSPARSE")

add_library(libModulesNonSparse STATIC ${moduleConstantsSource})
target_link_libraries(libModulesNonSparse ${allNonSparseModules})
target_link_libraries(libModulesNonSparse libMisc)

add_library(gTectonGeneralSparse ${libGtectonGeneralSource})
target_include_directories(gTectonGeneralSparse PUBLIC ${sparseModDir})
add_dependencies(gTectonGeneralSparse libMisc libModulesSparse)
target_compile_options(gTectonGeneralSparse PUBLIC "-DSPARSE")
target_link_libraries(gTectonGeneralSparse libMisc)


add_library(gTectonGeneralNonSparse ${libGtectonGeneralSource})
target_include_directories(gTectonGeneralNonSparse PUBLIC ${nonSparseModDir})
add_dependencies(gTectonGeneralNonSparse libMisc libModulesNonSparse)
target_link_libraries(gTectonGeneralNonSparse libMisc)
add_dependencies(gTectonGeneralNonSparse triangle)

# plot libraries
add_library(libplnplt STATIC ${plnpltSource})
target_include_directories(libplnplt PUBLIC ${nonSparseModDir})
target_link_libraries(libplnplt libModulesNonSparse)
add_dependencies(libplnplt libModulesNonSparse)
#add_dependencies(libplnplt triangle)

add_library(libopnplt STATIC ${opnpltSource})
target_include_directories(libopnplt PUBLIC ${nonSparseModDir})
add_dependencies(libopnplt libModulesNonSparse)
target_link_libraries(libopnplt libModulesNonSparse)
#add_dependencies(libopnplt triangle)

add_library(libaxiplt STATIC ${axipltSource})
target_include_directories(libaxiplt PUBLIC ${nonSparseModDir})
add_dependencies(libaxiplt libModulesNonSparse)
target_link_libraries(libaxiplt libModulesNonSparse)
#add_dependencies(libaxiplt triangle)

add_library(libMisc STATIC ${libMiscSource})

add_library(plot2DGeneral ${plt2DSource})
add_dependencies(plot2DGeneral libModulesNonSparse)
target_include_directories(plot2DGeneral PUBLIC ${nonSparseModDir})

add_library(plot3DGeneral ${plt3DSource})
add_dependencies(plot3DGeneral libModulesNonSparse)
target_include_directories(plot3DGeneral PUBLIC ${nonSparseModDir})


###### tools:

add_executable(picknps ${picknpsSource})
target_link_libraries(picknps m)

# triangle splits between the SPARSE and non-SPARSE sections of the code.
add_executable(triangle ${triangleSource})
target_link_libraries(triangle m)

add_executable(showme ${showmeSource})
target_link_libraries(showme ${X11_X11_LIB})
# note that the variable X11_LIBRARIES contains other unneeded stuff for showme, 
# that causes an error with -luuid I could not resolve.
# X11_X11_LIB contains just libX11.so


add_executable(elmside ${elmsideSource})

add_executable(setarea ${setAreaSource})
add_dependencies(setarea libMisc)
target_link_libraries(setarea libMisc)

add_executable(polysort ${polysortSource})
add_dependencies(polysort libMisc)
target_link_libraries(polysort libMisc)

add_executable(mergefiles ${mergefilesSource})
target_include_directories(mergefiles PUBLIC ${nonSparseModDir})
add_dependencies(mergefiles libModulesNonSparse libMisc)
target_link_libraries(mergefiles libModulesNonSparse libMisc)

add_executable(gmsh2tecton ${gmsh2tectonSource})
target_include_directories(gmsh2tecton PUBLIC ${nonSparseModDir})
add_dependencies(gmsh2tecton libMisc)
target_link_libraries(gmsh2tecton libMisc)

add_executable(partition2paraview ${partition2paraviewSource})

add_executable(partition ${partitionSource})
target_link_libraries(partition ${METIS_LIBRARIES} m)
target_include_directories(partition PUBLIC ${METIS_INCLUDE_DIR})

add_executable(tri2fe ${tri2feSource} ${libMiscSource})
add_dependencies(tri2fe libMisc)
target_link_libraries(tri2fe libMisc)

add_executable(domainUnity ${domainUnitySource})

add_executable(PartitionedGmshToTecton ${PartitionedGmshToTectonSource})

add_executable(setprem ${setpremSource})

add_executable(getnormals ${getnormalsSource})

add_executable(partitionInterfaceCount ${partitionInterfaceCountSource})

add_executable(relvel ${relvelSource})
add_dependencies(relvel libMisc)
target_link_libraries(relvel libMisc)

add_executable(regrid ${regridSource} ${plnlibSource})
add_dependencies(regrid libModulesSparse libMisc)
target_link_libraries(regrid gTectonGeneralSparse libModulesSparse libMisc)


#add_dependencies(regrid libModulesNonSparse)
#target_link_libraries(regrid libModulesNonSparse tools)


###### main executables:

# pln

add_executable(pln ${libGtectonMainSource} ${plnlibSource})
add_dependencies(pln libModulesSparse libMisc)
target_link_libraries(pln gTectonGeneralSparse libModulesSparse libMisc ${CMAKE_DL_LIBS})

add_executable(plnplt ${plnlibSource})
add_dependencies(plnplt libModulesNonSparse libMisc plot2DGeneral gTectonGeneralNonSparse)
target_link_libraries(plnplt plot2DGeneral libplnplt gTectonGeneralNonSparse libModulesNonSparse libMisc)
add_dependencies(plnplt triangle)

# axi

add_executable(axi ${libGtectonMainSource} ${axilibSource})
add_dependencies(axi libModulesSparse libMisc)
target_link_libraries(axi gTectonGeneralSparse libModulesSparse libMisc ${CMAKE_DL_LIBS})

add_executable(axiplt ${axilibSource})
add_dependencies(axiplt libModulesNonSparse libMisc plot2DGeneral gTectonGeneralNonSparse)
target_link_libraries(axiplt plot2DGeneral libaxiplt gTectonGeneralNonSparse libModulesNonSparse libMisc)
add_dependencies(axiplt triangle)

# opn

add_executable(opn ${libGtectonMainSource} ${opnlibSource})
add_dependencies(opn libModulesSparse libMisc)
target_compile_options(opn PUBLIC "-l${CMAKE_DL_LIBS}")
target_link_libraries(opn gTectonGeneralSparse libModulesSparse libMisc ${CMAKE_DL_LIBS})

add_executable(opnplt ${opnlibSource})
add_dependencies(opnplt libModulesNonSparse libMisc plot2DGeneral gTectonGeneralNonSparse)
target_link_libraries(opnplt plot2DGeneral libopnplt gTectonGeneralNonSparse libModulesNonSparse libMisc)
add_dependencies(opnplt triangle)

#f3d

add_executable(f3d ${libGtectonMainSource} ${f3dlibSource})
add_dependencies(f3d libModulesSparse libMisc)
target_link_libraries(f3d gTectonGeneralSparse libModulesSparse libMisc ${CMAKE_DL_LIBS})

add_executable(plt3d ${f3dlibSource})
add_dependencies(plt3d libModulesNonSparse libMisc plot3DGeneral gTectonGeneralNonSparse)
target_link_libraries(plt3d plot3DGeneral gTectonGeneralNonSparse libModulesNonSparse libMisc)
add_dependencies(plt3d triangle)

