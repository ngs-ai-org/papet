# compiler options
add_compile_options(-std=c++2a)
add_compile_options(-O3)
add_compile_options(-Wall)
add_compile_options(-Wextra)
add_compile_options(-Werror)
add_compile_options(-Wfatal-errors)
add_compile_options(-pedantic)
add_compile_options(-pg)

# include file locations
include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories(${INCLUDE_DIRECTORIES})
include_directories(${Boost_INCLUDE_DIRS})

# linking to external lib
# Boost, pbcopper, pbbam, gtest are there
# link_directories("/usr/local/lib")
link_directories(${LINK_DIRECTORY})

# list of src files for lib
set(FILES_CPP 
    "applications/ApplicationInterface.cpp"
    "applications/ApplicationKinetics.cpp"
    "applications/ApplicationKineticsKmer.cpp"
    "applications/ApplicationKineticsWig.cpp"
    "applications/ApplicationModelKinetic.cpp"
    "applications/ApplicationModelKineticTxt.cpp"
    "applications/ApplicationModelSequence.cpp"
    "applications/ApplicationModelSequenceTxt.cpp"
    "applications/ApplicationPredict.cpp"
    "applications/ApplicationPapet.cpp")


# make install, as set up by cmake, will erase the 
# RUNPATH from the ELF header of the executable.
# After calling make install, the executable won't be able
# to find the shared libraries.
# This forces cmake to set RUNPATH after make install to 
# the given value in the executable ELF header.
set(CMAKE_INSTALL_RPATH ${LINK_DIRECTORY})

# compile executable
## papet
set(EXE_PAPET "papet")
add_executable(${EXE_PAPET} 
               ${FILES_CPP})
target_link_libraries(${EXE_PAPET} ngsaipp
                                   pthread
                                   pbbam
                                   pbcopper
                                   boost_program_options
                                   boost_serialization)
set_target_properties(${EXE_PAPET} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${papet_SOURCE_DIR}/bin")

install(TARGETS ${EXE_PAPET}
        RUNTIME DESTINATION ${INSTALL_DIRECTORY})