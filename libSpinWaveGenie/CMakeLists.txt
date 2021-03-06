file(GLOB CONTAINERS_SRCS src/Containers/*.cpp)
file(GLOB GENIE_SRCS src/Genie/*.cpp)
file(GLOB INTERACTIONS_SRCS src/Interactions/*.cpp)
file(GLOB PLOT_SRCS src/Plot/*.cpp)
SOURCE_GROUP("Source Files\\Containers" FILES ${CONTAINERS_SRCS})
SOURCE_GROUP("Source Files\\Genie" FILES ${GENIE_SRCS})
SOURCE_GROUP("Source Files\\Interactions" FILES ${INTERACTIONS_SRCS})
SOURCE_GROUP("Source Files\\Plot" FILES ${PLOT_SRCS})
file(GLOB CONTAINERS_HDRS include/SpinWaveGenie/Containers/*.h)
file(GLOB GENIE_HDRS include/SpinWaveGenie/Genie/*.h)
file(GLOB INTERACTIONS_HDRS include/SpinWaveGenie/Interactions/*.h)
file(GLOB PLOT_HDRS include/SpinWaveGenie/Plot/*.h)
SOURCE_GROUP("Header Files\\Containers" FILES ${CONTAINERS_HDRS})
SOURCE_GROUP("Header Files\\Genie" FILES ${GENIE_HDRS})
SOURCE_GROUP("Header Files\\Interactions" FILES ${INTERACTIONS_HDRS})
SOURCE_GROUP("Header Files\\Plot" FILES ${PLOT_HDRS})
set(LIBRARY_SRCS ${EXTERNAL_SRCS} ${EXTERN_SRCS} ${CELL_SRCS} ${CONTAINERS_SRCS} ${GENIE_SRCS} ${INTERACTIONS_SRCS} ${PLOT_SRCS})
set(LIBRARY_HDRS ${EXTERNAL_HDRS} ${EXTERN_HDRS} ${CELL_HDRS} ${CONTAINERS_HDRS} ${GENIE_HDRS} ${INTERACTIONS_HDRS} ${PLOT_HDRS})

if (COVERALLS)
    foreach(loop_var ${LIBRARY_SRCS} ${LIBRARY_HDRS})
      set(COVERAGE_SRCS "${COVERAGE_SRCS}\n${loop_var}")
    endforeach(loop_var)
    string(SUBSTRING ${COVERAGE_SRCS} 1 -1 SRCS_FILE)
    set( SRCS_FILENAME "${CMAKE_CURRENT_BINARY_DIR}/sources.txt")
    file(WRITE ${SRCS_FILENAME} ${SRCS_FILE})
    # Create the coveralls target.
    coveralls_setup(
        ${SRCS_FILENAME} # The source files.
        ${COVERALLS_UPLOAD}    # If we should upload.
        ${COVERALLS_TRAVISCI}  # building on travis-ci
        )
endif()

include(GenerateExportHeader)
add_library(SpinWaveGenie SHARED ${LIBRARY_SRCS} ${LIBRARY_HDRS})

#borrow LTO settings from pybind11
if (PYBIND11_LTO_CXX_FLAGS)
  target_compile_options(SpinWaveGenie PRIVATE "$<$<NOT:$<CONFIG:Debug>>:${PYBIND11_LTO_CXX_FLAGS}>")
endif()
if (PYBIND11_LTO_LINKER_FLAGS)
  target_link_libraries(SpinWaveGenie PRIVATE "$<$<NOT:$<CONFIG:Debug>>:${PYBIND11_LTO_LINKER_FLAGS}>")
endif()

generate_export_header(SpinWaveGenie EXPORT_FILE_NAME SpinWaveGenie/Export.h)

#We need need the header files both to build SpinWaveGenie and for the interface.
target_include_directories( SpinWaveGenie PRIVATE include ${CMAKE_CURRENT_BINARY_DIR})               
set_property( TARGET SpinWaveGenie PROPERTY
              INTERFACE_INCLUDE_DIRECTORIES
              "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include;${CMAKE_CURRENT_BINARY_DIR}>"
              "$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>"
)
set_property( TARGET SpinWaveGenie PROPERTY CXX_VISIBILITY_PRESET hidden)
set_property( TARGET SpinWaveGenie PROPERTY VISIBILITY_INLINES_HIDDEN 1)

target_include_directories( SpinWaveGenie SYSTEM INTERFACE PUBLIC ${Boost_INCLUDE_DIR} )
target_include_directories( SpinWaveGenie PRIVATE src)
target_include_directories( SpinWaveGenie SYSTEM PRIVATE ${TBB_INCLUDE_DIRS})
target_link_libraries( SpinWaveGenie PUBLIC Eigen3::Eigen
                                     PRIVATE ${TBB_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

if(BUILD_TESTING)
    add_subdirectory(test)
endif (BUILD_TESTING)

include(GNUInstallDirs)

install(TARGETS SpinWaveGenie EXPORT SpinWaveGenieTargets 
                              ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
                              LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

#create package configuration files
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/SpinWaveGenie/SpinWaveGenieConfigVersion.cmake"
    VERSION ${SPINWAVEGENIE_VERSION}
    COMPATIBILITY AnyNewerVersion
)

export(EXPORT SpinWaveGenieTargets
    FILE "${CMAKE_CURRENT_BINARY_DIR}/SpinWaveGenie/SpinWaveGenieTargets.cmake"
    NAMESPACE SpinWaveGenie::
)

configure_file(SpinWaveGenieConfig.cmake
    "${CMAKE_CURRENT_BINARY_DIR}/SpinWaveGenie/SpinWaveGenieConfig.cmake"
    COPYONLY
)

#install package configuration files
set(ConfigPackageLocation ${CMAKE_INSTALL_LIBDIR}/cmake/SpinWaveGenie)
install(EXPORT SpinWaveGenieTargets
    FILE SpinWaveGenieTargets.cmake
    DESTINATION ${ConfigPackageLocation}
)

install(
    FILES
    "${CMAKE_CURRENT_BINARY_DIR}/SpinWaveGenie/SpinWaveGenieConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/SpinWaveGenie/SpinWaveGenieConfigVersion.cmake"
    DESTINATION ${ConfigPackageLocation}
)

