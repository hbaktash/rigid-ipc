# if(TARGET Boost::boost)
#     return()
# endif()

# message(STATUS "Third-party: creating targets 'Boost::boost'")

# # Add boost lib sources
# set(BOOST_INCLUDE_LIBRARIES numeric/interval)
# set(BOOST_ENABLE_CMAKE ON)

# # Download and extract the boost library from GitHub
# message(STATUS "Downloading and extracting boost library sources. This will take some time...")
# include(FetchContent)
# Set(FETCHCONTENT_QUIET FALSE) # Needed to print downloading progress
# FetchContent_Declare(
#     Boost
#     URL https://github.com/boostorg/boost/releases/download/boost-1.84.0/boost-1.84.0.7z # downloading a zip release speeds up the download
#     USES_TERMINAL_DOWNLOAD TRUE 
#     GIT_PROGRESS TRUE   
#     DOWNLOAD_NO_EXTRACT FALSE
# )
# FetchContent_MakeAvailable(Boost)

# # Add the boost library to the project
# add_library(Boost::boost INTERFACE IMPORTED)
# target_include_directories(Boost::boost INTERFACE ${boost_SOURCE_DIR})
# target_compile_definitions(Boost::boost INTERFACE BOOST_ALL_NO_LIB)
# target_link_libraries(Boost::boost INTERFACE Boost::numeric_interval Boost::core)

if(TARGET Boost::boost)
    return()
endif()

message(STATUS "Third-party: creating targets 'Boost::boost'")

include(FetchContent)
FetchContent_Declare(
    boost-cmake
    GIT_REPOSITORY https://github.com/Orphis/boost-cmake.git
    GIT_TAG 7f97a08b64bd5d2e53e932ddf80c40544cf45edf
)

set(PREVIOUS_CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
set(OLD_CMAKE_POSITION_INDEPENDENT_CODE ${CMAKE_POSITION_INDEPENDENT_CODE})
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# This guy will download boost using FetchContent
FetchContent_GetProperties(boost-cmake)
if(NOT boost-cmake_POPULATED)
    FetchContent_Populate(boost-cmake)
    # File lcid.cpp from Boost_locale.cpp doesn't compile on MSVC, so we exclude them from the default
    # targets being built by the project (only targets explicitly used by other targets will be built).
    add_subdirectory(${boost-cmake_SOURCE_DIR} ${boost-cmake_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

set(CMAKE_POSITION_INDEPENDENT_CODE ${OLD_CMAKE_POSITION_INDEPENDENT_CODE})
set(CMAKE_CXX_FLAGS "${PREVIOUS_CMAKE_CXX_FLAGS}")