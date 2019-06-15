# -- CMake build rules

# Boost stuff
set(Boost_USE_STATIC_LIBS ON CACHE BOOL "Use statically-linked boost libraries")
# '- Note: By default set Boost_USE_STATIC_LIBS to ON (to link them statically
#    to avoid possible problems during execution of the application (not finding
#    the .so files)), but let the user change it (useful for example to use
#    Fedora native boost libraries)
set(Boost_USE_MULTITHREADED ON)
# find_package(Boost) # <-- Use this simple line for using only headers
find_package(Boost
  COMPONENTS
    # lexicographically sorted
    date_time
    filesystem
    program_options
    system
)
if (Boost_FOUND)
  message(STATUS "Boost found:")
  include_directories(${Boost_INCLUDE_DIRS})
  message(STATUS "  Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
  message(STATUS "  Boost_LIBRARIES: ${Boost_LIBRARIES}
  Note: Boost_LIBRARIES can be empty if no COMPONENTS are used.") # Used below
else ()
  message(FATAL_ERROR "Boost not found - Hint(?): ENV{BOOST_ROOT} is $ENV{BOOST_ROOT}")
endif ()

# ITK stuff
FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

# Depending of the ITK version add -fpermissive. Delete when no longer needed
if (${ITK_VERSION} VERSION_LESS "4.9" AND "${CMAKE_SYSTEM_NAME}" MATCHES "Linux")
  message(STATUS "ITK < 4.9, (version ${ITK_VERSION}) on Linux found, adding -fpermissive")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive")
endif()

if (ITK_FOUND)
  message(STATUS "ITK found: ${ITK_USE_FILE}")
else ()
  message(FATAL_ERROR "ITK not found")
endif ()


set(SRC_LIST
  aux-raw-compiler-warnings-off++begin.h aux-raw-compiler-warnings-off++end.h
  main.cpp
)

add_executable(${PROJECT_NAME} ${SRC_LIST})

# Link the executable to the required libraries

target_link_libraries(${PROJECT_NAME}
  ${Boost_LIBRARIES}
  ${ITK_LIBRARIES}
)

set(TARGET_NAME ${PROJECT_NAME}Targets)
install(TARGETS ${PROJECT_NAME} # EXPORT ${TARGET_NAME}
# LIBRARY DESTINATION lib
# ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
# INCLUDES DESTINATION include
)
