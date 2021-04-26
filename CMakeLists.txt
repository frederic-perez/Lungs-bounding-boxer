cmake_minimum_required(VERSION 3.0)

set(PROJECT_NAME "Lungs-bounding-boxer")
set(PROJECT_TYPE "CXX")
project(${PROJECT_NAME} ${PROJECT_TYPE})

set(CMAKE_VERBOSE_MAKEFILE TRUE)

message(STATUS "CMAKE_C_PLATFORM_ID is ${CMAKE_C_PLATFORM_ID}")
message(STATUS "CMAKE_SYSTEM_NAME is ${CMAKE_SYSTEM_NAME}")

if("${CMAKE_SYSTEM_NAME}" MATCHES "Windows")
  set(MyUsername $ENV{username})
elseif("${CMAKE_SYSTEM_NAME}" MATCHES "Linux")
  set(MyUsername $ENV{USER})
elseif("${CMAKE_SYSTEM_NAME}" MATCHES "Darwin")
  set(MyUsername $ENV{USER})
else()
  set(MyUsername "Unknown")
endif()
message(STATUS "Your username (MyUsername variable) is ${MyUsername}")
set(def_CMAKE_MYUSERNAME CMAKE_MYUSERNAME=${MyUsername})
message(STATUS "> def_CMAKE_MYUSERNAME is ${def_CMAKE_MYUSERNAME}")

cmake_host_system_information(RESULT MyHostname QUERY HOSTNAME)
message(STATUS "Your host name (MyHostname variable) is ${MyHostname}")
set(def_CMAKE_MYHOSTNAME CMAKE_MYHOSTNAME=${MyHostname})
message(STATUS "> def_CMAKE_MYHOSTNAME is ${def_CMAKE_MYHOSTNAME}")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if("${CMAKE_SYSTEM_NAME}" MATCHES "Windows")
  if(MSVC)
    message(STATUS "Entering Windows/MSVC branch.")
    string(CONCAT CMAKE_CXX_FLAGS
      "${CMAKE_CXX_FLAGS} -D_SCL_SECURE_NO_WARNINGS -D_CRT_SECURE_NO_WARNINGS "
      "/W4 /WX /EHa /wd4127 /wd4351 "
      "/D${def_CMAKE_MYUSERNAME} /D${def_CMAKE_MYHOSTNAME}")
    # '- http://msdn.microsoft.com/en-US/library/thxezb7y.aspx:
    #    /wd n  Disables the compiler warning that is specified in n.
    #    For example, /wd4326 disables compiler warning C4351.
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /FORCE:MULTIPLE")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /LARGEADDRESSAWARE")
    message(STATUS "CMAKE_EXE_LINKER_FLAGS are: " ${CMAKE_EXE_LINKER_FLAGS})
  else()
    message(ERROR "Uncontemplated case.")
  endif()
elseif("${CMAKE_SYSTEM_NAME}" MATCHES "Linux")
  if(CMAKE_COMPILER_IS_GNUCXX)
    message(STATUS "Entering Linux/CMAKE_COMPILER_IS_GNUCXX branch.")
# TODO: Use string(CONCAT CMAKE_CXX_FLAGS after updating cmake-gui!
    set(CMAKE_CXX_FLAGS
      "${CMAKE_CXX_FLAGS} -Wall -std=c++11 -D${def_CMAKE_MYUSERNAME} -D${def_CMAKE_MYHOSTNAME}")
  else()
    message(STATUS "Entering Linux/else branch. We assume clang++.")
# TODO: Use string(CONCAT CMAKE_CXX_FLAGS after updating cmake-gui!
    set(CMAKE_CXX_FLAGS
      "${CMAKE_CXX_FLAGS} -Wall -std=c++11 -D ${def_CMAKE_MYUSERNAME} -D ${def_CMAKE_MYHOSTNAME}")
  endif()
elseif("${CMAKE_SYSTEM_NAME}" MATCHES "Darwin")
  message(STATUS "Entering Darwin branch. We assume clang++.")
  string(CONCAT CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} -Wall -std=c++11 "
    "-D${def_CMAKE_MYUSERNAME} -D${def_CMAKE_MYHOSTNAME}")
else()
  message(FATAL_ERROR "Unknown system")
endif()

message(STATUS "CMAKE_CXX_FLAGS are:" ${CMAKE_CXX_FLAGS})

add_subdirectory(src)
#add_subdirectory(data)

#include(CTest)
