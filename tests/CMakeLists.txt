cmake_minimum_required(VERSION 3.14)
project(CMTJ)

# GoogleTest requires at least C++11
set(CMAKE_CXX_STANDARD 20)

find_package (Eigen3 3.4 REQUIRED NO_MODULE)

######## FFTW DOWNLOADS ########

# download FFTW
# configure_file(downloadFindFFTW.cmake.in findFFTW-download/CMakeLists.txt)
# execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
#         RESULT_VARIABLE result
#         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-download )
# if(result)
#     message(FATAL_ERROR "CMake step for findFFTW failed: ${result}")
#     else()
#     message("CMake step for findFFTW completed (${result}).")
# endif()
# execute_process(COMMAND ${CMAKE_COMMAND} --build .
#         RESULT_VARIABLE result
#         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-download )
# if(result)
#     message(FATAL_ERROR "Build step for findFFTW failed: ${result}")
# endif()

# set(findFFTW_DIR ${CMAKE_CURRENT_BINARY_DIR}/findFFTW-src)

# set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${findFFTW_DIR}")
# find_package(FFTW REQUIRED)


include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)


enable_testing()

include(GoogleTest)

# source: https://cliutils.gitlab.io/modern-cmake/chapters/testing/googletest.html
macro(package_add_test TESTNAME)
    # create an exectuable in which the tests will be stored
    add_executable(${TESTNAME} ${ARGN}) # ${FFTW_INCLUDE_DIRS})
    # link the Google test infrastructure, mocking library, and a default main fuction to
    # the test executable.  Remove g_test_main if writing your own main function.
    target_link_libraries(${TESTNAME} gtest gmock gtest_main) #${FFTW_LIBRARIES} ${FFTW_DOUBLE_LIB} ${FFTW_DOUBLE_THREADS_LIB})
    # gtest_discover_tests replaces gtest_add_tests,
    # see https://cmake.org/cmake/help/v3.10/module/GoogleTest.html for more options to pass to it
    gtest_discover_tests(${TESTNAME}
        # set a working directory so your project root so that you can find test data via paths relative to the project root
        WORKING_DIRECTORY ${PROJECT_DIR}
        PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${PROJECT_DIR}"
    )
    set_target_properties(${TESTNAME} PROPERTIES FOLDER tests)
endmacro()

# package_add_test(test_result test_result.cpp)
# package_add_test(test_rng test_rng.cpp)
# package_add_test(test_fft test_fft.cpp)
# package_add_test(test_sot test_sot.cpp)
# package_add_test(test_vector test_vector.cpp)
# package_add_test(test_layer test_layer.cpp)
package_add_test(test_reservoir test_reservoir.cpp)
target_link_libraries(test_reservoir Eigen3::Eigen)
