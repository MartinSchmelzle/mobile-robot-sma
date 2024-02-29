
cmake_minimum_required(VERSION 3.14.5)  # Minimum CMake version required


project(slimemould_cpp)  # Replace 'MyProject' with your project name


add_executable(slime_alg_exe examples/example-2.cpp src/main_alg.cpp src/F00.cpp src/F00_functions.cpp src/custom_math_operations.cpp src/get_Energy_v2.cpp src/getData.cpp src/SMA.cpp src/processpath.cpp)

#the following are compiler commands that improve performance
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
#set the following ones when you know the code works
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-generate")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
#this one disables assertations
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DNDEBUG")

