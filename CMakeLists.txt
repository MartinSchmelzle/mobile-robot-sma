
cmake_minimum_required(VERSION 3.14.5)  # Minimum CMake version required


project(slimemould_cpp)  # Replace 'MyProject' with your project name


add_executable(slime_alg_exe CALL-ALG.cpp main_alg.cpp F00.cpp)

#the following are compiler commands that improve performance
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
#set the following ones when you know the code works
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-generate")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
#this one disables assertations
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DNDEBUG")

