CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(plsm)

FIND_PACKAGE(CUDA REQUIRED)
#FIND_PACKAGE(MPI REQUIRED)

INCLUDE(FindCUDA)

INCLUDE_DIRECTORIES(/usr/local/cuda-8.0/include /home/m2a02/NVIDIA_CUDA-8.0_Samples/common/inc)

#added by wasif
INCLUDE_DIRECTORIES(src/)
INCLUDE_DIRECTORIES(cmdline/)

FILE(GLOB CMDLINE_SRCS "cmdline/*.cpp")
#add_executable(plsm ${CMDLINE_SRCS})


include_directories(pure-include/)
#FILE(GLOB PLSMMODSRCS Array.cc plsm.cu)

FILE(GLOB SOURCES "*.cu" "*.cpp" "*.c" "*.cc" "*.h")
CUDA_ADD_EXECUTABLE(plsm ${SOURCES} ${CMDLINE_SRCS})
CUDA_ADD_LIBRARY(idiap-plsm SHARED ${SOURCES} ${CMDLINE_SRCS} )
#add_library(idiap-plsm SHARED ${SOURCES} ${CMDLINE_SRCS})

LIST(APPEND CMAKE_CXX_FLAGS "-std=c++0x -O3 -ffast-math -Wall")

LIST(APPEND CUDA_NVCC_FLAGS --compiler-options -fno-strict-aliasing -lineinfo -use_fast_math -Xptxas -dlcm=cg)
LIST(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_20,code=sm_20)
LIST(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_30,code=sm_30)
LIST(APPEND CUDA_NVCC_FLAGS -gencode arch=compute_35,code=sm_35)
IF(UNIX)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -D_FORTIFY_SOURCE=2 -march=core2 -mcx16 -msahf -mpopcnt -msse4.2 --param l1-cache-size=32 --param l1-cache-size=32 --param l1-cache-line-size=64 --param l2-cache-size=12288 -ffast-math -funsafe-loop-optimizations -fmodulo-sched -fmodulo-sched-allow-regmoves -fgcse-after-reload -O3")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -D_FORTIFY_SOURCE=2 -march=core2 -mcx16 -msahf -mpopcnt -msse4.2 --param l1-cache-size=32 --param l1-cache-size=32 --param l1-cache-line-size=64 --param l2-cache-size=12288 -ffast-math -funsafe-loop-optimizations -fmodulo-sched -fmodulo-sched-allow-regmoves -fgcse-after-reload -fipa-matrix-reorg -floop-parallelize-all -ftree-parallelize-loops=4 -floop-strip-mine -floop-block -ftree-loop-linear -ftree-loop-distribution -fipa-pta -O3")
# -fwhole-program  -flto
ENDIF(UNIX)
TARGET_LINK_LIBRARIES(plsm /usr/local/cuda/lib64/libcudart.so ${MPI_LIBRARIES})
TARGET_LINK_LIBRARIES(idiap-plsm /usr/local/cuda/lib64/libcudart.so ${MPI_LIBRARIES})
