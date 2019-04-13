cmake_minimum_required(VERSION 3.13)
project(NGS_Demo)

set(CMAKE_CXX_STANDARD 14)

add_executable(NGS_Demo code/htslib_sample.cpp code/correct/CDBG.h code/correct/CDBG.cpp code/correct/util/ThreadPool.h code/correct/util/Semaphore.h code/correct/RMtoDBG.h code/correct/RMtoDBG.cpp code/correct/util/Usage.hpp code/correct/KmerStat.cpp code/correct/KmerStat.h)

