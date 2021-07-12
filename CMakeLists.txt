cmake_minimum_required(VERSION 3.19)
project(tng2 C)

set(CMAKE_C_STANDARD 11)

include_directories(.)

add_executable(tng
        tng.c
        bwlzh.c
        bwlzh.h
        bwt.c
        bwt.h
        coder.c
        coder.h
        dict.c
        dict.h
        fixpoint.c
        fixpoint.h
        huffman.c
        huffman.h
        huffmem.c
        lz77.c
        lz77.h
        merge_sort.c
        merge_sort.h
        mtf.c
        mtf.h
        my64bit.h
        rle.c
        rle.h
        tng_compress.c
        tng_compress.h
        vals16.c
        vals16.h
        warnmalloc.c
        warnmalloc.h
        widemuldiv.c
        widemuldiv.h
        xtc2.c
        xtc3.c)
target_link_libraries(tng m)
