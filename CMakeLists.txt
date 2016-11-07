cmake_minimum_required(VERSION 3.6)
project(staci)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    Agelem.cpp
    Agelem.h
    AnyOption.cpp
    AnyOption.h
    BukoMutargy.cpp
    BukoMutargy.h
    Csatorna.cpp
    Csatorna.h
    Cso.cpp
    Cso.h
    Csomopont.cpp
    Csomopont.h
    data_io.cpp
    data_io.h
    JelleggorbesFojtas.cpp
    JelleggorbesFojtas.h
    KonstNyomas.cpp
    KonstNyomas.h
    lubksb.cpp
    ludcmp.cpp
    nr.h
    nrtypes.h
    nrtypes_nr.h
    nrutil.h
    nrutil_nr.h
    Staci.cpp
    Staci.h
    staci_calibrate.cpp
    staci_main.cpp
    StaciException.h
    Szivattyu.cpp
    Szivattyu.h
    Vegakna.cpp
    Vegakna.h
    VisszacsapoSzelep.cpp
    VisszacsapoSzelep.h
    xmlParser.cpp
    xmlParser.h staci_calibrate.cpp.backup1.cpp staci_split_graph.cpp)

add_executable(staci ${SOURCE_FILES})