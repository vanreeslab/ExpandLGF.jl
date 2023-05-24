################################################################################
# @copyright Copyright © MIT 2023
# 
# ExpandLGF++ is a header-only library for evaluating Lattice Green's Functions
# in domains with one unbounded direction.
# 
################################################################################
# C++ compiler
CXX := clang++

# optional flags
CXXFLAGS ?= -O3 -g -std=c++17

# GTest
GTEST_DIR ?= /soft/googletest
GTEST_INC ?= ${GTEST_DIR}/include
GTEST_LIB ?= ${GTEST_DIR}/lib
GTEST_LIBNAME ?= -lgtest
GTEST_MAIN ?= -lgtest_main

################################################################################
# sources
SRC_DIR := .
SRC_HPP := ${SRC_DIR}/lgf_one_unbounded.hpp

# test
TEST_DIR := .
TEST_HPP := ${TEST_DIR}/test.hpp
TEST_CPP := ${TEST_DIR}/test.cpp
TEST_NAME := lgf_test

# includes
INC := -I$(SRC_DIR) -I$(TEST_DIR) 
TEST_INC := -I$(GTEST_INC)

# linking
LIB := -pthread 
TEST_LIB := -L$(GTEST_LIB) $(GTEST_MAIN) $(GTEST_LIBNAME) -Wl,-rpath,$(GTEST_LIB) 

################################################################################
# Targets
test : $(SRC_HPP) $(TEST_HPP) $(TEST_CPP)
	$(CXX) $(CXXFLAGS) $(INC) $(TEST_INC) $(LIB) $(TEST_LIB) $(TEST_CPP) -o $(TEST_NAME)

clean:
	rm -f lgf_test
