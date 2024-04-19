#pragma once
#include <iostream>
#include <chrono>
#include <functional> 

// if defined show some debug info
// #define MATRIX_DEBUG

typedef struct Dimensions {
	int n,m;

    bool operator==(const Dimensions& other) const {
        return n == other.n && m == other.m;
    }

    bool operator!=(const Dimensions& other) const {
        return !(*this == other);
    }
} Dimensions;

