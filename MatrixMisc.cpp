#include "MatrixMisc.h"

namespace Mat {


int countOccurrences(const std::string& str, char ch) {
    int count = 0;
    for (char c : str) {
        if (c == ch) {
            count++;
        }
    }
    return count;
}

}