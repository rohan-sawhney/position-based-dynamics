#ifndef CLOTH_GENERATOR_H
#define CLOTH_GENERATOR_H

#include "Types.h"

class ClothGenerator {
public:
    // constructor
    ClothGenerator();
    
    // generate cloth
    bool generate(const double& extent, const double& base, const int& intervals, Mesh& mesh);
};

#endif
