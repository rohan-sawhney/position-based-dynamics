#ifndef SIMULATOR_H
#define SIMULATOR_H

#include "Types.h"
#include "Mesh.h"
#include "ConstraintHandler.h"

class Simulator {
public:
    // constructor
    Simulator(std::vector<Mesh>& meshes0);
    
    // initialize
    void initialize();
    
    // step
    void step();
    
    // reset
    void reset();
    
private:
    // member variables
    std::vector<Mesh>& meshes;
    ConstraintHandler constraintHandler;
};

#endif
