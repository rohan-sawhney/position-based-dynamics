#ifndef CONSTRAINT_HANDLER_H
#define CONSTRAINT_HANDLER_H

#include "Types.h"
class Constraint;

class ConstraintHandler {
public:
    // constructor
    ConstraintHandler();
    
    // destructor
    ~ConstraintHandler();
    
    // generate constraints
    void generateConstraints(const double& kStretch, const double& kBend, Mesh& mesh);
    
    // generate collisions
    void generateCollisions(std::vector<Mesh>& meshes);
    
    // project constraints
    void projectConstraints();
    
    // project collisions
    void updateVelocities();
    
    // clear constraints
    void clear();
    
private:
    // clear collisions
    void clearCollisions();
    
    // member variables
    std::vector<Constraint *> constraints;
    int count;
};

#endif
