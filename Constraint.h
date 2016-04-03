#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#include "Types.h"

class Constraint {
public:
    // constructor
    Constraint(const int& n0, const double& k0);
    
    // destructor
    virtual ~Constraint();
    
    // cardinality
    int n;
    
    // indices
    std::vector<VertexIter> vs;
    
    // stiffness parameter
    double k;
    
    // solves position constraint
    virtual void solve() = 0;
    
    // update velocity
    virtual void updateVelocity() = 0;
};

class StretchingConstraint: public Constraint {
public:
    // constructor
    StretchingConstraint(VertexIter v1, VertexIter v2, const double& k0);
    
    // destructor
    ~StretchingConstraint();
    
    // rest length
    double l;
    
    // solves stretching constraint
    void solve();
    
    // update velocity
    void updateVelocity();
};

class BendingConstraint: public Constraint {
public:
    // constructor
    BendingConstraint(VertexIter v1, VertexIter v2,
                      VertexIter v3, VertexIter v4,
                      const double& k0);
    
    // destructor
    ~BendingConstraint();
    
    // rest angle
    double theta;
    
    // solves bending constraint
    void solve();
    
    // update velocity
    void updateVelocity();
};

class CollisionConstraint: public Constraint {
public:
    // constructor
    CollisionConstraint(VertexIter v1, const Eigen::Vector3d& q0, const Eigen::Vector3d& normal0);
    
    // destructor
    ~CollisionConstraint();
    
    // collision point
    Eigen::Vector3d q;
    
    // normal
    Eigen::Vector3d normal;
    
    // flag
    bool didSolve;
    
    // solves collision constraint
    void solve();
    
    // update velocity
    void updateVelocity();
};

#endif 
