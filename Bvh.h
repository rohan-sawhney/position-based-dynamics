#ifndef BVH_H
#define BVH_H
#define RAY_INTERSECTION 0
#define NEAREST_POINT 1

#include "Types.h"
#include "BoundingBox.h"

class Node {
public:
    // member variables
    BoundingBox boundingBox;
    int startId, range, rightOffset;
};

class Bvh {
public:
    // constructor
    Bvh(const int& leafSize0 = 1);
    
    // builds the bvh
    void build(Mesh *mesh0);
    
    // returns face index
    int getIntersection(const int& mode, double& hit, Eigen::Vector3d& q,
                        const Eigen::Vector3d& o, const Eigen::Vector3d& d = Eigen::Vector3d::Zero()) const;
    
private:
    int nodeCount, leafCount, leafSize;
    std::vector<Node> flatTree;
    Mesh *mesh;
};

#endif
