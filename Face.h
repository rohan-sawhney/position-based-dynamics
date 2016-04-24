#ifndef FACE_H
#define FACE_H

#include "Types.h"
class BoundingBox;

class Face {
public:
    // one of the halfedges associated with this face
    HalfEdgeIter he;
    
    // id between 0 and |F|-1
    int index;
    
    // checks if this face lies on boundary
    bool isBoundary() const;
    
    // returns normal to face
    Eigen::Vector3d normal() const;
    
    // returns centroid
    Eigen::Vector3d centroid() const;
    
    // returns face area
    double area() const;
    
    // computes the bounding box of the face
    BoundingBox boundingBox() const;
    
    // get ray triangle intersection
    double intersect(const Eigen::Vector3d& o, const Eigen::Vector3d& d) const;
    
    // get point triangle distance
    double nearestPoint(const Eigen::Vector3d& p, Eigen::Vector3d& q) const;
        
    // checks whether face contains point
    bool containsPoint(const Eigen::Vector3d& p, double& dist) const;
};

#endif
