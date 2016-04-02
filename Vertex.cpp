#include "Vertex.h"
#include "Edge.h"
#include "Face.h"
#include "HalfEdge.h"

std::vector<HalfEdge> isolated;

Vertex::Vertex():
position(Eigen::Vector3d::Zero()),
nPosition(Eigen::Vector3d::Zero()),
velocity(Eigen::Vector3d::Zero()),
mass(1.0),
invMass(1.0 / mass),
index(-1)
{

}

bool Vertex::isIsolated() const
{
    return he == isolated.begin();
}

Eigen::Vector3d Vertex::normal() const
{
    Eigen::Vector3d normal;
    Eigen::Vector3d n;
    normal.setZero();
    
    if (isIsolated()) {
        return normal;
    }
    
    double angle;
    HalfEdgeCIter h = he;
    do {
        Eigen::Vector3d e1 = h->flip->vertex->position - position;
        Eigen::Vector3d e2 = h->flip->next->flip->vertex->position - position;
        
        double c = e1.dot(e2) / sqrt(e1.dot(e1) * e2.dot(e2));
        if (c < -1.0) c = -1.0;
        else if (c >  1.0) c = 1.0;
        angle = acos(c);
        
        n = h->face->normal();
        if (n.squaredNorm() == 0.0) {
            n.setZero();
        }
        
        normal += angle * n;
        
        h = h->flip->next;
    } while (h != he);
    
    normal.normalize();
    return normal;
}
