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
    Eigen::Vector3d normal = Eigen::Vector3d::Zero();
    if (isIsolated()) return normal;
    
    HalfEdgeCIter h = he;
    do {
        Eigen::Vector3d e1 = h->next->vertex->position - position;
        Eigen::Vector3d e2 = h->next->next->vertex->position - position;
        
        double d = e1.dot(e2) / sqrt(e1.squaredNorm() * e2.squaredNorm());
        if (d < -1.0) d = -1.0;
        else if (d >  1.0) d = 1.0;
        double angle = acos(d);
        
        Eigen::Vector3d n = h->face->normal();
        normal += angle * n;
        
        h = h->flip->next;
    } while (h != he);
    
    if (!normal.isZero()) normal.normalize();
    return normal;
}
