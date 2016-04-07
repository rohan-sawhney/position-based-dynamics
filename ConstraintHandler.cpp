#include "ConstraintHandler.h"
#include "Constraint.h"
#include "Mesh.h"

ConstraintHandler::ConstraintHandler():
count(0)
{
    
}

ConstraintHandler::~ConstraintHandler()
{
    clear();
}

void ConstraintHandler::generateConstraints(const double& kStretch, const double& kBend,
                                            Mesh& mesh)
{
    int rowVerts = (int)sqrt((double)mesh.vertices.size());
    for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        if (v->index < rowVerts || mesh.isRigid) v->invMass = 0.0;
    }
    
    if (!mesh.isRigid) {
        for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
            // add stretching constraint
            count++;
            constraints.push_back(new StretchingConstraint(e->he->vertex, e->he->next->vertex, kStretch));
            
            if (!e->isBoundary()) {
                // add bending constraint
                count++;
                constraints.push_back(new BendingConstraint(e->he->vertex,
                                                            e->he->next->vertex,
                                                            e->he->next->next->vertex,
                                                            e->he->flip->next->next->vertex,
                                                            kBend));
            }
        }
    }
}

void addStaticCollisionConstraints(Mesh& m1, Mesh& m2, std::vector<Constraint *>& constraints,
                                   const double& friction, const double& restitution)
{
    for (VertexIter v = m1.vertices.begin(); v != m1.vertices.end(); v++) {
        Eigen::Vector3d d = v->nPosition - v->position;
        Eigen::Vector3d q;
        double hit1 = d.norm(); d /= hit1;
        double hit2 = INFINITY;
        int index;
        if ((index = m2.bvh.getIntersection(RAY_INTERSECTION, hit1, q, v->position, d)) != -1 ||
            (index = m2.bvh.getIntersection(NEAREST_POINT, hit2, q, v->nPosition)) != -1) {
            
            Eigen::Vector3d n = m2.faces[index].normal().normalized();
            constraints.push_back(new StaticCollisionConstraint(v, q, n, friction, restitution));
        }
    }
}

void addTrianglePointCollisionConstraints(Mesh& m1, Mesh& m2, std::vector<Constraint *>& constraints,
                                          const double& friction, const double& restitution)
{
    for (FaceIter f = m1.faces.begin(); f != m1.faces.end(); f++) {
        Eigen::Vector3d n;
        double hit = INFINITY;
        int index;
        if ((index = m2.bvh.getIntersection(NEAREST_POINT_INV, hit, n, n, n, &(*f))) != -1) {
            VertexIter vj = m2.vertices.begin() + index;
            constraints.push_back(new TrianglePointCollisionConstraint(vj, f->he->vertex, f->he->next->vertex,
                                                                       f->he->next->next->vertex, n,
                                                                       friction, restitution));
        }
    }
}

void ConstraintHandler::generateCollisions(std::vector<Mesh>& meshes, const double& friction,
                                           const double& restitution)
{
    // clear any previous collisions
    clearCollisions();
    
    // add collision constraints
    for (size_t i = 0; i < meshes.size()-1; i++) {
        for (size_t j = i+1; j < meshes.size(); j++) {
            addStaticCollisionConstraints(meshes[i], meshes[j], constraints, friction, restitution);
            addTrianglePointCollisionConstraints(meshes[i], meshes[j], constraints, friction, restitution);
        }
    }
}

void ConstraintHandler::projectConstraints()
{
    // project constraints
    for (size_t i = 0; i < constraints.size(); i++) {
        constraints[i]->solve();
    }
}

void ConstraintHandler::updateVelocities()
{
    // update velocities due to collisions
    for (size_t i = count; i < constraints.size(); i++) {
        constraints[i]->updateVelocity();
    }
}

void ConstraintHandler::clearCollisions()
{
    // delete collision constraints
    for (size_t i = count; i < constraints.size(); i++) {
        delete constraints[i];
    }
    
    constraints.resize(count);
}

void ConstraintHandler::clear()
{
    count = 0;
    
    // delete constraints
    for (size_t i = count; i < constraints.size(); i++) {
        delete constraints[i];
    }
    
    constraints.clear();
}
