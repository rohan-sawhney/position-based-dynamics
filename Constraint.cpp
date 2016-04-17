#include "Constraint.h"
#include "Mesh.h"

Constraint::Constraint(const int& n0, const double& k0):
n(n0),
k(k0)
{
    
}

Constraint::~Constraint()
{
    
}

StretchingConstraint::StretchingConstraint(VertexIter v1, VertexIter v2,
                                           const double& k0):
Constraint(2, k0)
{
    vs = { v1, v2 };
    l = (v1->position-v2->position).norm();
}

StretchingConstraint::~StretchingConstraint()
{
    
}

void StretchingConstraint::solve()
{
    VertexIter& v1(vs[0]);
    VertexIter& v2(vs[1]);
    
    // compute and add correction
    Eigen::Vector3d q = v1->nPosition - v2->nPosition;
    double d = q.norm();
    q.normalize();
    
    double sum = v1->invMass + v2->invMass;
    if (sum < EPSILON) return;
        
    double lambda = k*(d-l) / sum;
    v1->nPosition -= v1->invMass*lambda*q;
    v2->nPosition += v2->invMass*lambda*q;
}

void StretchingConstraint::updateVelocity()
{
    // nothing to do
}

BendingConstraint::BendingConstraint(VertexIter v1, VertexIter v2,
                                     VertexIter v3, VertexIter v4,
                                     const double& k0):
Constraint(4, k0)
{
    vs = { v1, v2, v3, v4 };
    
    // compute dihedral angle
    Eigen::Vector3d n1 = (v1->position-v3->position).cross(v2->position-v3->position).normalized();
    Eigen::Vector3d n2 = (v2->position-v4->position).cross(v1->position-v4->position).normalized();
    double d = n1.dot(n2);
    if (d < -1.0) d = -1.0;
    else if (d >  1.0) d = 1.0;
    theta = acos(d);
}

BendingConstraint::~BendingConstraint()
{
    
}

void BendingConstraint::solve()
{
    VertexIter& v1(vs[0]);
    VertexIter& v2(vs[1]);
    VertexIter& v3(vs[2]);
    VertexIter& v4(vs[3]);
    
    // compute and add correction
    Eigen::Vector3d e = v2->nPosition-v1->nPosition;
    Eigen::Vector3d n1 = (v1->nPosition-v3->nPosition).cross(v2->nPosition-v3->nPosition);
    Eigen::Vector3d n2 = (v2->nPosition-v4->nPosition).cross(v1->nPosition-v4->nPosition);
    
    double de = e.norm();
    double dn1 = n1.norm(); n1 /= dn1;
    double dn2 = n2.norm(); n2 /= dn2;
    
    Eigen::Vector3d q1 = ((v3->nPosition - v2->nPosition).dot(e)*n1 / dn1 +
                          (v4->nPosition - v2->nPosition).dot(e)*n2 / dn2) / de;
    Eigen::Vector3d q2 = ((v1->nPosition - v3->nPosition).dot(e)*n1 / dn1 +
                          (v1->nPosition - v4->nPosition).dot(e)*n2 / dn2) / de;
    Eigen::Vector3d q3 = de * n1 / dn1;
    Eigen::Vector3d q4 = de * n2 / dn2;
    
    double sum = v1->invMass*q1.squaredNorm() +
                 v2->invMass*q2.squaredNorm() +
                 v3->invMass*q3.squaredNorm() +
                 v4->invMass*q4.squaredNorm();
    if (sum < EPSILON) return;
    
    double d = n1.dot(n2);
    if (d < -1.0) d = -1.0;
    else if (d >  1.0) d = 1.0;
    
    double lambda = k*(acos(d)-theta) / sum;
    if (n1.cross(n2).dot(e) > 0.0) lambda = -lambda;
    
    v1->nPosition -= v1->invMass*lambda*q1;
    v2->nPosition -= v2->invMass*lambda*q2;
    v3->nPosition -= v3->invMass*lambda*q3;
    v4->nPosition -= v4->invMass*lambda*q4;
}

void BendingConstraint::updateVelocity()
{
    // nothing to do
}

StaticCollisionConstraint::StaticCollisionConstraint(VertexIter v1, const Eigen::Vector3d& q0,
                                                     const Eigen::Vector3d& normal0,
                                                     const double& friction0, const double& restitution0):
Constraint(1, 1.0),
q(q0),
normal(normal0),
friction(friction0),
restitution(restitution0),
didSolve(false)
{
    vs = { v1 };
}

StaticCollisionConstraint::~StaticCollisionConstraint()
{
    
}

void StaticCollisionConstraint::solve()
{
    VertexIter& v1(vs[0]);
    
    // compute and add correction
    double lambda = normal.dot(v1->nPosition-q);
    if (lambda < 0) {
        v1->nPosition -= lambda*normal;
        didSolve = true;
    }
}

void StaticCollisionConstraint::updateVelocity()
{
    if (didSolve) {
        VertexIter& v1(vs[0]);
        Eigen::Vector3d vn = v1->velocity.dot(normal)*normal;
        Eigen::Vector3d vt = v1->velocity - vn;
        v1->velocity = friction*vt + restitution*vn;
    }
}

TrianglePointCollisionConstraint::TrianglePointCollisionConstraint(VertexIter v1, VertexIter v2,
                                                                   VertexIter v3, VertexIter v4,
                                                                   const double& friction0,
                                                                   const double& restitution0):
Constraint(4, 1.0),
friction(friction0),
restitution(restitution0),
didSolve(false)
{
    vs = { v1, v2, v3, v4 };
    normal = (v3->nPosition-v2->nPosition).cross(v4->nPosition-v2->nPosition);
}

TrianglePointCollisionConstraint::~TrianglePointCollisionConstraint()
{
    
}

void TrianglePointCollisionConstraint::solve()
{
    /*
    VertexIter& v1(vs[0]);
    VertexIter& v2(vs[1]);
    VertexIter& v3(vs[2]);
    VertexIter& v4(vs[3]);
    
    // compute and add correction
    Eigen::Vector3d u = v1->nPosition-v2->nPosition;
    Eigen::Vector3d v = v3->nPosition-v2->nPosition;
    Eigen::Vector3d w = v4->nPosition-v2->nPosition;
    double dn = normal.norm();
    normal /= dn;
    
    double lambda = u.dot(normal);
    if (lambda < 0) {
        Eigen::Vector3d q1 = normal;
        Eigen::Vector3d q3 = (w.cross(u) + normal.cross(w)*lambda) / dn;
        Eigen::Vector3d q4 = -(v.cross(u) + normal.cross(v)*lambda) / dn;
        Eigen::Vector3d q2 = -(q1 + q3 + q4);
        
        double sum = v1->invMass*q1.squaredNorm() +
                     v2->invMass*q2.squaredNorm() +
                     v3->invMass*q3.squaredNorm() +
                     v4->invMass*q4.squaredNorm();
        if (sum < EPSILON) return;
        lambda /= sum;
        
        v2->nPosition -= v2->invMass*lambda*q2;
        v3->nPosition -= v3->invMass*lambda*q3;
        v4->nPosition -= v4->invMass*lambda*q4;
        
        didSolve = true;
    }
    */
}

void TrianglePointCollisionConstraint::updateVelocity()
{
    /*
    if (didSolve) {
        // incorrect
        VertexIter& v2(vs[1]);
        VertexIter& v3(vs[2]);
        VertexIter& v4(vs[3]);
        
        Eigen::Vector3d vn = v2->velocity.dot(normal)*normal;
        Eigen::Vector3d vt = v2->velocity - vn;
        v2->velocity = friction*vt + restitution*vn;
        
        vn = v3->velocity.dot(normal)*normal;
        vt = v3->velocity - vn;
        v3->velocity = friction*vt + restitution*vn;

        vn = v4->velocity.dot(normal)*normal;
        vt = v4->velocity - vn;
        v4->velocity = friction*vt + restitution*vn;
    }
    */
}

