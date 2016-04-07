#include "Face.h"
#include "HalfEdge.h"
#include "Vertex.h"
#include "BoundingBox.h"

bool Face::isBoundary() const
{
    return he->onBoundary;
}

Eigen::Vector3d Face::normal() const
{
    const Eigen::Vector3d& a(he->vertex->position);
    const Eigen::Vector3d& b(he->next->vertex->position);
    const Eigen::Vector3d& c(he->next->next->vertex->position);
    
    return (a-b).cross(c-b);
}

Eigen::Vector3d Face::centroid() const
{
    const Eigen::Vector3d& a(he->vertex->position);
    const Eigen::Vector3d& b(he->next->vertex->position);
    const Eigen::Vector3d& c(he->next->next->vertex->position);
    
    return (a + b + c) / 3.0;
}

double Face::area() const
{
    if (isBoundary()) {
        return 0;
    }
    
    return 0.5 * normal().norm();
}

BoundingBox Face::boundingBox() const
{
    if (isBoundary()) {
        return BoundingBox(he->vertex->position, he->next->vertex->position);
    }
    
    const Eigen::Vector3d& a(he->vertex->position);
    const Eigen::Vector3d& b(he->next->vertex->position);
    const Eigen::Vector3d& c(he->next->next->vertex->position);
    
    Eigen::Vector3d min = a;
    Eigen::Vector3d max = a;
    
    if (b.x() < min.x()) min.x() = b.x();
    if (c.x() < min.x()) min.x() = c.x();
    
    if (b.y() < min.y()) min.y() = b.y();
    if (c.y() < min.y()) min.y() = c.y();
    
    if (b.z() < min.z()) min.z() = b.z();
    if (c.z() < min.z()) min.z() = c.z();
    
    if (b.x() > max.x()) max.x() = b.x();
    if (c.x() > max.x()) max.x() = c.x();
    
    if (b.y() > max.y()) max.y() = b.y();
    if (c.y() > max.y()) max.y() = c.y();
    
    if (b.z() > max.z()) max.z() = b.z();
    if (c.z() > max.z()) max.z() = c.z();
    
    return BoundingBox(min, max);
}

double Face::intersect(const Eigen::Vector3d& o, const Eigen::Vector3d& d) const
{
    // Möller–Trumbore intersection algorithm
    const Eigen::Vector3d& a(he->vertex->position);
    const Eigen::Vector3d& b(he->next->vertex->position);
    const Eigen::Vector3d& c(he->next->next->vertex->position);
    
    Eigen::Vector3d e1 = b - a;
    Eigen::Vector3d e2 = c - a;
    Eigen::Vector3d n = d.cross(e2);
    double det = e1.dot(n);
    
    // ray does not lie in the plane
    if (std::abs(det) < EPSILON) {
        return INFINITY;
    }
    
    double invDet = 1.0 / det;
    Eigen::Vector3d t = o - a;
    double u = t.dot(n) * invDet;
    
    // ray lies outside triangle
    if (u < 0.0 || u > 1.0) {
        return INFINITY;
    }
    
    Eigen::Vector3d q = t.cross(e1);
    double v = d.dot(q) * invDet;
    
    // ray lies outside the triangle
    if (v < 0.0 || v + u > 1.0) {
        return INFINITY;
    }
    
    // check for intersection
    double s = e2.dot(q) * invDet;
    if (s > 0) {
        return s;
    }
    
    // no hit
    return INFINITY;
}

double Face::nearestPoint(const Eigen::Vector3d& p, Eigen::Vector3d& q) const
{
    const Eigen::Vector3d& a(he->vertex->position);
    const Eigen::Vector3d& b(he->next->vertex->position);
    const Eigen::Vector3d& c(he->next->next->vertex->position);
    
    Eigen::Vector3d e1 = b - a;
    Eigen::Vector3d e2 = c - a;
    Eigen::Vector3d e3 = c - b;
    
    // check if p is outside vertex region a
    Eigen::Vector3d v1 = p - a;
    double d1 = e1.dot(v1), d2 = e2.dot(v1);
    if (d1 <= 0 && d2 <= 0) {
        q = a;
        return (p-q).norm();
    }
    
    // check if p is outside vertex region b
    Eigen::Vector3d v2 = p - b;
    double d3 = e1.dot(v2), d4 = e2.dot(v2);
    if (d3 >= 0 && d4 <= d3) {
        q = b;
        return (p-q).norm();
    }
    
    // check if p is in edge region e1, if so return projection of p onto e1
    double vc = d1*d4 - d3*d2;
    if (vc <= 0 && d1 >= 0 && d3 <= 0) {
        double v = d1 / (d1-d3);
        q = a + v * e1;
        return (p-q).norm();
    }
    
    // check if p in vertex region outside c
    Eigen::Vector3d v3 = p - c;
    double d5 = e1.dot(v3), d6 = e2.dot(v3);
    if (d6 >= 0 && d5 <= d6) {
        q = c;
        return (p-q).norm();
    }
    
    // check if p is in edge region e2, if so return projection of p onto e2
    double vb = d5*d2 - d1*d6;
    if (vb <= 0 && d2 >= 0 && d6 <= 0) {
        double w = d2 / (d2-d6);
        q = a + w * e2;
        return (p-q).norm();
    }
    
    // check if p is in edge region e3, if so return projection of p onto e3
    double va = d3*d6 - d5*d4;
    if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) {
        double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        q = b + w * e3;
        return (p-q).norm();
    }
    
    // p inside face region. Compute point through its barycentric coordinates (u,v,w)
    double d = 1.0 / (va + vb + vc);
    double v = vb * d;
    double w = vc * d;
    q = a + e1 * v + e2 * w;
    return (p-q).norm();
}

bool containsPoint(const Eigen::Vector3d& a, const Eigen::Vector3d& b,
                   const Eigen::Vector3d& c, const Eigen::Vector3d& p, double& dist)
{
    Eigen::Vector3d v1 = b - a;
    Eigen::Vector3d v2 = c - a;
    Eigen::Vector3d v3 = p - a;
    
    double dot11 = v1.dot(v1);
    double dot12 = v1.dot(v2);
    double dot13 = v1.dot(v3);
    double dot22 = v2.dot(v2);
    double dot23 = v2.dot(v3);
    
    // compute barycentric coordinates
    double invDen = 1 / (dot11*dot22 - dot12*dot12);
    double v = (dot22*dot13 - dot12*dot23) * invDen;
    double w = (dot11*dot23 - dot12*dot13) * invDen;
    double u = 1 - v - w;
    
    if ((u >= 0) && (v >= 0) && (u + v <= 1)) {
        dist = (p - (u*a + v*b + w*c)).norm();
        return true;
    }
    
    dist = INFINITY;
    return false;
}

bool Face::containsPointInPrism(const Eigen::Vector3d& p, double& dist) const
{
    bool hit1 = containsPoint(he->vertex->position,
                              he->next->vertex->position,
                              he->next->next->vertex->position,
                              p, dist);
    bool hit2 = containsPoint(he->vertex->nPosition,
                              he->next->vertex->nPosition,
                              he->next->next->vertex->nPosition,
                              p, dist);
    
    return hit1 && hit2;
}
