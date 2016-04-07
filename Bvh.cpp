#include "Bvh.h"
#include "Mesh.h"
#include <stack>

struct NodeEntry {
    // parent id
    int parentId;
    
    // range of objects covered by the node
    int startId, endId;
};

struct TraversalEntry {
    // constructor
    TraversalEntry(const int id0, const double d0) : id(id0), d(d0) {}
    
    // id
    int id;
    
    // distance
    double d;
};

Bvh::Bvh(const int& leafSize0):
leafSize(leafSize0)
{
    
}

void Bvh::build(Mesh *mesh0)
{
    mesh = mesh0;
    nodeCount = 0;
    leafCount = 0;
    flatTree.clear();
    
    int faceCount = (int)mesh->faces.size();
    
    NodeEntry nodeEntry;
    nodeEntry.parentId = -1;
    nodeEntry.startId = 0;
    nodeEntry.endId = faceCount;
    
    std::stack<NodeEntry> stack;
    stack.push(nodeEntry);
    
    Node node;
    flatTree.reserve(faceCount * 2);
    
    while (!stack.empty()) {
        // pop item off the stack and create a node
        nodeEntry = stack.top();
        stack.pop();
        int startId = nodeEntry.startId;
        int endId = nodeEntry.endId;
        
        nodeCount++;
        node.startId = startId;
        node.range = endId - startId;
        node.rightOffset = 2;
        
        // calculate bounding box
        BoundingBox boundingBox(mesh->faces[startId].boundingBox());
        BoundingBox boundingCentroid(mesh->faces[startId].centroid());
        for (int i = startId+1; i < endId; i++) {
            boundingBox.expandToInclude(mesh->faces[i].boundingBox());
            boundingCentroid.expandToInclude(mesh->faces[i].centroid());
        }
        node.boundingBox = boundingBox;
        
        // if node is a leaf
        if (node.range <= leafSize) {
            node.rightOffset = 0;
            leafCount++;
        }
        
        flatTree.push_back(node);
        
        // compute parent's rightOffset
        if (nodeEntry.parentId != -1) {
            flatTree[nodeEntry.parentId].rightOffset--;
            
            if (flatTree[nodeEntry.parentId].rightOffset == 0) {
                flatTree[nodeEntry.parentId].rightOffset = nodeCount - 1 - nodeEntry.parentId;
            }
        }
        
        // if a leaf, no need to subdivide
        if (node.rightOffset == 0) {
            continue;
        }
        
        // find the center of the longest dimension
        int maxDimension = boundingCentroid.maxDimension();
        double splitCoord = 0.5 * (boundingCentroid.min[maxDimension] +
                                   boundingCentroid.max[maxDimension]);

        
        // partition faces
        int mid = startId;
        for (int i = startId; i < endId; i++) {
            if (mesh->faces[i].centroid()[maxDimension] < splitCoord) {
                std::swap(mesh->faces[i], mesh->faces[mid]);
                mid ++;
            }
        }
        
        // in case of a bad split
        if (mid == startId || mid == endId) {
            mid = startId + (endId - startId) / 2;
        }
        
        // push right child
        nodeEntry.startId = mid;
        nodeEntry.endId = endId;
        nodeEntry.parentId = nodeCount - 1;
        stack.push(nodeEntry);
        
        // push left child
        nodeEntry.startId = startId;
        nodeEntry.endId = mid;
        nodeEntry.parentId = nodeCount - 1;
        stack.push(nodeEntry);
    }
}

int Bvh::getIntersection(const int& mode, double& hit, Eigen::Vector3d& q,
                         const Eigen::Vector3d& o, const Eigen::Vector3d& d,
                         const Face *f) const
{
    int id = 0;
    int closer, further;
    bool hit0 = false;
    bool hit1 = false;
    double dist1 = 0.0, dist2 = 0.0;
    int index = -1;
    
    if (hit < EPSILON) return -1;
    
    TraversalEntry t(id, -INFINITY);
    BoundingBox bbox;
    
    if (mode == RAY_INTERSECTION) flatTree[id].boundingBox.intersect(o, d, t.d);
    else if (mode == NEAREST_POINT) flatTree[id].boundingBox.intersect(o, t.d);
    else if (mode == NEAREST_POINT_INV) {
        bbox = f->boundingBox();
        flatTree[id].boundingBox.intersect(bbox, t.d);
    }
    
    std::stack<TraversalEntry> stack;
    stack.push(t);
    
    while (!stack.empty()) {
        TraversalEntry t = stack.top();
        id = t.id;
        stack.pop();
        
        if (hit <= t.d) continue;
        
        const Node &node(flatTree[id]);
        // node is a leaf
        if (node.rightOffset == 0) {
            for (int i = 0; i < node.range; i++) {
                const Face& face(mesh->faces[node.startId+i]);
                if (face.isBoundary()) continue;
                
                if (mode == RAY_INTERSECTION) {
                    double dist = face.intersect(o, d);
                    if (dist < hit) {
                        index = node.startId + i;
                        q = o + dist*d;
                        hit = dist;
                    }
                    
                } else if (mode == NEAREST_POINT) {
                    Eigen::Vector3d p;
                    double dist = face.nearestPoint(o, p);
                    if (dist < hit) {
                        index = node.startId + i;
                        q = p;
                        hit = dist;
                    }
                
                } else if (mode == NEAREST_POINT_INV) {
                    HalfEdgeIter h = face.he;
                    do {
                        if (index != h->vertex->index) {
                            double dist;
                            if (f->containsPointInPrism(h->vertex->position, dist)) {
                                if (dist < hit) {
                                    q = face.normal().normalized();
                                    index = h->vertex->index;
                                    hit = dist;
                                }
                            }
                        }
                        
                        h = h->next;
                    } while (h != face.he);
                }
            }
            
        } else { // not a leaf
            // check for intersection with bounding box
            if (mode == RAY_INTERSECTION) {
                hit0 = flatTree[id+1].boundingBox.intersect(o, d, dist1);
                hit1 = flatTree[id+node.rightOffset].boundingBox.intersect(o, d, dist2);
                
            } else if (mode == NEAREST_POINT) {
                hit0 = flatTree[id+1].boundingBox.intersect(o, dist1);
                hit1 = flatTree[id+node.rightOffset].boundingBox.intersect(o, dist2);
            
            } else if (mode == NEAREST_POINT_INV) {
                hit0 = flatTree[id+1].boundingBox.intersect(bbox, dist1);
                hit1 = flatTree[id+node.rightOffset].boundingBox.intersect(bbox, dist2);
            }
            
            // hit both bounding boxes
            if (hit0 && hit1) {
                closer = id + 1;
                further = id + node.rightOffset;
                
                if (dist2 < dist1) {
                    std::swap(dist1, dist2);
                    std::swap(closer, further);
                }
                
                // push farther node first
                stack.push(TraversalEntry(further, dist2));
                stack.push(TraversalEntry(closer, dist1));
                
            } else if (hit0) {
                stack.push(TraversalEntry(id+1, dist1));
                
            } else if (hit1) {
                stack.push(TraversalEntry(id+node.rightOffset, dist2));
            }
        }
    }
    
    return index;
}
