#include "ClothGenerator.h"
#include "Mesh.h"
#include "MeshIO.h"

ClothGenerator::ClothGenerator()
{
    
}

bool ClothGenerator::generate(const double& extent, const double& base,
                              const int& intervals, Mesh& mesh)
{
    MeshData data;
    const double inc = extent / intervals;
    const double start = -extent / 2.0;
    const double end = extent / 2.0;
    
    // add vertices
    double shift = 0.0;
    for (double y = end; y >= start; y -= inc) {
        for (double x = start; x <= end; x += inc) {
            data.positions.push_back(Eigen::Vector3d(x, y / 2.0 + shift, base + shift));
        }
        shift += 0.05;
    }
    
    // add indices
    int col = intervals+1;
    for (int i = 0; i < intervals; i++) {
        for (int j = 0; j < intervals; j++) {
            
            int x = 0;
            if (i%2 == j%2) x = 1;
            std::vector<Index> indices = {Index(i*col+j      , -1, -1),
                                          Index(i*col+j+1    , -1, -1),
                                          Index((i+1)*col+j+x, -1, -1)};
            data.indices.push_back(indices);
            
            indices[0].position = (i+1)*col+j+1;
            indices[1].position = (i+1)*col+j;
            indices[2].position = i*col+j+1-x;
            data.indices.push_back(indices);
        }
    }
    
    return MeshIO::buildMesh(data, mesh);
}

