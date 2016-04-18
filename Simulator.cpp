#include "Simulator.h"

#define TIME_STEP 0.01
#define K_STRETCH 1.0
#define K_BEND 1.0
#define FRICTION 1.0
#define RESTITUTION 1.0
#define MAX_ITERS 5
const Eigen::Vector3d gravity(0.0, -9.8, 0.0);

Simulator::Simulator(std::vector<Mesh>& meshes0):
meshes(meshes0)
{
    
}

void Simulator::initialize()
{
    // initialize constraints
    double p = 1.0 / MAX_ITERS;
    const double ks = 1 - pow(1-K_STRETCH, p);
    const double kb = 1 - pow(1-K_BEND, p);
    for (size_t i = 0; i < meshes.size(); i++) {
        constraintHandler.generateConstraints(ks, kb, meshes[i]);
    }
}

void stepExplicitEuler(VertexIter v)
{
    v->velocity += gravity * v->invMass * TIME_STEP;
    v->nPosition = v->position + v->velocity * TIME_STEP;
}

void updatePositionAndVelocity(VertexIter v)
{
    v->velocity = (v->nPosition - v->position) / TIME_STEP;
    v->position = v->nPosition;
}

void Simulator::step()
{
    // TODO:
    // 1) Fix Triangle Point collision constraint
    // 2) Damp Velocities
    
    // step explicit euler
    for (size_t i = 0; i < meshes.size(); i++) {
        for (VertexIter v = meshes[i].vertices.begin(); v != meshes[i].vertices.end(); v++) {
            stepExplicitEuler(v);
        }
    }
    
    // TODO: damp velocity
    
    // generate collisions
    constraintHandler.generateCollisions(meshes, FRICTION, RESTITUTION);
    
    // project constraints
    for (int i = 0; i < MAX_ITERS; i++) {
        constraintHandler.projectConstraints();
    }
    
    // set positions and velocities
    for (size_t i = 0; i < meshes.size(); i++) {
        for (VertexIter v = meshes[i].vertices.begin(); v != meshes[i].vertices.end(); v++) {
            updatePositionAndVelocity(v);
        }
    }
    
    // update velocities from collisions
    constraintHandler.updateVelocities();
}

void Simulator::reset()
{
    constraintHandler.clear();
}