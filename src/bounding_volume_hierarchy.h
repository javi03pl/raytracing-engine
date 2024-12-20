#pragma once
#include "common.h"
#include <array>
#include <framework/ray.h>
#include <vector>

// Forward declaration.
struct Scene;
struct Triangle {
    Vertex v0;
    Vertex v1;
    Vertex v2;

    Triangle();
    Triangle(Vertex& v0, Vertex& v1, Vertex& v2);
};


struct Node {
    bool isLeaf;
    AxisAlignedBox bounds;
    int leftChild;
    int rightChild;
    //std::vector<Triangle> triangles;
    std::vector<int> triangle_indices;
    int level;
    int axis; // 0 for x, 1 for y, 2 for z when using %

    Node(std::vector<int>& trianangle_indices, int level, int axis, std::vector<Node>& nodes, Scene& scene);
};

class BoundingVolumeHierarchy {
public:
    std::vector<Node> nodes;

    // Constructor. Receives the scene and builds the bounding volume hierarchy.
    BoundingVolumeHierarchy(Scene* pScene);

    // Return how many levels there are in the tree that you have constructed.
    [[nodiscard]] int numLevels() const;

    // Return how many leaf nodes there are in the tree that you have constructed.
    [[nodiscard]] int numLeaves() const;

    // Visual Debug 1: Draw the bounding boxes of the nodes at the selected level.
    void debugDrawLevel(int level);

    // Visual Debug 2: Draw the triangles of the i-th leaf
    void debugDrawLeaf(int leafIdx);

    // Return true if something is hit, returns false otherwise.
    // Only find hits if they are closer than t stored in the ray and the intersection
    // is on the correct side of the origin (the new t >= 0).
    bool intersect(Ray& ray, HitInfo& hitInfo, const Features& features) const;


private:
    int m_numLevels;
    int m_numLeaves;
    Scene* m_pScene;
};




