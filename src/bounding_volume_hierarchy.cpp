#include "bounding_volume_hierarchy.h"
#include "draw.h"
#include "interpolate.h"
#include "intersect.h"
#include "scene.h"
#include "texture.h"
#include <glm/glm.hpp>
#include <iostream>



Triangle get_triangle(int index, Scene& sc)
{
    for (Mesh& mesh : sc.meshes) {
        if (mesh.triangles.size() <= index) {
            //index = index - mesh.triangles.size();
            return Triangle();
        } else {
            glm::uvec3 triangle = mesh.triangles.at(index);
            return Triangle(mesh.vertices.at(triangle.x), mesh.vertices.at(triangle.y), mesh.vertices.at(triangle.z));
        }
    }
    
}

AxisAlignedBox computeAABB(std::vector<int>& triangle_indices, Scene& scene)
{
    AxisAlignedBox res;
    float minX = get_triangle(triangle_indices[0], scene).v0.position.x;
    float minY = get_triangle(triangle_indices[0], scene).v0.position.y;
    float minZ = get_triangle(triangle_indices[0], scene).v0.position.z;
    float maxX = get_triangle(triangle_indices[0], scene).v0.position.x;
    float maxY = get_triangle(triangle_indices[0], scene).v0.position.y;
    float maxZ = get_triangle(triangle_indices[0], scene).v0.position.z;

    for (int t : triangle_indices) {
        if (get_triangle(t, scene).v0.position.x > maxX)
            maxX = get_triangle(t, scene).v0.position.x;
        if (get_triangle(t, scene).v1.position.x > maxX)
            maxX = get_triangle(t, scene).v1.position.x;
        if (get_triangle(t, scene).v2.position.x > maxX)
            maxX = get_triangle(t, scene).v2.position.x;
    }

    for (int t : triangle_indices) {
        if (get_triangle(t, scene).v0.position.x < minX)
            minX = get_triangle(t, scene).v0.position.x;
        if (get_triangle(t, scene).v1.position.x < minX)
            minX = get_triangle(t, scene).v1.position.x;
        if (get_triangle(t, scene).v2.position.x < minX)
            minX = get_triangle(t, scene).v2.position.x;
    }

    for (int t : triangle_indices) {
        if (get_triangle(t, scene).v0.position.y > maxY)
            maxY = get_triangle(t, scene).v0.position.y;
        if (get_triangle(t, scene).v1.position.y > maxY)
            maxY = get_triangle(t, scene).v1.position.y;
        if (get_triangle(t, scene).v2.position.y > maxY)
            maxY = get_triangle(t, scene).v2.position.y;
    }

    for (int t : triangle_indices) {
        if (get_triangle(t, scene).v0.position.y < minY)
            minY = get_triangle(t, scene).v0.position.y;
        if (get_triangle(t, scene).v1.position.y < minY)
            minY = get_triangle(t, scene).v1.position.y;
        if (get_triangle(t, scene).v2.position.y < minY)
            minY = get_triangle(t, scene).v2.position.y;
    }
    for (int t : triangle_indices) {
        if (get_triangle(t, scene).v0.position.z > maxZ)
            maxZ = get_triangle(t, scene).v0.position.z;
        if (get_triangle(t, scene).v1.position.z > maxZ)
            maxZ = get_triangle(t, scene).v1.position.z;
        if (get_triangle(t, scene).v2.position.z > maxZ)
            maxZ = get_triangle(t, scene).v2.position.z;
    }
    for (int t : triangle_indices) {
        if (get_triangle(t, scene).v0.position.z < minZ)
            minZ = get_triangle(t, scene).v0.position.z;
        if (get_triangle(t, scene).v1.position.z < minZ)
            minZ = get_triangle(t, scene).v1.position.z;
        if (get_triangle(t, scene).v2.position.z < minZ)
            minZ = get_triangle(t, scene).v2.position.z;
    }
    res.lower = glm::vec3(minX, minY, minZ);
    res.upper = glm::vec3(maxX, maxY, maxZ);
    return res;
}

void sortAccordingToAxis(std::vector<int>& triangle_indices, int axis, Scene& scene)
{
    switch (axis % 3) {
    case 0:
        std::sort(triangle_indices.begin(), triangle_indices.end(), [&](int lhs, int rhs) {
            float minLeft = std::min(std::min(get_triangle(lhs, scene).v0.position.x, get_triangle(lhs, scene).v1.position.x), get_triangle(lhs, scene).v2.position.x);
            float minRight = std::min(std::min(get_triangle(rhs, scene).v0.position.x, get_triangle(rhs, scene).v1.position.x), get_triangle(rhs, scene).v2.position.x);
            return minLeft < minRight;
        });
        break;
    case 1:
        std::sort(triangle_indices.begin(), triangle_indices.end(), [&](int lhs, int rhs) {
            float minLeft = std::min(std::min(get_triangle(lhs, scene).v0.position.y, get_triangle(lhs, scene).v1.position.y), get_triangle(lhs, scene).v2.position.y);
            float minRight = std::min(std::min(get_triangle(rhs, scene).v0.position.y, get_triangle(rhs, scene).v1.position.y), get_triangle(rhs, scene).v2.position.y);
            return minLeft < minRight;
        });
        break;
    case 2:
        std::sort(triangle_indices.begin(), triangle_indices.end(), [&](int lhs, int rhs) {
            float minLeft = std::min(std::min(get_triangle(lhs, scene).v0.position.z, get_triangle(lhs, scene).v1.position.z), get_triangle(lhs, scene).v2.position.z);
            float minRight = std::min(std::min(get_triangle(rhs, scene).v0.position.z, get_triangle(rhs, scene).v1.position.z), get_triangle(rhs, scene).v2.position.z);
            return minLeft < minRight;
        });
        break;
    }
}

Triangle::Triangle(Vertex& v0, Vertex& v1, Vertex& v2)
{
    this->v0 = v0;
    this->v1 = v1;
    this->v2 = v2;
    
}

Triangle::Triangle() {

}

Node::Node(std::vector<int>& triangle_indices, int level, int axis, std::vector<Node>& nodes, Scene& scene)
{
    this->axis = axis;
    this->level = level;
    this->bounds = computeAABB(triangle_indices, scene);
    if (triangle_indices.size() > 10) {
        sortAccordingToAxis(triangle_indices, axis, scene);
        int i = triangle_indices.size() / 2;
        //printf("%d\n", i);
        std::vector<int> left_triangles = std::vector<int>(triangle_indices.begin(), triangle_indices.begin() + i);
        std::vector<int> right_triangles = std::vector<int>(triangle_indices.begin() + i + 1, triangle_indices.end());

        Node leftNode = Node(left_triangles, level + 1, axis + 1, nodes, scene);
        Node rightNode = Node(right_triangles, level + 1, axis + 1, nodes, scene);

        leftChild = nodes.size();
        nodes.push_back(leftNode);
        rightChild = nodes.size();
        nodes.push_back(rightNode);
        this->isLeaf = false;
    } else {
        this->level = level;
        this->triangle_indices = triangle_indices;
        this->isLeaf = true;
    }
}

BoundingVolumeHierarchy::BoundingVolumeHierarchy(Scene* pScene)
    : m_pScene(pScene)
{
    std::vector<int> triangle_indices;
    int counter = 0;
    for (Mesh& mesh : this->m_pScene->meshes) {
        for (const auto& tri : mesh.triangles) {
            //triangles.push_back(Triangle(mesh.vertices[tri[0]], mesh.vertices[tri[1]], mesh.vertices[tri[2]]));
            triangle_indices.push_back(counter);
            counter++;
        }
    }
    int level = 0;
    int axis = 0; // x
    Node root = Node(triangle_indices, level, axis, this->nodes, *(m_pScene));
    nodes.insert(nodes.begin(), root);
    //int kur = 0; // debug, delete later
    //printf("%d ------", nodes.size());
}

// Return the depth of the tree that you constructed. This is used to tell the
// slider in the UI how many steps it should display for Visual Debug 1.
int BoundingVolumeHierarchy::numLevels() const
{
    int i = 0;
    for (const Node& n : nodes) {
        i = (n.level > i ? n.level : i);
    }
    return i + 1;
}

// Return the number of leaf nodes in the tree that you constructed. This is used to tell the
// slider in the UI how many steps it should display for Visual Debug 2.
int BoundingVolumeHierarchy::numLeaves() const
{
    int i = 0;
    for (const Node& n : nodes) {
        i += (n.isLeaf ? 1 : 0);
    }
    return i;
}

// Use this function to visualize your BVH. This is useful for debugging. Use the functions in
// draw.h to draw the various shapes. We have extended the AABB draw functions to support wireframe
// mode, arbitrary colors and transparency.
void BoundingVolumeHierarchy::debugDrawLevel(int level)
{
    for (const Node& n : nodes) {
        if (n.level == level) {
            drawAABB(n.bounds, DrawMode::Wireframe);
        }
    }
}

// Use this function to visualize your leaf nodes. This is useful for debugging. The function
// receives the leaf node to be draw (think of the ith leaf node). Draw the AABB of the leaf node and all contained triangles.
// You can draw the triangles with different colors. NoteL leafIdx is not the index in the node vector, it is the
// i-th leaf node in the vector.
void BoundingVolumeHierarchy::debugDrawLeaf(int leafIdx)
{
    int i = 0;
    for (const Node& n : nodes) {
        if (n.isLeaf) {
            i++;
            if (i == leafIdx) {
                drawAABB(n.bounds, DrawMode::Wireframe);
                // break;
            }
        }
    }
}

// Return true if something is hit, returns false otherwise. Only find hits if they are closer than t stored
// in the ray and if the intersection is on the correct side of the origin (the new t >= 0). Replace the code
// by a bounding volume hierarchy acceleration structure as described in the assignment. You can change any
// file you like, including bounding_volume_hierarchy.h.
bool BoundingVolumeHierarchy::intersect(Ray& ray, HitInfo& hitInfo, const Features& features) const
{
    // If BVH is not enabled, use the naive implementation.
    if (!features.enableAccelStructure) {
        bool hit = false;
        Vertex triangleV0, triangleV1, triangleV2;
        // Intersect with all triangles of all meshes.
        for (const auto& mesh : m_pScene->meshes) {
            for (const auto& tri : mesh.triangles) {
                const auto v0 = mesh.vertices[tri[0]];
                const auto v1 = mesh.vertices[tri[1]];
                const auto v2 = mesh.vertices[tri[2]];
                //std::cout << intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo) << std::endl;
                if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {
                    hitInfo.material = mesh.material;

                    //Used for shading -> Set hitInfo normal to cross product of vectors between vertices
                    glm::vec3 vector1 = glm::normalize(v0.position - v1.position);
                    glm::vec3 vector2 = glm::normalize(v0.position - v2.position);
                    hitInfo.normal = glm::normalize(glm::cross(vector1, vector2));

                    triangleV0 = v0;
                    triangleV1 = v1;
                    triangleV2 = v2;

                    //Used for texture mapping -> Calculate barycentric coordinates of intersection point
                    glm::vec3 intersectionPosition = ray.origin + ray.direction * ray.t;
                    hitInfo.barycentricCoord = computeBarycentricCoord(v0.position, v1.position, v2.position, intersectionPosition);

                    if (features.enableTextureMapping) {
                        Image img = *mesh.material.kdTexture.get();          //Get texture image defined in material
                        //Get texture coordinates at intersectionPosition
                        glm::vec2 texCoords = interpolateTexCoord(v0.texCoord, v1.texCoord, v2.texCoord, hitInfo.barycentricCoord);
                        //Get kd at texture coordinates & set kd of hitInfo
                        hitInfo.material.kd = acquireTexel(img, texCoords, features);
                    }

                    hit = true;
                }
            }
        }
        // Intersect with spheres.
        for (const auto& sphere : m_pScene->spheres)
            hit |= intersectRayWithShape(sphere, ray, hitInfo);
        interpolationVisualDebug(triangleV0, triangleV1, triangleV2, ray, features);
        return hit;
    } else {

        // TODO: implement here the bounding volume hierarchy traversal.
        // Please note that you should use `features.enableNormalInterp` and `features.enableTextureMapping`
        // to isolate the code that is only needed for the normal interpolation and texture mapping features.
        return false;
    }
}
