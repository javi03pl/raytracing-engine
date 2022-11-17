#include "interpolate.h"
#include <glm/geometric.hpp>

glm::vec3 computeBarycentricCoord (const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& p)
{
    // TODO: implement this function.
    float triangleArea = glm::length(glm::cross(v0 - v2, v1 - v2));
    float alpha = glm::length(glm::cross(v1 - p, v2 - p)) / triangleArea;
    float beta = glm::length(glm::cross(v0 - p, v2 - p)) / triangleArea;
    float gamma = glm::length(glm::cross(v0 - p, v1 - p)) / triangleArea;
    return glm::vec3(alpha, beta, gamma);
}

glm::vec3 interpolateNormal (const glm::vec3& n0, const glm::vec3& n1, const glm::vec3& n2, const glm::vec3 barycentricCoord)
{
    return glm::normalize(n0 * barycentricCoord.x + n1 * barycentricCoord.y + n2 * barycentricCoord.z);
}

//Returns interolated texture coordinate using weights inside barycentricCoord
glm::vec2 interpolateTexCoord (const glm::vec2& t0, const glm::vec2& t1, const glm::vec2& t2, const glm::vec3 barycentricCoord)
{
    return glm::vec2 {
        t0.x * barycentricCoord.x + t1.x * barycentricCoord.y + t2.x * barycentricCoord.z,
        t0.y * barycentricCoord.x + t1.y * barycentricCoord.y + t2.y * barycentricCoord.z,
    };
}
