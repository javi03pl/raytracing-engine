#include "texture.h"
#include <cmath>
#include <glm/geometric.hpp>
#include <shading.h>

const glm::vec3 computeShading(const glm::vec3& lightPosition, const glm::vec3& lightColor, const Features& features, Ray ray, HitInfo hitInfo)
{
    // If shading not enabled, return kd 
    if(!features.enableShading) return hitInfo.material.kd;

    //Get variables
    glm::vec3 normal = glm::normalize(hitInfo.normal);
    glm::vec3 intersectionPos = ray.origin + ray.t * ray.direction;
    glm::vec3 viewDir = glm::normalize(ray.origin - intersectionPos);
    glm::vec3 lightDir = glm::normalize(lightPosition - intersectionPos);
    glm::vec3 reflectionDir = glm::normalize(2 * glm::dot(lightDir, normal) * normal - lightDir);

    //Diffuse Component (Formula in lectures
    glm::vec3 diffuse = hitInfo.material.kd * glm::max(glm::dot(lightDir, normal), 0.0f) * lightColor;
    
    //Specular Component (Formula in lectures)
    glm::vec3 specular = glm::vec3(0.0f, 0.0f, 0.0f);
    if (glm::dot(lightDir, normal) > 0) {
        specular = hitInfo.material.ks * glm::pow(glm::max(glm::dot(reflectionDir, viewDir), 0.0f), hitInfo.material.shininess) * lightColor;
    }
    //Compute sum of both components
    return diffuse + specular;
}


const Ray computeReflectionRay (Ray ray, HitInfo hitInfo)
{
    // Do NOT use glm::reflect!! write your own code.
    // TODO: implement the reflection ray computation.

    glm::vec3 normal = hitInfo.normal;  // the normal at the intersection point
    glm::vec3 newOrigin = ray.origin + ray.direction * ray.t;   // intersection point - will be the origin of the reflection
    glm::vec3 incomingVector = ray.direction;
    // compute the direction of the reflected ray (reflect the direction vector of the old ray against the normal vector))
    glm::vec3 reflectedVector = incomingVector - 2 * glm::dot(incomingVector, normal) * normal;

    Ray reflectionRay {newOrigin, reflectedVector};
    return reflectionRay;
}