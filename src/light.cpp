#include "light.h"
#include "config.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
DISABLE_WARNINGS_POP()
#include <cmath>
#include <draw.h>

//Calculates a random point and its interpolated color inside a partititon of a segment
void sampleSegmentLight(const SegmentLight& segmentLight, glm::vec3& position, glm::vec3& color)
{
    //Get weights to calculate interpolated position -> Random points are calculated by casting a ray and choosing a random t (inside a min & max)
    float minWeight = position.x;
    float maxWeight = position.y;
    float randomWeight = minWeight + static_cast<float>(rand()) / (static_cast<float>(RAND_MAX / (maxWeight - minWeight)));   //Get a random weight between min and max
    float w1 = 1.0f - randomWeight;
    float w2 = randomWeight;
    position = w1 * segmentLight.endpoint0 + w2 * segmentLight.endpoint1;           //Multiply ednpoints by weights to get inerpolated position
    //printf("SECOND: minW = %f, maxW = %f, randW = %f ---------", minWeight, maxWeight, randomWeight);

    // Get interpolated color at position by using distances
    float vecLength = glm::length(segmentLight.endpoint0 - segmentLight.endpoint1);
    float dist1 = glm::length(segmentLight.endpoint0 - position);
    float dist2 = glm::length(segmentLight.endpoint1 - position);
    color = w1 * segmentLight.color0 + w2 * segmentLight.color1;
    // printf("Color: w(endpoint0) = %f, w(endpoint1) = %f", w1, w2);
}

// Calculates a random point inside a gid in a parallelogram & gets its interpolated color
void sampleParallelogramLight(const ParallelogramLight& parallelogramLight, glm::vec3& position, glm::vec3& color)
{
    // Get weights to calculate interpolated position -> Random points are calculated by casting a ray and choosing a random t (inside a min & max)
    float minWeight01 = position.x;
    float maxWeight01 = position.y;
    float minWeight02 = color.x;
    float maxWeight02 = color.y;
    float randomWeight01 = minWeight01 + static_cast<float>(rand()) / (static_cast<float>(RAND_MAX / (maxWeight01 - minWeight01))); // Get a random weight between min and max (Random x inside grid)
    float randomWeight02 = minWeight02 + static_cast<float>(rand()) / (static_cast<float>(RAND_MAX / (maxWeight02 - minWeight02))); // Get a random weight between min and max (Random y inside grid)
    position = parallelogramLight.v0 + parallelogramLight.edge01 * randomWeight01 + parallelogramLight.edge02 * randomWeight02;
    //printf("SECOND: minW01 = %f, minW02 = %f, randW01 = %f, randW02: %f ---------", minWeight01, minWeight02, randomWeight01, randomWeight02);

    //Calculate interpolated color using bilinear interpolation formula from slides
    float alpha = randomWeight02;
    float beta = randomWeight01;
    color = ((1.0f - alpha) * parallelogramLight.color0 + alpha * parallelogramLight.color2) * (1.0f - beta) + ((1.0f - alpha) * parallelogramLight.color1 + alpha * parallelogramLight.color3) * beta;
}

// test the visibility at a given light sample
// returns 1.0 if sample is visible, 0.0 otherwise
float testVisibilityLightSample(const glm::vec3& samplePos, const glm::vec3& debugColor, const BvhInterface& bvh, const Features& features, Ray ray, HitInfo hitInfo)
{
    glm::vec3 intersectionPosition = ray.origin + ray.direction * ray.t;
    // Shadow raydirection = from the position of light - to the position of intersection point
    glm::vec3 shadowRayDirection = samplePos - intersectionPosition;

    float epsilon = 10e-4f;
    Ray shadowRay = Ray(intersectionPosition + epsilon * shadowRayDirection, shadowRayDirection);

    bool inShadow = bvh.intersect(shadowRay, hitInfo, features);

    if (inShadow && shadowRay.t > 0.0f && shadowRay.t < 1.0f) {
        drawRay(shadowRay, glm::vec3(1.0f, 0.0f, 0.0f)); // draw a red ray from the light to the point in the color of the pixel
        return 0.0f;
    } else {
        shadowRay.t = 1.0f; // needed for the visual debug - to make the ray connect the light and the intersection point
        drawRay(shadowRay, debugColor); // draw a ray from the light to the point in the color of the pixel
        return 1.0f;
    }
}


// given an intersection, computes the contribution from all light sources at the intersection point
// in this method you should cycle the light sources and for each one compute their contribution
// don't forget to check for visibility (shadows!)
//
// Regarding the soft shadows for **other** light sources **extra** feature:
// To add a new light source, define your new light struct in scene.h and modify the Scene struct (also in scene.h)
// by adding your new custom light type to the lights std::variant. For example:
// std::vector<std::variant<PointLight, SegmentLight, ParallelogramLight, MyCustomLightType>> lights;
//
// You can add the light sources programmatically by creating a custom scene (modify the Custom case in the
// loadScene function in scene.cpp). Custom lights will not be visible in rasterization view.
glm::vec3 computeLightContribution(const Scene& scene, const BvhInterface& bvh, const Features& features, Ray ray, HitInfo hitInfo)
{
    if (features.enableShading) {
        // If shading is enabled, compute the contribution from all lights.

        glm::vec3 color = { 0.0f, 0.0f, 0.0f };

        for (const auto& light : scene.lights) {
            if (std::holds_alternative<PointLight>(light)) {
                const PointLight pointLight = std::get<PointLight>(light);
                if (features.enableHardShadow) {
                    // Perform your calculations for a point light.
                    if (testVisibilityLightSample(pointLight.position, pointLight.color, bvh, features, ray, hitInfo) == 1.0f) {      //Only if point light is visible, compute color
                        color += computeShading(pointLight.position, pointLight.color, features, ray, hitInfo);                
                    }
                } else {
                    color += computeShading(pointLight.position, pointLight.color, features, ray, hitInfo);                           //If hard shadows disabled, point light should contribute to all positions
                }
            } else if (std::holds_alternative<SegmentLight>(light)) {
                if (features.enableSoftShadow) {
                    int segmentLightSamples = 50;                  //Number of samples used to calculate light -> More samples = More render time
                    const SegmentLight segmentLight = std::get<SegmentLight>(light);
                    glm::vec3 sum = { 0.0f, 0.0f, 0.0f };
                    // Partition segment in equal lengths -> Used for jittered sampling (Ch 13 Book)
                    std::vector<glm::vec2> weightsList;
                    float weightLength = 1.0f / segmentLightSamples;           // How much distance will each partition of segment have
                    // Populate list of weights (Sizes of segments for jittered sampling)
                    for (float i = 0; i < segmentLightSamples; i++) {
                        float minWeight = i / segmentLightSamples;
                        float maxWeight = minWeight + weightLength;
                        weightsList.push_back({ minWeight, maxWeight });
                    }
                    // For each partition in segment, sample a random point and compute color
                    for (glm::vec2 weights : weightsList) {
                        glm::vec3 samplePos = { weights.x, weights.y, 0.0f };
                        glm::vec3 sampleColor = { 0.0f, 0.0f, 0.0f };
                        sampleSegmentLight(segmentLight, samplePos, sampleColor);    //Get random point inside parititon
                        //visualDebugSoftShadows(samplePos, sampleColor, bvh, features, ray, hitInfo);
                        if (testVisibilityLightSample(samplePos, sampleColor, bvh, features, ray, hitInfo) == 1) {   //Can also be used as visual debug
                            glm::vec3 shadedColor = computeShading(samplePos, sampleColor, features, ray, hitInfo); // If visible, compute shading for sample & add to sum
                            sum += glm::vec3 { shadedColor.x, shadedColor.y, shadedColor.z };
                        }
                    }
                    color += glm::vec3 { sum.x / segmentLightSamples, sum.y / segmentLightSamples, sum.z / segmentLightSamples };     //Final Color = Average of all colors
                }
            } else if (std::holds_alternative<ParallelogramLight>(light)) {
                if (features.enableSoftShadow) {
                    int parallelogramLightSamples = 5;
                    const ParallelogramLight parallelogramLight = std::get<ParallelogramLight>(light);
                    glm::vec3 sum = { 0.0f, 0.0f, 0.0f };
                    std::vector<glm::vec2> edge01Weights;
                    std::vector<glm::vec2> edge02Weights;
                    float weightLength = 1.0f / parallelogramLightSamples;
                    // Populate list of weights (Sizes of grid for jittered sampling)
                    for (float i = 0; i < parallelogramLightSamples; i++) {
                        float minWeight = i / parallelogramLightSamples;
                        float maxWeight = minWeight + weightLength;
                        edge01Weights.push_back({ minWeight, maxWeight });
                        edge02Weights.push_back({ minWeight, maxWeight });
                    }
                    // For each grid in parallelogram, sample a point and compute color
                    for (glm::vec2 edge01Weight : edge01Weights) {
                        for (glm::vec2 edge02Weight : edge02Weights) {
                            glm::vec3 samplePos = glm::vec3 { edge01Weight.x, edge01Weight.y, 0.0f };
                            glm::vec3 sampleColor = { edge02Weight.x, edge02Weight.y, 0.0f };
                            sampleParallelogramLight(parallelogramLight, samplePos, sampleColor);      //Get random point inside grid
                            // visualDebugSoftShadows(samplePos, sampleColor, bvh, features, ray, hitInfo);
                            if (testVisibilityLightSample(samplePos, sampleColor, bvh, features, ray, hitInfo) == 1) {    //Can also be used as visual debug
                                glm::vec3 shadedColor = computeShading(samplePos, sampleColor, features, ray, hitInfo); // If visible, compute shading for sample & add to sum
                                sum += glm::vec3 { shadedColor.x, shadedColor.y, shadedColor.z };
                            }
                        }
                    }
                    color += glm::vec3 { sum.x / parallelogramLightSamples, sum.y / parallelogramLightSamples, sum.z / parallelogramLightSamples };   //Final Color = Average of all colors
                }
            }
        }
        return color;
    } else {
        // If shading is disabled, return the albedo of the material. ->
        return hitInfo.material.kd;
    }
}

