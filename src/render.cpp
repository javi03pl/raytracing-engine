#include "render.h"
#include "intersect.h"
#include "light.h"
#include "screen.h"
#include <framework/trackball.h>
#ifdef NDEBUG
#include <omp.h>
#endif

glm::vec3 getFinalColor(const Scene& scene, const BvhInterface& bvh, Ray ray, const Features& features, int rayDepth)
{
    HitInfo hitInfo;
    if (bvh.intersect(ray, hitInfo, features)) {

        glm::vec3 Lo = computeLightContribution(scene, bvh, features, ray, hitInfo);

        if (features.enableRecursive && rayDepth <= 10) {
            Ray reflection = computeReflectionRay(ray, hitInfo);
            // TODO: put your own implementation of recursive ray tracing here.
            // move origin a bit along the direction of the ray, so it is not behind the plane it originates from
            float epsilon = 10e-4f;
            reflection = Ray(reflection.origin + epsilon * reflection.direction, reflection.direction);

            if (hitInfo.material.ks.x != 0.0 || hitInfo.material.ks.y != 0.0 || hitInfo.material.ks.x != 0.0) {
                bvh.intersect(reflection, hitInfo, features);
                glm::vec3 additionalColor = getFinalColor(scene, bvh, reflection, features, rayDepth + 1);
                Lo += additionalColor;
            }
        }

        // Visual debug for shading and recursive ray-tracing
        drawRay(ray, Lo);
        
        // Set the color of the pixel to white if the ray hits.
        return Lo;

    } else {
        // Draw a red debug ray if the ray missed.
        drawRay(ray, glm::vec3(1.0f, 0.0f, 0.0f));
        // Set the color of the pixel to black if the ray misses.
        return glm::vec3(0.0f);
    }
}

//Checks if a set of coordinates are inside the bounds of the screen
bool insideScreen(glm::vec2 coords, Screen& screen) {
    if (coords.x < 0 || coords.x >= screen.resolution().x || coords. y < 0 || coords.y >= screen.resolution().y) return false;
    return true;
}

//Calculates a final color by using the gaussian distribution values for a kernel size of 21 (Weights pre computed for performance reasons)
glm::vec3 blurFragment(glm::vec2 blurTextureCoords[], Screen& screen, int filterSize)
{
    glm::vec3 res = { 0.0f, 0.0f, 0.0f };
    for (int i = 0; i < filterSize; i++) {
        if (!insideScreen(blurTextureCoords[i], screen)) {              //If a pixel is outside the screen, its color will be set to that of the center pixel in the image
            blurTextureCoords[i] = blurTextureCoords[(filterSize - 1) / 2];
        }
    }
    //Multiply the colors of each coordinate in texCoords by its corresponding weight inside gaussian distribution
    res += screen.pixels()[screen.indexAt(blurTextureCoords[0].x, blurTextureCoords[0].y)] * 0.0005f;
    res += screen.pixels()[screen.indexAt(blurTextureCoords[1].x, blurTextureCoords[1].y)] * 0.0015f;
    res += screen.pixels()[screen.indexAt(blurTextureCoords[2].x, blurTextureCoords[2].y)] * 0.0039f;
    res += screen.pixels()[screen.indexAt(blurTextureCoords[3].x, blurTextureCoords[3].y)] * 0.0089f;
    res += screen.pixels()[screen.indexAt(blurTextureCoords[4].x, blurTextureCoords[4].y)] * 0.0183f;
    res += screen.pixels()[screen.indexAt(blurTextureCoords[5].x, blurTextureCoords[5].y)] * 0.0334f;
    res += screen.pixels()[screen.indexAt(blurTextureCoords[6].x, blurTextureCoords[6].y)] * 0.0549f;
    res += screen.pixels()[screen.indexAt(blurTextureCoords[7].x, blurTextureCoords[7].y)] * 0.0807f;
    res += screen.pixels()[screen.indexAt(blurTextureCoords[8].x, blurTextureCoords[8].y)] * 0.1063f;
    res += screen.pixels()[screen.indexAt(blurTextureCoords[9].x, blurTextureCoords[9].y)] * 0.1253f;
    res += screen.pixels()[screen.indexAt(blurTextureCoords[10].x, blurTextureCoords[10].y)] * 0.1324f;
    res += screen.pixels()[screen.indexAt(blurTextureCoords[11].x, blurTextureCoords[11].y)] * 0.1253f;
    res += screen.pixels()[screen.indexAt(blurTextureCoords[12].x, blurTextureCoords[12].y)] * 0.1063f;
    res += screen.pixels()[screen.indexAt(blurTextureCoords[13].x, blurTextureCoords[13].y)] * 0.0807f;
    res += screen.pixels()[screen.indexAt(blurTextureCoords[14].x, blurTextureCoords[14].y)] * 0.0549f;
    res += screen.pixels()[screen.indexAt(blurTextureCoords[15].x, blurTextureCoords[15].y)] * 0.0334f;
    res += screen.pixels()[screen.indexAt(blurTextureCoords[16].x, blurTextureCoords[16].y)] * 0.0183f;
    res += screen.pixels()[screen.indexAt(blurTextureCoords[17].x, blurTextureCoords[17].y)] * 0.0089f;
    res += screen.pixels()[screen.indexAt(blurTextureCoords[18].x, blurTextureCoords[18].y)] * 0.0039f;
    res += screen.pixels()[screen.indexAt(blurTextureCoords[19].x, blurTextureCoords[19].y)] * 0.0015f;
    res += screen.pixels()[screen.indexAt(blurTextureCoords[20].x, blurTextureCoords[20].y)] * 0.0005f;
    return res;
}

//Returns a horizontaly blurred image using gaussian blur
Screen hBlur(Screen& screen)
{
    const int filterSize = 21; //Size of kernel for gaussian distribution
    glm::vec2 blurTextureCoords[filterSize]; 
    Screen res = Screen(screen.resolution(), false);
    int arrayLimits = filterSize / 2;
    for (int y = 0; y < screen.resolution().y; y++) {
        for (int x = 0; x != screen.resolution().x; x++) {
            for (int i = -arrayLimits; i <= arrayLimits; i++) {
                blurTextureCoords[i + arrayLimits] = glm::vec2 { x + i, y }; //Add texture coordinates that will be used for filter
            }
            res.setPixel(x, y, blurFragment(blurTextureCoords, screen, filterSize)); //Calculate filter for each pixel inside original screen
        }
    }
    return res;
}

// Returns a vertically blurred image using gaussian blur
Screen vBlur(Screen& screen)
{
    const int filterSize = 21; //Size of kernel for gaussian distribution
    glm::vec2 blurTextureCoords[filterSize];
    Screen res = Screen(screen.resolution(), false);
    int arrayLimits = filterSize / 2;
    for (int y = 0; y < screen.resolution().y; y++) {
        for (int x = 0; x != screen.resolution().x; x++) {
            for (int i = -arrayLimits; i <= arrayLimits; i++) {
                blurTextureCoords[i + arrayLimits] = glm::vec2 { x, y + i }; //Add texture coordinates that will be used for filter
            }
            res.setPixel(x, y, blurFragment(blurTextureCoords, screen, filterSize)); //Calculate filter for each pixel inside original screen
        }
    }
    return res;
}

int bloomIntensity; //Intensity of bloom fiter -> How many blur passes should be done
bool bloomDebug; //If true, only the bloom layer will appear when rendering images

void setBloomDebug(bool value) {
    bloomDebug = value;
}

void setBloomIntensity(int value) {
    bloomIntensity = value;
}

void renderRayTracing(const Scene& scene, const Trackball& camera, const BvhInterface& bvh, Screen& screen, const Features& features)
{
    glm::ivec2 windowResolution = screen.resolution();
    Screen thresholdedImage = Screen(windowResolution, false);

    // Enable multi threading in Release mode
#ifdef NDEBUG
#pragma omp parallel for schedule(guided)
#endif
    for (int y = 0; y < windowResolution.y; y++) {
        for (int x = 0; x != windowResolution.x; x++) {
            // NOTE: (-1, -1) at the bottom left of the screen, (+1, +1) at the top right of the screen.
            const glm::vec2 normalizedPixelPos {
                float(x) / float(windowResolution.x) * 2.0f - 1.0f,
                float(y) / float(windowResolution.y) * 2.0f - 1.0f
            };
            const Ray cameraRay = camera.generateRay(normalizedPixelPos);
            glm::vec3 computedColor = getFinalColor(scene, bvh, cameraRay, features);
            screen.setPixel(x, y, computedColor);

            //If bloom enabled, get the thresholded image used for blurring
            if (features.extra.enableBloomEffect) {
                float brightness = (computedColor.x * 0.2126f) + (computedColor.y * 0.7152f) + (computedColor.z * 0.0722f); //Using luma brightness formula to calculate brightness (Uses relative luminance insted of a threshold)
                thresholdedImage.setPixel(x, y, computedColor * brightness); //Each color will be multiplied by its brightness in the thresholded image
            }
        }
    }
    
    if (features.extra.enableBloomEffect) {
        // Calculate blurred thresholded image
        Screen horizontalBlur = hBlur(thresholdedImage);
        Screen verticalBlur = hBlur(horizontalBlur);
        for (int i = 0; i < bloomIntensity; i++) { //More intensity -> More blur passes -> More rendering time
            horizontalBlur = hBlur(verticalBlur);
            verticalBlur = vBlur(horizontalBlur);
        }


        //Sum both images
        for (int y = 0; y < windowResolution.y; y++) {
            for (int x = 0; x != windowResolution.x; x++) {
                if (bloomDebug) {
                    screen.setPixel(x, y, verticalBlur.pixels()[verticalBlur.indexAt(x, y)]);    //Only return bloom layer for debugging
                } else {
                    if (bloomIntensity > 0) {
                        screen.setPixel(x, y, screen.pixels()[screen.indexAt(x, y)] + verticalBlur.pixels()[verticalBlur.indexAt(x, y)]); //Sum original image & bloom layer to get final image
                    }  
                }
            }
        }
    }
}