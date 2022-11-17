#include "texture.h"
#include <framework/image.h>

glm::vec3 acquireTexel(const Image& image, const glm::vec2& texCoord, const Features& features)
{
    if (!features.enableTextureMapping)
        return image.pixels[0];
    
    int maxIndex = image.pixels.size() - 1;
    float xCoord = texCoord.x * image.width - 0.5f;
    float yCoord = (1.0f - texCoord.y) * image.height - 0.5f;
    if (features.extra.enableBilinearTextureFiltering) {
        //Get closest coordinates
        int leftXcoord = floor(xCoord);
        int rightXCoord = (leftXcoord + 1);
        int lowerYCoord = floor(yCoord);
        int upperYCoord = (lowerYCoord + 1);
        //Get weights for closest coordinates
        float alphaX = (rightXCoord - xCoord);
        float betaX = 1.0f - alphaX;
        float alphaY = (upperYCoord - yCoord);
        float betaY = 1.0f - alphaY;
        //Return bilinear interpolation of coordinates using calculated weights using bilinear interpolation formula
        return alphaX * alphaY * image.pixels[std::clamp((lowerYCoord * image.width + leftXcoord), 0, maxIndex)] + alphaX * betaY * image.pixels[std::clamp((upperYCoord * image.width + leftXcoord), 0, maxIndex)]
            + betaX * alphaY * image.pixels[std::clamp((lowerYCoord * image.width + rightXCoord), 0, maxIndex)] + betaX * betaY * image.pixels[std::clamp((upperYCoord * image.width + rightXCoord), 0, maxIndex)];
    }
    int i = std::round(xCoord);
    int j = std::round(yCoord);
    return image.pixels[std::clamp((j * image.width + i), 0, maxIndex)];
    // Note, the center of the first pixel is at image coordinates (0.5, 0.5)
}