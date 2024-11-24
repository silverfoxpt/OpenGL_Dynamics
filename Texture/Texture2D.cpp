#include "Texture2D.h"

Texture2DCustom::Texture2DCustom(const std::string& imageFilePath, int unitIndex) {
    this->textureUnitIndex = unitIndex;

    // Generate texture
    glGenTextures(1, &this->textureId);
    glActiveTexture(GL_TEXTURE0 + this->textureUnitIndex);
    glBindTexture(GL_TEXTURE_2D, this->textureId);

    // Load images with stb_image.h functions
    int width, height, nrChannels;
    unsigned char *data = stbi_load(imageFilePath.c_str(), &width, &height, &nrChannels, 0);

    // Check data alright
    if (data) {
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, data);
        glGenerateMipmap(GL_TEXTURE_2D);
    } 
    else {
        std::cerr << "Error: Texture not loadable! Please recheck.";
        return;
    }
    stbi_image_free(data);

    // Set parameters for texture - maybe separate function?
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
}

void Texture2DCustom::UseTexture() {
    glActiveTexture(GL_TEXTURE0 + this->textureUnitIndex);
    glBindTexture(GL_TEXTURE_2D, this->textureId);
}