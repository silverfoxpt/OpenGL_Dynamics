#ifndef _TEXTURE_2D_H_
#define _TEXTURE_2D_H_

#include <iostream>
#include <fstream>
#include <sstream>

#include <glad/glad.h> 
#include <GLFW/glfw3.h>

#include "stb_image.h"

class Texture2DCustom {
    public:
        unsigned int textureId;
        int textureUnitIndex = 0;

        Texture2DCustom(const std::string& imageFilePath, int unitIndex);\
        Texture2DCustom() { }

        void UseTexture();

    private:
};

#endif