#version 330 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aColor;
layout (location = 2) in vec2 aTex;

out vec4 fragColor;
out vec2 texCoord;

void main() {
    gl_Position = vec4(aPos.x, aPos.y, aPos.z, 1.0);
    fragColor = vec4(aColor.x, aColor.y, aColor.z, 1.0f);
    texCoord = aTex;
}