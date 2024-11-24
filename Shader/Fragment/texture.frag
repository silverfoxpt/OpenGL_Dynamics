#version 330 core
in vec4 fragColor;
in vec2 texCoord;

out vec4 FragColor;

uniform sampler2D ourTexture;
uniform sampler2D smileTexture;

void main() {
    //FragColor = vec4(0.15, 0.87, 0.63, 1.0);
    //FragColor = fragColor;
    //FragColor = texture(ourTexture, texCoord) * fragColor;
    FragColor = mix(texture(ourTexture, texCoord), texture(smileTexture, texCoord), 0.5) * fragColor;
}