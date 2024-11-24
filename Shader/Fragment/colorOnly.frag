#version 330 core
in vec4 fragColor;
out vec4 FragColor;

void main() {
    //FragColor = vec4(0.15, 0.87, 0.63, 1.0);
    FragColor = fragColor;
}