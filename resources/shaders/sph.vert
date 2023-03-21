#version 330 core
layout (location = 0) in vec2 Pos;   // 位置变量的属性位置值为 0
layout (location = 1) in float velocity; // 速度变量的属性位置值为 1
layout (location = 2) in float type; // 类型变量的属性位置值为 1

uniform mat4 model;

vec3 black = vec3(0.0, 0.0, 0.0);
vec3 green = vec3(0.0, 1.0, 0.0);
vec3 white = vec3(1.0, 1.0, 1.0);
vec3 red = vec3(1.0, 0.0, 0.0);

out vec3 ourColor; // 向片段着色器输出一个颜色

void main()
{
    vec4 homo_pos = vec4(Pos, 0.0, 1.0);
    gl_Position = model * homo_pos;

    ourColor = mix(mix(mix(black, green, step(1., type)), white, step(2., type)), red, velocity); // 将ourColor设置为我们从顶点数据那里得到的输入颜色
}