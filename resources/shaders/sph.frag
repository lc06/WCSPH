#version 330 core

uniform vec2 u_resolution;
uniform float u_time;

in vec3 ourColor;

float circle(vec2 st, float radius) {
    vec2 distVec = st - vec2(0.5);

    float distValue = 1. - step(radius, dot(distVec, distVec) * 4.);
    return distValue;
}

void main() {
//    vec2 st = gl_FragCoord.xy / u_resolution;
//
//    float radius = 0.25;
//    vec3 color = mix(ourColor, vec3(0., 0., 0.), circle(st, radius));

    gl_FragColor = vec4(ourColor, 1.0f);
}