#version 450 core
precision mediump float;

layout (location = 0) in vec3 position;
layout (location = 1) in vec3 normal;
layout (location = 2) in vec2 uv;
layout (location = 3) in vec4 color;

out vec4 v_color;
out vec3 v_position;
out vec3 v_normal;
out vec2 v_uv;

struct LightIntencity {
    float constant;
    float linear;
    float quadratic;
};

struct Light {
    vec3 pos;
    vec4 color;
    LightIntencity intencity;
};

layout(std140, binding = 0) uniform Lights {
    Light u_lights[3];
};

uniform mat4 u_mvp;
uniform vec4 u_intencities;
uniform sampler2D u_texture_ambient;
uniform sampler2D u_texture_diffuse;
uniform sampler2D u_texture_specular;
uniform float u_time;
uniform vec3 u_camPosition;

vec3 normalise(vec3 v)
{
	return v / length(v);
}

void main()
{
	v_normal = normal;
	v_position = position;
	v_uv = uv;
	gl_Position = (u_mvp * vec4(v_position, 1.0));
	v_color = color;
}
