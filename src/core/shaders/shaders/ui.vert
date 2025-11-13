#version 450 core
precision mediump float;

layout (location = 0) in vec2 position;
layout (location = 1) in vec2 uv;

out vec2 v_position;
out vec2 v_uv;

uniform sampler2D u_texture;

void main()
{
	v_position = position;
    v_uv = uv;
    gl_Position = vec4(v_position, 0.0, 1.0);
}
