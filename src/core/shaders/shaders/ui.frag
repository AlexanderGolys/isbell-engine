#version 450 core
precision mediump float;

in vec2 v_position;
in vec2 v_uv;

layout(location = 0) out vec4 color;

uniform sampler2D u_texture;

void main()
{
    color = texture(u_texture, v_uv);
}
