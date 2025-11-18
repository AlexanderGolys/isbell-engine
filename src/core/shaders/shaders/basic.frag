#version 450 core
precision mediump float;

in vec3 v_position;
in vec3 v_normal;
in vec2 v_uv;
in vec4 v_color;

layout(location = 0) out vec4 color;

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

layout(std140, binding = 1) uniform Camera {
    vec3 u_camPosition;
    mat4 u_vp;
};
layout(std140, binding = 2) uniform Model {
    mat4 u_m;
};

layout(std140, binding = 3) uniform Time {
    float u_t;
    float u_dt;
};
layout(std140, binding = 4) uniform Material {
    vec4 u_intencities;
};

uniform sampler2D u_texture_ambient;
uniform sampler2D u_texture_diffuse;
uniform sampler2D u_texture_specular;


vec4 saturate(vec4 c) {
    return clamp(c, 0., 1.);
}

float saturate(float x) {
    return clamp(x, 0., 1.);
}



vec4 lightFactor(int i, vec3 position) {
    float distance = length(u_lights[i].pos - position);
    LightIntencity inten = u_lights[i].intencity;
    return u_lights[i].color / (inten.constant + inten.linear * distance + inten.quadratic * distance * distance);
}



vec4 colorFromPointlight(vec3 pos, vec3 n, vec2 uv)
{
    vec4 result = vec4(0.);
    for (int i = 0; i < u_lights.length(); i++) {
        vec3 light_direction = normalize(u_lights[i].pos - pos);
        float cosTheta = max(dot(n, light_direction), 0.);

        vec3 reflectDirection = normalize(reflect(light_direction, n));
        vec3 viewDirection = normalize(u_camPosition - pos);
        float cosAlpha = max(dot(viewDirection, reflectDirection), 0.);
        if (cosTheta <= 0.0001)
            cosAlpha = 0.;

        result += saturate(texture(u_texture_ambient,  uv) * u_intencities.x +
                texture(u_texture_diffuse, uv) * u_intencities.y * cosTheta * lightFactor(i, pos) +
                texture(u_texture_specular, uv) * u_intencities.z  * pow(cosAlpha, u_intencities.w) * lightFactor(i, pos)) ;
    }
    return saturate(result);
}



void main()
{
	vec3 normal = normalize(v_normal);
	color = colorFromPointlight(v_position, normal, v_uv);
    color.a = 1;
    color.r = .9;
}
