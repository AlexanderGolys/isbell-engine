#version 450 core
precision mediump float;

in vec3 v_position;
in vec3 v_normal;
in vec2 v_uv;
in vec4 v_color;

layout(location = 0) out vec4 color;

uniform float u_time;
uniform mat4 u_mvp;
uniform vec3[3] u_lightPos;
uniform vec4[3] u_lightColor;
uniform vec3[3] u_lightIntencities;
uniform vec4 u_intencities;
uniform sampler2D u_texture_ambient;
uniform sampler2D u_texture_diffuse;
uniform sampler2D u_texture_specular;
uniform vec3 u_camPosition;


vec4 saturate(vec4 color) {
    return clamp(color, 0., 1.);
}

float saturate(float x) {
    return clamp(x, 0., 1.);
}



vec4 lightFactor(int i, vec3 position) {
    vec3 p = u_lightPos[i];
    float distance = length(p - position);
    vec4 col = u_lightColor[i];
    float constant = u_lightIntencities[i].x;
    float linear = u_lightIntencities[i].y;
    float quadratic = u_lightIntencities[i].z;
    return col / (constant + linear * distance + quadratic * distance * distance);
}



vec4 colorFromPointlight(int i, vec3 camPos, vec3 pos, vec3 n, vec2 uv)
{
    vec3 lightPosition = u_lightPos[i];
	vec3 light_direction = normalize(lightPosition - pos);
//    float cosTheta = abs(dot(n, light_direction));
    float cosTheta = max(dot(n, -light_direction), 0.);

	vec3 reflectDirection = normalize(reflect(-light_direction, n));
	vec3 viewDirection = normalize(camPos - pos);
	float cosAlpha = max(dot(viewDirection, reflectDirection), 0.);
//    float cosAlpha = abs(dot(viewDirection, reflectDirection));
    vec4 light = u_lightColor[i];

    return  texture(u_texture_ambient,  uv) * u_intencities.x +
            texture(u_texture_diffuse,  uv) * u_intencities.y * cosTheta * lightFactor(i, pos) +
            texture(u_texture_specular,  uv) * u_intencities.z  * pow(cosAlpha, u_intencities.w) * lightFactor(i, pos);
}



void main()
{
	vec3 normal = v_normal/length(v_normal);
	color = saturate(colorFromPointlight(0, u_camPosition, v_position, normal, v_uv)) +
              saturate(colorFromPointlight(1, u_camPosition, v_position, normal, v_uv)) +
              saturate(colorFromPointlight(2, u_camPosition, v_position, normal, v_uv));
    color.a = 1;
}
