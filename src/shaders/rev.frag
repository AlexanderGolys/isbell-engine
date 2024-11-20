#version 330 core

in vec3 v_position;
in vec3 v_normal;
in vec2 v_uv;
in vec4 v_color;

layout(location = 0) out vec4 color;

uniform mat4 mvp;
uniform mat4 light1;
uniform mat4 light2;
uniform mat4 light3;
uniform vec3 camPosition;
uniform vec4 intencities;
uniform sampler2D texture_ambient;
uniform sampler2D texture_diffuse;
uniform sampler2D texture_specular;
uniform float time;
uniform vec4 bdColor;


float TAU = 6.28318530718;

float mod(float a, float b)
{
    return a - b * floor(a / b);
}

vec3 lightPosition(mat4 m)
{
	return (m*vec4(1., 0., 0., 0.)).xyz;
}

float lightPower(mat4 m)
{
	return (m*vec4(1., 0., 0., 0.)).w;
}

vec4 lightColor(mat4 m)
{
	return m*vec4(0., 1., 0., 0.);
}

vec4 colorFactorFromSingleLight(mat4 light, vec3 camPos, vec3 pos, vec3 n, vec2 uv)
{
	vec3 light_direction = normalize(lightPosition(light) - pos);
//    float cosTheta = abs(dot(n, light_direction));

    float cosTheta = max(abs(dot(n, light_direction)), 0.);
	float distanceSquared = dot(lightPosition(light) - pos, lightPosition(light) - pos);
	vec3 reflectDirection = normalize(reflect(-light_direction, n));
	vec3 viewDirection = normalize(camPos - pos);
	float cosAlpha = max(abs(dot(viewDirection, reflectDirection)), 0.);

    return  texture(texture_ambient, uv) * intencities.x +
    texture(texture_diffuse, uv) * intencities.y * cosTheta * lightPower(light) * lightColor(light) / distanceSquared+
    texture(texture_specular, uv) * intencities.z  * pow(cosAlpha, intencities.w) * lightColor(light) * lightPower(light) / distanceSquared;

}



void main()
{
    float t = time;
	vec3 normal = normalize(v_normal);
	if (v_uv.y >= t)
    {
       discard;
    }
	color = colorFactorFromSingleLight(light1, camPosition, v_position, normal, v_uv) +
			colorFactorFromSingleLight(light2, camPosition, v_position, normal, v_uv) +
			colorFactorFromSingleLight(light3, camPosition, v_position, normal, v_uv);

    color.a = 1.;
    float collar = smoothstep(-.1, 0, v_uv.y -t);
    color = mix(color, bdColor, collar);
//    collar = 1-smoothstep(0, .1, v_uv.y);
//    color = mix(color, bdColor, collar);

}
