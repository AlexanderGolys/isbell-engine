#pragma once

#include <glm/glm.hpp>

namespace glm {
    inline vec3 xyz(const vec4& v) {
        return vec3(v.x, v.y, v.z);
    }
    inline vec3 yzw(const vec4& v) {
        return vec3(v.y, v.z, v.w);
    }
    inline vec2 xy(const vec4& v) {
        return vec2(v.x, v.y);
    }
    inline vec2 yz(const vec4& v) {
        return vec2(v.y, v.z);
    }
    inline vec2 zw(const vec4& v) {
        return vec2(v.z, v.w);
    }
    inline vec2 xy(const vec3& v) {
        return vec2(v.x, v.y);
    }
    inline vec2 yz(const vec3& v) {
        return vec2(v.y, v.z);
    }

    inline vec4 xyz1(const vec3& v) {
        return vec4(v, 1);
    }


}