#pragma once

#include <exception>
#include <string>
#include "macros.hpp"
#include <glm/glm.hpp>

inline std::string plural (std::string word) {
    if (word.size() == 0)
      return word;
    bool caps = word.back() >= 'A' && word.back() <= 'Z' && word.size() > 1;
    bool es = word.back() == 's'||
              word.back() == 'x' ||
              word.back() == 'z' ||
              (word.size() > 1 && word[word.size()-2] == 'c' && word.back() == 'h') ||
              (word.size() > 1 && word[word.size()-2] == 's' && word.back() == 'h') ||
              word.back() == 'S'||
              word.back() == 'X' ||
              word.back() == 'Z' ||
              (word.size() > 1 && word[word.size()-2] == 'C' && word.back() == 'H') ||
              (word.size() > 1 && word[word.size()-2] == 'S' && word.back() == 'H');
    if (es && !caps)
      return word + "es";
    if (es)
      return word + "ES";
    if (caps)
      return word + "S";
    return word + "s";
}

class NotImplementedError : public std::exception {
    std::string msg_;
public:
    explicit NotImplementedError(const std::string& notImplementedMethodName, const std::string& lackingType="Method")
        : msg_(lackingType + " " + notImplementedMethodName + " is not implemented yet.") {}

    const char* what() const noexcept override {
        return msg_.c_str();
    }
};

class NotImplementedMethodError : public NotImplementedError {
public:
  explicit NotImplementedMethodError(const std::string& methodName)
      :NotImplementedError(methodName, "Method") {}
};

class NotImplementedFunctionError : public NotImplementedError {
public:
  explicit NotImplementedFunctionError(const std::string& name)
      :NotImplementedError(name, "Function") {}
};

class NotImplementedVariantError : public NotImplementedError {
public:
  NotImplementedVariantError(const std::string& variant, const std::string& ofWhat)
      :NotImplementedError(variant + " of " + ofWhat, "Variant") {}
};

class UnknownVariantError : public std::exception {
  std::string msg_;
public:
  UnknownVariantError(const std::string& variant, const std::string& ofWhat)
      : msg_("Variant " + variant + " of  " + ofWhat + " is an unknown type in this context.") {}
  explicit UnknownVariantError(const std::string& msg)
      : msg_(msg) {}


  const char* what() const noexcept override {
    return msg_.c_str();
  }
};

class IllegalVariantError : public std::exception {
  std::string msg_;
public:
  IllegalVariantError(const std::string& variant, const std::string& ofWhat, const std::string& rejectingMethod)
      : msg_(plural(ofWhat) + " in variant " + variant + " are considered invalid by method " + rejectingMethod + ".") {}
  explicit IllegalVariantError(const std::string& msg)
      : msg_(msg) {}

  const char* what() const noexcept override {
    return msg_.c_str();
  }
};

class COLOR_PALETTE {
public:
    glm::vec4 mainColor;
    glm::vec4 second;
    glm::vec4 third;
    glm::vec4 accent;
    glm::vec4 accent2;

    COLOR_PALETTE(glm::vec4 mainColor, glm::vec4 second, glm::vec4 third, glm::vec4 accent, glm::vec4 accent2);
    COLOR_PALETTE(glm::vec3 mainColor, glm::vec3 second, glm::vec3 third, glm::vec3 accent, glm::vec3 accent2);
    COLOR_PALETTE(glm::ivec3 mainColor, glm::ivec3 second, glm::ivec3 third, glm::ivec3 accent, glm::ivec3 accent2);
    std::vector<glm::vec4> colors();
    glm::vec4 operator[] (int i);

};

class COLOR_PALETTE10 {
public:
    std::array<glm::vec4, 10> cls;

    explicit COLOR_PALETTE10(std::array<glm::vec4, 10> colors) : cls(colors) {}
    COLOR_PALETTE10(COLOR_PALETTE p1, COLOR_PALETTE p2) : cls({p1[0], p1[1], p1[2], p1[3], p1[4], p2[0], p2[1], p2[2], p2[3], p2[4]}) {}
    COLOR_PALETTE10(glm::ivec3 c1, glm::ivec3 c2, glm::ivec3 c3, glm::ivec3 c4, glm::ivec3 c5, glm::ivec3 c6, glm::ivec3 c7, glm::ivec3 c8, glm::ivec3 c9, glm::ivec3 c10);
    std::vector<glm::vec4> colors() const { return std::vector<glm::vec4>(cls.begin(), cls.end()); }
    glm::vec4 operator[] (int i) const { return cls[i]; }
};


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

    template<typename T>
    void printVector(std::vector<T> v, std::string title="vector")
    {
        std::cout << title << " [" << v.size() << "]: ";
        for (int i = 0; i < v.size(); i++)
        {
            std::cout << v[i] << ", ";
        }
        std::cout << std::endl;
    }


const PolyGroupID DEFAULT_POLY_GROUP_ID = PolyGroupID(0);
constexpr RP1 inf = std::nullopt;
constexpr RP1 unbounded = std::nullopt;
const glm::vec3 e1 = glm::vec3(1, 0, 0);
const glm::vec3 e2 = glm::vec3(0, 1, 0);
const glm::vec3 e3 = glm::vec3(0, 0, 1);
const glm::vec3 ORIGIN = glm::vec3(0, 0, 0);
const glm::vec2 PLANE_ORIGIN = glm::vec2(0, 0);
const PolyGroupID DFLT_CURV = PolyGroupID(420);
