#include "dynamicMesh.hpp"

Vertex3D::Vertex3D(): position(vec3(0, 0, 0)), normal(vec3(0, 0, 1)), uv(vec2(0, 0)), color(vec4(1, 1, 1, 1)) {}

Vertex3D::Vertex3D(const vec3& position, const vec3& normal, const vec2& uv, const vec4& color): position(position), normal(normal), uv(uv), color(color) {}

Vertex2D::Vertex2D(): position(vec2(0, 0)), uv(vec2(0, 0)), color(vec4(1, 1, 1, 1)) {}

Vertex2D::Vertex2D(const vec2& position, const vec2& uv, const vec4& color): position(position), uv(uv), color(color) {}
