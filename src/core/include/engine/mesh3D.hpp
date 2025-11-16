#pragma once
#include "indexedMesh.hpp"
#include "mat.hpp"
#include "smoothParametric.hpp"


struct Vertex3D {
	vec3 position;
	vec3 normal;
	vec2 uv;
	vec4 color;

	Vertex3D(vec3 position, vec3 normal, vec2 uv, vec4 color);

	raw_data_ptr data() const;

	static VertexBufferLayout layout();
};

class Object3D {
	PROPERTY(SE3, transform);
public:
	void translate(const vec3& t);
	void apply(const SE3& T);
};

class Mesh3D : public IndexedMesh<Vertex3D>, public Object3D {
public:
	using IndexedMesh::IndexedMesh;
	Mesh3D(const SmoothParametricSurface& surf, int uRes, int vRes);
};
