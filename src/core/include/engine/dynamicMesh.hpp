#pragma once

#include "buffers.hpp"



struct Vertex3D {
	vec3 position;
	vec3 normal;
	vec2 uv;
	vec4 color;

	Vertex3D();
	Vertex3D(const vec3& position, const vec3& normal, const vec2& uv, const vec4& color);
};

struct Vertex2D {
	vec2 position;
	vec2 uv;
	vec4 color;

	Vertex2D();
	Vertex2D(const vec2& position, const vec2& uv, const vec4& color);
};


template<typename V>
struct Triangle {
	V a, b, c;
	Triangle(const V& a, const V& b, const V& c) : a(a), b(b), c(c) {}
};

template<typename VertexStruct>
class DynamicMesh {
	vector<VertexStruct> vertices;
	vector<ivec3> indices;
public:
	DynamicMesh();
	virtual ~DynamicMesh();
	DynamicMesh(const vector<VertexStruct>& vertices, const vector<ivec3>& indices);

	void addVertex(const VertexStruct& vertex);
	void addTriangle(int v0, int v1, int v2);
	void addTriangle(const VertexStruct& a, const VertexStruct& b, const VertexStruct& c);
	void addTriangles(const Triangle<VertexStruct>& triangle);
	void reserveVertices(size_t count);
	void reserveIndices(size_t count);
	VertexStruct getVertex(int index) const;
	Triangle<VertexStruct> getTriangle(int index) const;
	ivec3 getTriangleIndices(int index) const;
	void merge(const DynamicMesh& other);
	void transformVertices(const HOM(const VertexStruct&, VertexStruct)& f);
	const vector<ivec3>& getIndexArray() const;
	const void* getVertexData() const;
	size_t size() const;
};

using DynamicMesh3D = DynamicMesh<Vertex3D>;
using Triangle3D = Triangle<Vertex3D>;
using DynamicMesh2D = DynamicMesh<Vertex2D>;
using Triangle2D = Triangle<Vertex2D>;



// -------------------------




template <typename VertexStruct>
DynamicMesh<VertexStruct>::DynamicMesh(const vector<VertexStruct>& vertices, const vector<ivec3>& indices): vertices(vertices), indices(indices) {}

template <typename VertexStruct>
void DynamicMesh<VertexStruct>::addVertex(const VertexStruct& vertex) {
	vertices.push_back(vertex);
}

template <typename VertexStruct>
void DynamicMesh<VertexStruct>::addTriangle(int v0, int v1, int v2) {
	indices.push_back(ivec3(v0, v1, v2));
}

template <typename VertexStruct>
void DynamicMesh<VertexStruct>::addTriangle(const VertexStruct& a, const VertexStruct& b, const VertexStruct& c) {
	addVertex(a);
	int idx_a = vertices.size() - 1;
	addVertex(b);
	int idx_b = vertices.size() - 1;
	addVertex(c);
	int idx_c = vertices.size() - 1;

	addTriangle(idx_a, idx_b, idx_c);
}

template <typename VertexStruct>
void DynamicMesh<VertexStruct>::addTriangles(const Triangle<VertexStruct>& triangle) {
	addTriangle(triangle.a, triangle.b, triangle.c);
}

template <typename VertexStruct>
void DynamicMesh<VertexStruct>::reserveVertices(size_t count) {
	vertices.reserve(count);
}

template <typename VertexStruct>
void DynamicMesh<VertexStruct>::reserveIndices(size_t count) {
	indices.reserve(count);
}

template <typename VertexStruct>
VertexStruct DynamicMesh<VertexStruct>::getVertex(int index) const {
	return vertices.at(index);
}

template <typename VertexStruct>
Triangle<VertexStruct> DynamicMesh<VertexStruct>::getTriangle(int index) const {
	ivec3 idx = indices.at(index);
	return Triangle<VertexStruct>{vertices[idx.x], vertices[idx.y], vertices[idx.z]};
}

template <typename VertexStruct>
ivec3 DynamicMesh<VertexStruct>::getTriangleIndices(int index) const {
	return indices.at(index);
}

template <typename VertexStruct>
void DynamicMesh<VertexStruct>::merge(const DynamicMesh& other) {
	int vertexOffset = vertices.size();
	vertices.insert(vertices.end(), other.vertices.begin(), other.vertices.end());
	for (const auto& idx : other.indices) {
		indices.emplace_back(idx.x + vertexOffset, idx.y + vertexOffset, idx.z + vertexOffset);
	}
}

template <typename VertexStruct>
void DynamicMesh<VertexStruct>::transformVertices(const std::function<VertexStruct(const VertexStruct&)>& f) {
	for (auto& v : vertices)
		v = f(v);
}

template <typename VertexStruct>
const vector<ivec3>& DynamicMesh<VertexStruct>::getIndexArray() const {
	return indices;
}

template <typename VertexStruct>
const void* DynamicMesh<VertexStruct>::getVertexData() const {
	return vertices.data();
}

template <typename VertexStruct>
size_t DynamicMesh<VertexStruct>::size() const { return vertices.size() * sizeof(VertexStruct); }
