#pragma once
#include "exceptions.hpp"
#include "shaderTypes.hpp"

template<typename T>
concept attribute_struct = requires(T a) {
	{ a.data() } -> std::same_as<raw_data_ptr>;
	{ T::layout() } -> std::same_as<VertexBufferLayout>;
};

template<attribute_struct vertex>
struct Triangle {
	vertex v0, v1, v2;
	Triangle(vertex v0, vertex v1, vertex v2) : v0(v0), v1(v1), v2(v2) {}
};

class GeometricData {
	DIRTY_FLAG;
public:
	virtual ~GeometricData() = default;
	virtual array_len numberOfVertices() const = 0;
	virtual array_len numberOfTriangles() const = 0;
	virtual byte_size vertexBufferSize() const = 0;
	virtual byte_size indexBufferSize() const = 0;
	virtual byte_size totalSize() const = 0;
	virtual raw_data_ptr vertexBufferData() const = 0;
	virtual raw_data_ptr indexBufferData() const = 0;
	virtual VertexBufferLayout layout() const = 0;

};

template<attribute_struct vertex>
class IndexedMesh : public GeometricData {
	vector<vertex> vertices;
	vector<ivec3> indices;

public:
	IndexedMesh() = default;
	IndexedMesh(const vector<vertex>& vertices, const vector<ivec3>& indices) : vertices(vertices), indices(indices) {}
	explicit IndexedMesh(const vector<Triangle<vertex>>& triangles);

	vector<vertex>::const_iterator begin() { return vertices.begin(); }
	vector<vertex>::const_iterator end() { return vertices.end(); }

	void setVertices(const vector<vertex>& verts);
	void setIndices(const vector<ivec3>& inds);
	void addVertex(const vertex& v);

	void addTriangleIndices(const ivec3& ind);
	void reserveVertices(array_len n) { vertices.reserve(n); }
	void reserveIndices(array_len n) { indices.reserve(n); }

	template<typename... Args>
	void emplaceVertex(Args&&... args) { vertices.emplace_back(std::forward<Args>(args)...); }
	void emplaceFace(int a, int b, int c) { indices.emplace_back(a, b, c); }

	array_len numberOfVertices() const override { return vertices.size(); }
	array_len numberOfTriangles() const override { return indices.size(); }
	byte_size vertexBufferSize() const override { return vertices.size() * sizeof(vertex);  }
	byte_size indexBufferSize() const override { return indices.size() * sizeof(ivec3);  }
	byte_size totalSize() const override { return vertexBufferSize() + indexBufferSize(); }
	raw_data_ptr vertexBufferData() const override { return vertices[0].data();  }
	raw_data_ptr indexBufferData() const override { return &indices[0]; }
	VertexBufferLayout layout() const override { return vertex::layout(); }

	const vertex& getVertex(array_index i) const;
	void setVertex(array_index i, const vertex& v);
	const ivec3& getTriangleIndices(array_index i) const;
	Triangle<vertex> getTriangle(array_index i) const;
};





// Implementation






template <attribute_struct vertex>
IndexedMesh<vertex>::IndexedMesh(const vector<Triangle<vertex>>& triangles) {
	vertices.reserve(triangles.size() * 3);
	indices.reserve(triangles.size());
	for (const auto& tri : triangles) {
		array_len baseIndex = vertices.size();
		vertices.push_back(tri.v0);
		vertices.push_back(tri.v1);
		vertices.push_back(tri.v2);
		indices.push_back(ivec3(baseIndex, baseIndex + 1, baseIndex + 2));
	}
}

template <attribute_struct vertex>
void IndexedMesh<vertex>::setVertices(const vector<vertex>& verts) {
	vertices = verts;
	markDirty();
}

template <attribute_struct vertex>
void IndexedMesh<vertex>::setIndices(const vector<ivec3>& inds) {
	indices = inds;
	markDirty();
}

template <attribute_struct vertex>
void IndexedMesh<vertex>::addVertex(const vertex& v) {
	vertices.push_back(v);
	markDirty();
}

template <attribute_struct vertex>
void IndexedMesh<vertex>::addTriangleIndices(const ivec3& ind) {
	indices.push_back(ind);
	markDirty();
}

template <attribute_struct vertex>
const vertex& IndexedMesh<vertex>::getVertex(array_index i) const {
	CHECK_OUT_OF_BOUNDS(i, vertices.size());
	return vertices[i];
}

template <attribute_struct vertex>
void IndexedMesh<vertex>::setVertex(array_index i, const vertex& v) {
	CHECK_OUT_OF_BOUNDS(i, vertices.size());
	vertices[i] = v;
	markDirty();
}

template <attribute_struct vertex>
const ivec3& IndexedMesh<vertex>::getTriangleIndices(array_index i) const {
	CHECK_OUT_OF_BOUNDS(i, indices.size());
	return indices[i];
}

template <attribute_struct vertex>
Triangle<vertex> IndexedMesh<vertex>::getTriangle(array_index i) const {
	CHECK_OUT_OF_BOUNDS(i, indices.size());
	const ivec3& inds = indices[i];
	return Triangle<vertex>(vertices[inds.x], vertices[inds.y], vertices[inds.z]);
}
