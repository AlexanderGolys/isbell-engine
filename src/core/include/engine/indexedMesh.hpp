#pragma once
#include "exceptions.hpp"
#include "shaderTypes.hpp"

template<typename T>
concept attribute_struct = requires(T a) {
	{ T::layout() } -> std::same_as<VertexBufferLayout>;
};

template<attribute_struct vertex>
struct Triangle {
	vertex v0, v1, v2;
	Triangle(vertex v0, vertex v1, vertex v2) : v0(v0), v1(v1), v2(v2) {}
};

class GeometricData : public DirtyFlag {
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
	IndexedMesh(const vector<vertex>& vertices, const vector<ivec3>& indices);
	explicit IndexedMesh(const vector<Triangle<vertex>>& triangles);

	IndexedMesh(const IndexedMesh&) = delete;
	IndexedMesh& operator=(const IndexedMesh&) = delete;

	vector<vertex>::const_iterator begin();
	vector<vertex>::const_iterator end();
	vector<vertex>::iterator begin_mut();
	vector<vertex>::iterator end_mut();

	void setVertices(const vector<vertex>& verts);
	void setIndices(const vector<ivec3>& inds);
	void addVertex(const vertex& v);

	void addTriangleIndices(const ivec3& ind);
	void reserveVertices(array_len n);
	void reserveIndices(array_len n);

	template<typename... Args>
	void emplaceVertex(Args&&... args);

	void emplaceFace(int a, int b, int c);

	array_len numberOfVertices() const final;
	array_len numberOfTriangles() const final;
	byte_size vertexBufferSize() const final;
	byte_size indexBufferSize() const final;
	byte_size totalSize() const final;
	raw_data_ptr vertexBufferData() const final;
	raw_data_ptr indexBufferData() const final;
	VertexBufferLayout layout() const final;

	const vertex& getVertex(array_index i) const;
	void setVertex(array_index i, const vertex& v);
	const ivec3& getTriangleIndices(array_index i) const;
	Triangle<vertex> getTriangle(array_index i) const;

	void transformVertices(const HOM(const vertex&, vertex)& f);
};





// Implementation


template <attribute_struct vertex>
IndexedMesh<vertex>::IndexedMesh(const vector<vertex>& vertices, const vector<ivec3>& indices)
: vertices(vertices), indices(indices) {}

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
typename vector<vertex>::const_iterator IndexedMesh<vertex>::begin() { return vertices.begin(); }

template <attribute_struct vertex>
typename vector<vertex>::const_iterator IndexedMesh<vertex>::end() { return vertices.end(); }

template <attribute_struct vertex>
typename vector<vertex>::iterator IndexedMesh<vertex>::begin_mut() { return vertices.begin(); }

template <attribute_struct vertex>
typename vector<vertex>::iterator IndexedMesh<vertex>::end_mut() { return vertices.end(); }

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
void IndexedMesh<vertex>::reserveVertices(array_len n) {
	vertices.reserve(n);
}

template <attribute_struct vertex>
void IndexedMesh<vertex>::reserveIndices(array_len n) {
	indices.reserve(n);
}

template <attribute_struct vertex>
template <typename ... Args>
void IndexedMesh<vertex>::emplaceVertex(Args&&... args) {
	vertices.emplace_back(std::forward<Args>(args)...);
	markDirty();
}

template <attribute_struct vertex>
void IndexedMesh<vertex>::emplaceFace(int a, int b, int c) {
	indices.emplace_back(a, b, c);
	markDirty();
}

template <attribute_struct vertex>
array_len IndexedMesh<vertex>::numberOfVertices() const {
	return vertices.size();
}

template <attribute_struct vertex>
array_len IndexedMesh<vertex>::numberOfTriangles() const {
	return indices.size();
}

template <attribute_struct vertex>
byte_size IndexedMesh<vertex>::vertexBufferSize() const {
	return vertices.size() * sizeof(vertex);
}

template <attribute_struct vertex>
byte_size IndexedMesh<vertex>::indexBufferSize() const {
	return indices.size() * sizeof(ivec3);
}

template <attribute_struct vertex>
byte_size IndexedMesh<vertex>::totalSize() const {
	return vertexBufferSize() + indexBufferSize();
}

template <attribute_struct vertex>
raw_data_ptr IndexedMesh<vertex>::vertexBufferData() const {
	return vertices[0].data();
}

template <attribute_struct vertex>
raw_data_ptr IndexedMesh<vertex>::indexBufferData() const {
	return &indices[0];
}

template <attribute_struct vertex>
VertexBufferLayout IndexedMesh<vertex>::layout() const {
	return vertex::layout();
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

template <attribute_struct vertex>
void IndexedMesh<vertex>::transformVertices(const std::function<vertex(const vertex&)>& f) {
	for (auto& v : vertices)
		v = f(v);
	markDirty();
}
