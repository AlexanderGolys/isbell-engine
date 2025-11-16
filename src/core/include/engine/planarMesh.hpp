#pragma once
#include "indexedMesh.hpp"
#include "planarGeometry.hpp"


struct PlanarVertex {
	vec2 position;
	vec2 uv;

	static VertexBufferLayout layout();
	raw_data_ptr data() const;
	PlanarVertex(vec2 position, vec2 uv) : position(position), uv(uv) {}
};

class BasicPlanarMesh : public IndexedMesh<PlanarVertex> {
	PROPERTY(SE2, transform);
public:
	using IndexedMesh::IndexedMesh;
	BasicPlanarMesh(const PlanarSurface& surface, int u_res, int v_res);

	void translate(const vec2& t);
	void rotate(float angle);
	void apply(const SE2& T);
	BasicPlanarMesh& operator*=(const SE2& T);
};
