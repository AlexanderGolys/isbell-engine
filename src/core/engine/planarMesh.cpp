#include "planarMesh.hpp"

VertexBufferLayout PlanarVertex::layout() {
	return VertexBufferLayout({GLSLPrimitive::VEC2, GLSLPrimitive::VEC2 });
}

raw_data_ptr PlanarVertex::data() const {
	return &position[0];
}

BasicPlanarMesh::BasicPlanarMesh(const PlanarSurface& surface, int u_res, int v_res) {
	auto surf = surface;
	surf.normaliseDomain();
	reserveVertices(u_res * v_res);
	reserveIndices(2 * u_res * v_res);
	for (int i_u = 0; i_u < u_res; i_u++) {
		for (int i_v = 0; i_v < v_res; i_v++) {
			float u = 1.f * i_u / (u_res - (int)surf.get_uPeriodic());
			float v = 1.f * i_v / (v_res - (int)surf.get_vPeriodic());
			emplaceVertex(surf(u, v), vec2(u, v));
			if (i_u > 0 && i_v > 0) {
				int i0 = (i_u - 1) * v_res + i_v - 1;
				int i1 = (i_u - 1) * v_res + i_v;
				int i2 = i_u * v_res + i_v - 1;
				int i3 = i_u * v_res + i_v;
				emplaceFace(i0, i1, i2);
				emplaceFace(i2, i1, i3);
			}
			if (surf.get_uPeriodic() && i_u == u_res - 1 && i_v > 0) {
				int i0 = i_v - 1;
				int i1 = i_v;
				int i2 = (u_res - 1) * v_res + i_v - 1;
				int i3 = (u_res - 1) * v_res + i_v;
				emplaceFace(i0, i1, i2);
				emplaceFace(i2, i1, i3);
			}
			if (surf.get_vPeriodic() && i_v == v_res - 1 && i_u > 0) {
				int i0 = (i_u - 1) * v_res;
				int i1 = (i_u - 1) * v_res + v_res - 1;
				int i2 = i_u * v_res;
				int i3 = i_u * v_res + v_res - 1;
				emplaceFace(i0, i1, i2);
				emplaceFace(i2, i1, i3);
			}
			if (surf.get_uPeriodic() && surf.get_vPeriodic() && i_u == u_res - 1 && i_v == v_res - 1) {
				int i0 = 0;
				int i1 = v_res - 1;
				int i2 = i_u * v_res;
				int i3 = i_u * v_res + i_v;
				emplaceFace(i0, i1, i2);
				emplaceFace(i2, i1, i3);
			}
		}
	}
}

void BasicPlanarMesh::translate(const vec2& t) {
	transform = SE2(t) * transform;
}

void BasicPlanarMesh::rotate(float angle) {
	transform = SE2(angle) * transform;
}

void BasicPlanarMesh::apply(const SE2& T) {
	transform = T * transform;
}

BasicPlanarMesh& BasicPlanarMesh::operator*=(const SE2& T) {
	apply(T);
	return *this;
}
