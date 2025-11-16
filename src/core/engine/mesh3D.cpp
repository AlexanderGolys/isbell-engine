#include "mesh3D.hpp"

VertexBufferLayout Vertex3D::layout() {
	return VertexBufferLayout({GLSLPrimitive::VEC3, GLSLPrimitive::VEC3, GLSLPrimitive::VEC2, GLSLPrimitive::VEC4 });
}

void Object3D::translate(const vec3& t) {
	transform = SE3(t) * transform;
}

void Object3D::apply(const SE3& T) {
	transform = T * transform;
}

Mesh3D::Mesh3D(const SmoothParametricSurface& surf, int uRes, int vRes) {
	reserveVertices(uRes * vRes);
	reserveIndices(2 * uRes * vRes);

	for (int i = 0; i < uRes; i++)
		for (int j = 0; j < vRes; j++) {
			float t_ = 1.f * i / (uRes - 1);
			if (surf.get_periodicT())
				t_ = 1.f * i / (1.f * uRes);

			float u_ = 1.f * j / (vRes - 1);
			if (surf.get_periodicU())
				u_ = 1.f * j / (1.f * vRes);

			float t = lerp(surf.t0(), surf.t1(), t_);
			float u = lerp(surf.u0(), surf.u1(), u_);

			emplaceVertex(surf(t, u), surf.normal(t, u), vec2(t_, u_), vec4(t, u, 0, 1));
		}

	ivec2 size = ivec2(uRes, vRes);

	for (int i = 1; i < uRes; i++)
		for (int j = 1; j < vRes; j++) {
			int i0 = flattened2DVectorIndex(i, j, size);
			int i1 = flattened2DVectorIndex(i - 1, j, size);
			int i2 = flattened2DVectorIndex(i - 1, j - 1, size);
			int i3 = flattened2DVectorIndex(i, j - 1, size);

			emplaceFace(i0, i1, i2);
			emplaceFace(i0, i3, i2);
		}


	if (surf.get_periodicT())
		for (int j = 1; j < vRes; j++) {
			int i0 = flattened2DVectorIndex(0, j, size);
			int i1 = flattened2DVectorIndex(-1, j, size);
			int i2 = flattened2DVectorIndex(-1, j - 1, size);
			int i3 = flattened2DVectorIndex(0, j - 1, size);
			emplaceFace(i0, i1, i2);
			emplaceFace(i0, i3, i2);
		}
	//
	if (surf.get_periodicU())
		for (int i = 1; i < uRes; i++) {
			int i0 = flattened2DVectorIndex(i, 0, size);
			int i1 = flattened2DVectorIndex(i - 1, 0, size);
			int i2 = flattened2DVectorIndex(i - 1, -1, size);
			int i3 = flattened2DVectorIndex(i, -1, size);
			emplaceFace(i0, i1, i2);
			emplaceFace(i0, i3, i2);
		}

	if (surf.get_periodicT() && surf.get_periodicU()) {
		int i0 = flattened2DVectorIndex(0, 0, size);
		int i1 = flattened2DVectorIndex(-1, 0, size);
		int i2 = flattened2DVectorIndex(-1, -1, size);
		int i3 = flattened2DVectorIndex(0, -1, size);
		emplaceFace(i0, i1, i2);
		emplaceFace(i0, i3, i2);
	}
}


raw_data_ptr Vertex3D::data() const {
	return &position[0];
}

Vertex3D::Vertex3D(vec3 position, vec3 normal, vec2 uv, vec4 color)
: position(position), normal(normal), uv(uv), color(color) {}
