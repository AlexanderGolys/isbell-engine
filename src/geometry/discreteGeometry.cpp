 #include "discreteGeometry.hpp"

using namespace glm;
 glm::mat3 TriangulatedManifold::faceVertices(int i) { return mat3(triangles[id][i].getVertex(0).getPosition(),
																   triangles[id][i].getVertex(1).getPosition(), triangles[id][i].getVertex(2).getPosition()); }

 mat2x3 TriangulatedManifold::orthoFaceTangents(int i) const {
	 mat3 frame = tangentNormalFrameOfFace(i);
	 return mat2x3(frame[0], frame[1]);
 }

 Discrete1Form DiscreteRealFunction::exteriorDerivative() const {
	 vector<float> values = {};
	 values.reserve(domain->numEdges());
	 for (int i = 0; i < domain->numEdges(); i++)
		 values.emplace_back((values[domain->edgeVerticesIndices(i).x] - values[domain->edgeVerticesIndices(i).y])/edgeLength(i));
	 return Discrete1Form(values, domain);
 }
