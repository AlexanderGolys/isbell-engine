#include "fluidsSpecific.hpp"

#include <utility>

using namespace glm;
using std::vector, std::string, std::shared_ptr, std::unique_ptr, std::pair, std::make_unique, std::make_shared, std::sin, std::cos, std::exp, std::log;


//HexVolumetricMeshWithBoundary volumetricRing(vec3 center, vec3 normal, float radiusSmall, float radiusBig, float height, int radial_res, int vertical_res, int extrude_res) {
//
//}

HexVolumetricMeshWithBoundary volumetricCylinder(vec3 center, vec3 normal, float radius, float height, int radial_res, int vertical_res, int extrude_res) {
	throw std::logic_error("Not implemented");
}

HexVolumetricMeshWithBoundary primitiveHex(std::array<vec3, 8> vertices) { throw std::logic_error("Not implemented"); }

HexVolumetricMeshWithBoundary evenGrid(int n, int m, int l, float lenx, float leny, float lenz) {
	vector<Vertex> vertices = {};
	vector<ivec4> faces = {};
	vector<ivec6> hexahedra = {};
	std::vector<int> bdVertices = {};

	for (int i = 0; i <= n; i++)
		for (int j = 0; j <= m; j++)
			for (int k = 0; k <= l; k++) {
				vec3 normal = vec3(1, 0, 0);
				if (i == n) normal = vec3(-1, 0, 0);
				if (j == 0) normal = vec3(0, 1, 0);
				if (j == m) normal = vec3(0, -1, 0);
				if (k == 0) normal = vec3(0, 0, 1);
				if (k == l) normal = vec3(0, 0, -1);
				vertices.emplace_back(vec3(lenx * i / n, leny * j / m, lenz * k / l), normal);
				if (i == 0 || i == n || j == 0 || j == m || k == 0 || k == l)
					bdVertices.push_back(vertices.size() - 1);
			}
	ivec3 dims = ivec3(n+1, m+1, l+1);
	ivec4 dimsFace = ivec4(n, m, l, 3);
	ivec3 dimsCells= ivec3(n, m, l);

	for (int i = 0; i <= n; i++)
		for (int j = 0; j < m; j++)
			for (int k = 0; k < l; k++)
				faces.emplace_back(flattened3DVectorIndex(i, j, k, dims),
								   flattened3DVectorIndex(i, j+1, k, dims),
								   flattened3DVectorIndex(i, j+1, k+1, dims),
								   flattened3DVectorIndex(i, j, k+1, dims));

	for (int i = 0; i < n; i++)
		for (int j = 0; j <= m; j++)
			for (int k = 0; k < l; k++)
				faces.emplace_back(flattened3DVectorIndex(i, j, k, dims),
								   flattened3DVectorIndex(i+1, j, k, dims),
								   flattened3DVectorIndex(i+1, j, k+1, dims),
								   flattened3DVectorIndex(i, j, k+1, dims));
	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			for (int k = 0; k <= l; k++)
				faces.emplace_back(flattened3DVectorIndex(i, j, k, dims),
								   flattened3DVectorIndex(i+1, j, k, dims),
								   flattened3DVectorIndex(i+1, j+1, k, dims),
								   flattened3DVectorIndex(i, j+1, k, dims));

	vector<ivec3> faceDims = {ivec3(n+1, m, l), ivec3(n, m+1, l), ivec3(n, m, l+1)};
	ivec2 shifts = ivec2((n+1)*m*l, n*(m+1)*l+(n+1)*m*l);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			for (int k = 0; k < l; k++) {
				hexahedra.emplace_back(
				flattened3DVectorIndex(i, j, k, faceDims[0]),
				flattened3DVectorIndex(i+1, j, k, faceDims[0]),
				flattened3DVectorIndex(i, j, k, faceDims[1]) + shifts[0],
				flattened3DVectorIndex(i, j+1, k, faceDims[1]) + shifts[0],
				flattened3DVectorIndex(i, j, k, faceDims[2]) + shifts[1],
				flattened3DVectorIndex(i, j, k+1, faceDims[2]) + shifts[1]
				);
			}

	return HexVolumetricMeshWithBoundary(hexahedra, faces, vertices, bdVertices);

}

CellMesh gridCellMesh(glm::ivec3 dims, vec3 len, const std::map<std::string, HOM(glm::ivec3, float)> &attributes) {
	std::vector<std::string> attributeNames = keys(attributes);
	vector<Cell> cells = {};
	for (int i=0; i<dims[0]; i++)
		for (int j=0; j<dims[1]; j++)
			for (int k=0; k<dims[2]; k++) {
				std::map<FACE, int> neighbours = {
					{DOWN, k == 0? -1 : flattened3DVectorIndex(i, j, k-1, dims)},
					{UP, k == dims[2]-1? -1 : flattened3DVectorIndex(i, j, k+1, dims)},
					{LEFT, i == 0? -1 : flattened3DVectorIndex(i-1, j, k, dims)},
					{RIGHT, i == dims[0]-1? -1 : flattened3DVectorIndex(i+1, j, k, dims)},
					{FRONT, j == 0? -1 : flattened3DVectorIndex(i, j-1, k, dims)},
					{BACK, j == dims[1]-1? -1 : flattened3DVectorIndex(i, j+1, k, dims)}
				};
				std::map<std::array<FACE, 3>, vec3> corners = {
					{{DOWN, LEFT, FRONT}, vec3(i*len.x/dims.x, j*len.y/dims.y, k*len.z/dims.z)},
					{{DOWN, RIGHT, FRONT}, vec3((i+1)*len.x/dims.x, j*len.y/dims.y, k*len.z/dims.z)},
					{{DOWN, LEFT, BACK}, vec3(i*len.x/dims.x, (j+1)*len.y/dims.y, k*len.z/dims.z)},
					{{DOWN, RIGHT, BACK}, vec3((i+1)*len.x/dims.x, (j+1)*len.y/dims.y, k*len.z/dims.z)},
					{{UP, LEFT, FRONT}, vec3(i*len.x/dims.x, j*len.y/dims.y, (k+1)*len.z/dims.z)},
					{{UP, RIGHT, FRONT}, vec3((i+1)*len.x/dims.x, j*len.y/dims.y, (k+1)*len.z/dims.z)},
					{{UP, LEFT, BACK}, vec3(i*len.x/dims.x, (j+1)*len.y/dims.y, (k+1)*len.z/dims.z)},
					{{UP, RIGHT, BACK}, vec3((i+1)*len.x/dims.x, (j+1)*len.y/dims.y, (k+1)*len.z/dims.z)}
				};
				std::vector<FACE> bdFaces = {};
				for (int f = 0; f < 6; f++)
					if (neighbours[static_cast<FACE>(f)] == -1)
						bdFaces.push_back(static_cast<FACE>(f));

				std::map<std::string, float> cellAttributes = {};
				for (const std::string & name : attributeNames)
					cellAttributes[name] = attributes.at(name)(glm::ivec3(i, j, k));

				cells.emplace_back(cells.size(), neighbours, corners, bdFaces, cellAttributes);
			}

	return CellMesh(cells);
}

float HeatFlow::changeOverTime(const Cell &c) {
	float a = (temperatureFluxPerFace(c, LEFT)*c.getFaceArea(LEFT) + temperatureFluxPerFace(c, RIGHT)*c.getFaceArea(RIGHT))	+
		   (temperatureFluxPerFace(c, FRONT)*c.getFaceArea(FRONT) + temperatureFluxPerFace(c, BACK)*c.getFaceArea(BACK))	+
		   (temperatureFluxPerFace(c, DOWN)*c.getFaceArea(DOWN) + temperatureFluxPerFace(c, UP)*c.getFaceArea(UP));
	return a;
}

VectorFieldR3Dynamic CouetteFlow(float h, float v0) {
	return VectorFieldR3Dynamic([h, v0](float t, vec3 x) {
		return vec3(v0*x.y/h, 0, 0);
	});
}


CellMesh flatGridWithHeight(glm::ivec2 dims, vec2 len, std::function<float(glm::ivec2)> h, const std::map<std::string, std::function<float(glm::ivec2)>> &attributes) {
	std::vector<std::string> attributeNames = keys(attributes);
	vector<Cell> cells = {};
	for (int i=0; i<dims[0]; i++)
		for (int j=0; j<dims[1]; j++){
				std::map<FACE, int> neighbours = {
					{LEFT, i == 0? -1 : flattened2DVectorIndex(i-1, j, dims)},
					{RIGHT, i == dims[0]-1? -1 : flattened2DVectorIndex(i+1, j, dims)},
					{FRONT, j == 0? -1 : flattened2DVectorIndex(i, j-1, dims)},
					{BACK, j == dims[1]-1? -1 : flattened2DVectorIndex(i, j+1, dims)},
						{DOWN, -1}, {UP, -1}
				};
				std::map<std::array<FACE, 3>, vec3> corners = {
					{{DOWN, LEFT, FRONT}, vec3(i*len.x/dims.x, j*len.y/dims.y, 0)},
					{{DOWN, RIGHT, FRONT}, vec3((i+1)*len.x/dims.x, j*len.y/dims.y, 0)},
					{{DOWN, LEFT, BACK}, vec3(i*len.x/dims.x, (j+1)*len.y/dims.y, 0)},
					{{DOWN, RIGHT, BACK}, vec3((i+1)*len.x/dims.x, (j+1)*len.y/dims.y, 0)},
					{{UP, LEFT, FRONT}, vec3(i*len.x/dims.x, j*len.y/dims.y, h(ivec2(i, j)))},
					{{UP, RIGHT, FRONT}, vec3((i+1)*len.x/dims.x, j*len.y/dims.y, h(ivec2(i+1, j)))},
					{{UP, LEFT, BACK}, vec3(i*len.x/dims.x, (j+1)*len.y/dims.y, h(ivec2(i, j+1)))},
					{{UP, RIGHT, BACK}, vec3((i+1)*len.x/dims.x, (j+1)*len.y/dims.y, h(ivec2(i+1, j+1)))}
				};
				std::vector<FACE> bdFaces = {};
				for (int f = 0; f < 6; f++)
					if (neighbours[static_cast<FACE>(f)] == -1)
						bdFaces.push_back(static_cast<FACE>(f));

				std::map<std::string, float> cellAttributes = {};
				for (const std::string & name : attributeNames)
					cellAttributes[name] = attributes.at(name)(glm::ivec2(i, j));

				cells.emplace_back(cells.size(), neighbours, corners, bdFaces, cellAttributes);
			}

	return CellMesh(cells);
}


CellMesh flatCustomGrid(glm::ivec2 dims, const std::function<vec3(glm::ivec2)>& position_lower, const std::function<vec3(glm::ivec2)>& position_higher, const std::map<std::string, std::function<float(glm::ivec2)>> &attributes) {
	std::vector<std::string> attributeNames = keys(attributes);
	vector<Cell> cells = {};
	for (int i=0; i<dims[0]; i++)
		for (int j=0; j<dims[1]; j++){
			std::map<FACE, int> neighbours = {
				{LEFT, i == 0? -1 : flattened2DVectorIndex(i-1, j, dims)},
				{RIGHT, i == dims[0]-1? -1 : flattened2DVectorIndex(i+1, j, dims)},
				{FRONT, j == 0? -1 : flattened2DVectorIndex(i, j-1, dims)},
				{BACK, j == dims[1]-1? -1 : flattened2DVectorIndex(i, j+1, dims)},
				{DOWN, -1},
				{UP, -1}
			};
			std::map<std::array<FACE, 3>, vec3> corners = {
				{{DOWN, LEFT, FRONT}, vec3(position_lower(ivec2(i, j)))},
				{{DOWN, RIGHT, FRONT}, vec3(position_lower(ivec2(i+1, j)))},
				{{DOWN, LEFT, BACK}, vec3(position_lower(ivec2(i, j+1)))},
				{{DOWN, RIGHT, BACK}, vec3(position_lower(ivec2(i+1, j+1)))},
				{{UP, LEFT, FRONT}, vec3(position_higher(ivec2(i, j)))},
				{{UP, RIGHT, FRONT}, vec3(position_higher(ivec2(i+1, j)))},
				{{UP, LEFT, BACK}, vec3(position_higher(ivec2(i, j+1)))},
				{{UP, RIGHT, BACK}, vec3(position_higher(ivec2(i+1, j+1)))}
			};
			std::vector<FACE> bdFaces = {};
			for (auto f: keys(neighbours))
				if (neighbours[f] == -1)
					bdFaces.push_back(f);

			std::map<std::string, float> cellAttributes = {};
			for (const std::string & name : attributeNames)
				cellAttributes[name] = attributes.at(name)(glm::ivec2(i, j));

			cells.emplace_back(cells.size(), neighbours, corners, bdFaces, cellAttributes);
		}

	return CellMesh(cells);
}

CellMesh flatCustomSubgrid(glm::ivec2 dims, const HOM(glm::ivec2, vec3) &position_lower, const HOM(glm::ivec2, vec3) &position_higher,
						   const std::vector<glm::ivec2>& removed, const std::map<std::string, std::function<float(glm::ivec2)>> &attributes) {

	CellMesh grid = flatCustomGrid(dims, position_lower, position_higher, attributes);
	std::vector<glm::ivec2> toRemove = removed;
	std::ranges::sort(toRemove, [dims](const glm::ivec2 &a, const glm::ivec2 &b) {
		return flattened2DVectorIndex(a.x, a.y, dims) > flattened2DVectorIndex(b.x, b.y, dims); });
	for (const glm::ivec2 & r : toRemove)
		grid.removeCell(flattened2DVectorIndex(r.x, r.y, dims));
	return grid;
}


CellMesh flatSubgrid(glm::ivec2 dims, vec3 len, const std::vector<glm::ivec2> &removed, const std::map<std::string, std::function<float(glm::ivec2)>> &attributes) {
	return flatCustomSubgrid(
		dims,
		[dims, len](ivec2 i) { return vec3((len.x*i.x)/dims.x, (len.y*i.y)/dims.y, 0); },
		[dims, len](ivec2 i) { return vec3((len.x*i.x)/dims.x, (len.y*i.y)/dims.y, len.z); },
		removed,
		attributes);
}

CellMesh flatCylinder(glm::ivec2 dims, float r, float R, float h, const std::map<std::string, std::function<float(glm::ivec2)>> &attributes) {
	std::vector<std::string> attributeNames = keys(attributes);
	vector<Cell> cells = {};
	for (int i=0; i<dims[0]; i++)
		for (int j=0; j<dims[1]; j++){
				std::map<FACE, int> neighbours = {
					{DOWN, -1},
					{UP, -1},
					{LEFT, i==0? flattened2DVectorIndex(dims[0]-1, j, dims) : flattened2DVectorIndex(i-1, j, dims)},
					{RIGHT, i==dims[0]-1 ? flattened2DVectorIndex(0, j, dims) : flattened2DVectorIndex(i+1%dims[0], j, dims)},
					{FRONT, j == 0? -1 : flattened2DVectorIndex(i, j-1, dims)},
					{BACK, j == dims[1]-1? -1 : flattened2DVectorIndex(i, j+1, dims)}
				};
				float theta = TAU*i/dims[0];
				std::map<std::array<FACE, 3>, vec3> corners = {
					{{DOWN, LEFT, FRONT}, vec3(cos(theta)*r, sin(theta)*r, h*j/dims[1])},
					{{DOWN, RIGHT, FRONT}, vec3(cos(theta+TAU/dims[0])*r, sin(theta+TAU/dims[0])*r, h*j/dims[1])},
					{{DOWN, LEFT, BACK}, vec3(cos(theta)*r, sin(theta)*r, h*(j+1)/dims[1])},
					{{DOWN, RIGHT, BACK}, vec3(cos(theta+TAU/dims[0])*r, sin(theta+TAU/dims[0])*r, h*(j+1)/dims[1])},
					{{UP, LEFT, FRONT}, vec3(cos(theta)*R, sin(theta)*R, h*j/dims[1])},
					{{UP, RIGHT, FRONT}, vec3(cos(theta+TAU/dims[0])*R, sin(theta+TAU/dims[0])*R, h*j/dims[1])},
					{{UP, LEFT, BACK}, vec3(cos(theta)*R, sin(theta)*R, h*(j+1)/dims[1])},
					{{UP, RIGHT, BACK}, vec3(cos(theta+TAU/dims[0])*R, sin(theta+TAU/dims[0])*R, h*(j+1)/dims[1])}
				};
				std::vector<FACE> bdFaces = {};
				for (int f = 0; f < 6; f++)
					if (neighbours[static_cast<FACE>(f)] == -1)
						bdFaces.push_back(static_cast<FACE>(f));

				std::map<std::string, float> cellAttributes = {};
				for (const std::string & name : attributeNames)
					cellAttributes[name] = attributes.at(name)(glm::ivec2(i, j));

				cells.emplace_back(cells.size(), neighbours, corners, bdFaces, cellAttributes);
			}

	return CellMesh(cells);
}
