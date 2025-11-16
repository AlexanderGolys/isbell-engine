#include "logging.hpp"
#include "SDFObjects2.hpp"



int main() {
	logging::Logger::init();
	GLSLValidType boxParams = GLSLValidType(GLSLParameterType(GLSLPrimitiveType::GLSL_VEC3));

	GLSLValidType mergeParams = GLSLValidType(GLSLParameterType(GLSLPrimitiveType::GLSL_FLOAT));

	GLSLValidType torusParams = GLSLValidType(make_shared<GLSLStruct>("TorusParams", std::map<string, GLSLValidType>{
		{"r", GLSLParameterType(GLSLPrimitiveType::GLSL_FLOAT)},
		{"R", GLSLParameterType(GLSLPrimitiveType::GLSL_FLOAT)}
	}));


	SDFFunctionParameterised torus_center = SDFFunctionParameterised("sdf_torus_centered", torusParams, "", "length(vec2(length(x.xz) - p.r, x.y)) - p.R;");
	SDFFunctionParameterised torus = torus_center.motionWrapper("sdf_torus");
	SDFFunctionParameterised box_centered = SDFFunctionParameterised("sdf_box_centered", boxParams, "", "vec3 q = abs(x) - p; length(max(q, 0.0)) + min(max(q.x, max(q.y, q.z)), 0.0);");
	SDFFunctionParameterised box = box_centered.motionWrapper("sdf_box");
	// GLSLSmoothBinaryOperator smoothMin = GLSLSmoothBinaryOperator("smooth_min", mergeParams, "float h = clamp(0.5 + 0.5*(y - z)/p.a, 0.0, 1.0);", "mix(y, z, h) - p.a*h*(1.0 - h);");

	LOG_PURE(box_centered.generateCode());
	LOG_PURE(box.generateCode());
	LOG_PURE(torus_center.getParameterType().declarationCode());
	LOG_PURE(torus_center.generateCode());
	LOG_PURE(torus.getParameterType().declarationCode());
	LOG_PURE(torus.generateCode());

	// LOG_PURE(smoothMin(torus, box, "sdf").getParameterType().declarationCode());
	// LOG_PURE(smoothMin(torus, box, "sdf").generateCode());
}
