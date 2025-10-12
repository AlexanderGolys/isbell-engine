#include "SDFObjects2.hpp"


string GLSLPrimitiveType::typeName() const {
	switch (type) {
		case GLSL_FLOAT:	return "float";
		case GLSL_INT:  	return "int";
		case GLSL_BOOL: 	return "bool";
		case GLSL_VEC2: 	return "vec2";
		case GLSL_VEC3:		return "vec3";
		case GLSL_VEC4:		return "vec4";
		case GLSL_MAT2:		return "mat2";
		case GLSL_MAT3:		return "mat3";
		case GLSL_MAT4:		return "mat4";
		case GLSL_IVEC2:	return "ivec2";
		case GLSL_IVEC3:	return "ivec3";
		case GLSL_IVEC4:	return "ivec4";
		case GLSL_VOID:		return "void";
		case SE3_STRUCT:	return "SE3";
		default:
			THROW(IllegalVariantError, "Unrecognised GLSLPrimitiveType");
	}
}

GLSLPrimitiveType::GLSLPrimitiveType(BuiltInType t): type(t) {}

GLSLValidType::GLSLValidType(shared_ptr<GLSLStruct> t): type(t) {}

GLSLValidType::GLSLValidType(const GLSLParameterType &t): type(t) {}

GLSLValidType::GLSLValidType(const GLSLPrimitiveType &t): type(GLSLParameterType(t)) {
}

GLSLValidType::GLSLValidType(GLSLPrimitiveType::BuiltInType t): type(GLSLParameterType(t)) {
}

SDFFunctionParameterised::SDFFunctionParameterised(const string &name, const GLSLValidType &parameters, const string &body, const string &returnLine):
GLSLFunction(name, GLSLValidType(GLSLPrimitiveType::GLSL_FLOAT), { {parameters, "p"}, {GLSLValidType(GLSLPrimitiveType::GLSL_VEC3), "x"}}, body, returnLine),
_parameterType(parameters) {}





bool GLSLPrimitiveType::operator==(const GLSLPrimitiveType &other) const {
	return type == other.type;
}


GLSLParameterType::GLSLParameterType(GLSLPrimitiveType type, int arraySize): type(type), arraySize(arraySize) {
	if (arraySize < 1)
		throw IllegalArgumentError("Array size must be at least 1.", __FILE__, __LINE__);
}

GLSLParameterType::GLSLParameterType(GLSLPrimitiveType::BuiltInType type, int arraySize): type(type), arraySize(arraySize) {}

string GLSLParameterType::typeName() const { return glslParameterTypeName(*this); }

bool GLSLParameterType::operator==(const GLSLParameterType &other) const {
	return type.typeName() == other.type.typeName();
}

string glslPrimitiveTypeName(GLSLPrimitiveType t) {
	return t.typeName();
}

string glslParameterTypeName(GLSLParameterType t) {
	if (t.arraySize == 1)
		return glslPrimitiveTypeName(t.type);
	return glslPrimitiveTypeName(t.type) + "[" + std::to_string(t.arraySize) + "]";
}

bool GLSLValidType::operator==(const GLSLValidType &other) const {
	if (isStruct() != other.isStruct())
		return false;
	if (isStruct())
		return *get<shared_ptr<GLSLStruct>>(type) == *get<shared_ptr<GLSLStruct>>(other.type);
	return get<GLSLParameterType>(type) == get<GLSLParameterType>(other.type);
}

bool GLSLValidType::operator==(const GLSLStruct &other) const {
	return isStruct() && other == *get<shared_ptr<GLSLStruct>>(type);
}

bool GLSLValidType::operator==(const GLSLParameterType &other) const {
	return !isStruct() && other == get<GLSLParameterType>(type);
}

bool GLSLValidType::operator==(const GLSLPrimitiveType &other) const {
	return !isStruct() && other == get<GLSLParameterType>(type).type;
}

bool GLSLValidType::operator==(const GLSLPrimitiveType::BuiltInType &other) const {
	return !isStruct() && other == get<GLSLParameterType>(type).type.type;
}

string GLSLValidType::declarationCode() const {
	if (isStruct())
		return get<shared_ptr<GLSLStruct>>(type)->declarationCode();
	return "";
}


GLSLStruct::GLSLStruct(const string &name, const dict(string, GLSLValidType) &members): name(name), members(members) {}

string GLSLStruct::declarationCode() const {
	string s = "struct " + name + " {\n";
	for (const auto& [memberName, memberType] : members)
		s += "    " + memberType.typeName() + " " + memberName + ";\n";
	s += "};\n\n";
	return s;
}

string GLSLStruct::typeName() const {
	return name;
}

bool GLSLStruct::operator==(const GLSLStruct &other) const {
	return declarationCode() == other.declarationCode();
}

GLSLValidType GLSLStruct::operator[](const string &memberName) const {
	return members.at(memberName);
}

string GLSLValidType::typeName() const {
	if (holds_alternative<GLSLParameterType>(type))
		return get<GLSLParameterType>(type).typeName();
	return get<shared_ptr<GLSLStruct>>(type)->typeName();
}

bool GLSLValidType::isStruct() const {
	return not holds_alternative<GLSLParameterType>(type);
}



GLSLFunctionSignature::GLSLFunctionSignature(const string &name, const GLSLValidType &returnType, const vector<pair<GLSLValidType, string>> &arguments): name(name), returnType(returnType), arguments(arguments) {}

GLSLValidType GLSLFunctionSignature::getReturnType() const { return returnType; }

string GLSLFunctionSignature::getName() const { return name; }

vector<pair<GLSLValidType, string>> GLSLFunctionSignature::getArguments() const { return arguments; }

vector<GLSLValidType> GLSLFunctionSignature::getArgumentTypes() const {
	return map<pair<GLSLValidType, string>, GLSLValidType>(arguments, [](const pair<GLSLValidType, string> &p) { return p.first; });
}

vector<string> GLSLFunctionSignature::getArgumentNames() const {
	return map<pair<GLSLValidType, string>, string>(arguments, [](const pair<GLSLValidType, string> &p) { return p.second; });
}

string GLSLFunctionSignature::generateDeclarationCode() const {
	return returnType.typeName() + " " + name + "(" + join(", ", map<pair<GLSLValidType, string>, string>(arguments, [](const pair<GLSLValidType, string> &p) { return p.first.typeName() + " " + p.second; })) + ")";
}

string GLSLFunctionSignature::generateCall(const vector<string> &argNames) const {
	if (argNames.size() != arguments.size())
		THROW(ValueError, "Number of arguments does not match in function call to " + name);
	return name + "(" + join(", ", argNames) + ")";
}

GLSLFunction::GLSLFunction(const GLSLFunctionSignature &signature, const string &body, const string &returnLine): signature(signature), body(body), returnLine(returnLine) {}

GLSLFunction::GLSLFunction(const string &name, const GLSLValidType &returnType, const vector<pair<GLSLValidType, string>> &arguments, const string &body, const string &returnLine): GLSLFunction(GLSLFunctionSignature(name, returnType, arguments), body, returnLine) {}

string GLSLFunction::generateCode() const {
	string code = signature.generateDeclarationCode() + "{\n";
	code += body + "\n";
	code += "\treturn " + returnLine + ";\n";
	code += "}\n";
	return replaceAll(code, "\n\n", "\n");
}

string GLSLFunction::getName() const {
	return signature.getName();
}

GLSLValidType GLSLFunction::getReturnType() const {
	return signature.getReturnType();
}

vector<pair<GLSLValidType, string>> GLSLFunction::getArguments() const {
	return signature.getArguments();
}

vector<GLSLValidType> GLSLFunction::getArgumentTypes() const {
	return signature.getArgumentTypes();
}

vector<string> GLSLFunction::getArgumentNames() const {
	return signature.getArgumentNames();
}

SDFFunctionParameterised SDFFunctionParameterised::motionWrapper(const string &name) const {
	return SDFFunctionParameterised(name, GLSLValidType(make_shared<SDFStateStruct>(_parameterType)), "", this->signature.generateCall({"p.params", "act_SE3(inv_SE3(p.g), x)"}));
}

SDFFunction SDFFunctionParameterised::addIsometriesAndSubstitute(const string &name, const string &stateStructName) const {
	return SDFFunction(name, "", this->signature.generateCall({stateStructName + ".params", "act_SE3(inv_SE3(" + stateStructName + ".g), x)"}));
}

SDFFunction SDFFunctionParameterised::substituteParameters(const string &name, const string &paramStructName) const {
	return SDFFunction(name, "", this->signature.generateCall({paramStructName, "x"}));
}

GLSLValidType SDFFunctionParameterised::getParameterType() const {
	return _parameterType;
}

SDFFunction::SDFFunction(const string &name, const string &body, const string &returnLine)
: GLSLFunction(name, GLSLValidType(GLSLPrimitiveType::GLSL_FLOAT), { {GLSLValidType(GLSLPrimitiveType::GLSL_VEC3), "x"} }, body, returnLine) {}




SDFStateStruct::SDFStateStruct(const GLSLValidType &base):
GLSLStruct("State_" + base.typeName(),
			dict(string, GLSLValidType){
				{"g", GLSLValidType(GLSLPrimitiveType::SE3_STRUCT)},
				{"params", base}
			}) {}
