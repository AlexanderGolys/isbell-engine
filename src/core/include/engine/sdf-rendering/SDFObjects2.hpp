#pragma once
#include "func.hpp"

struct GLSLPrimitiveType {
	enum BuiltInType {
		GLSL_VOID,
		GLSL_FLOAT,
		GLSL_INT,
		GLSL_BOOL,
		GLSL_VEC2,
		GLSL_VEC3,
		GLSL_VEC4,
		GLSL_MAT2,
		GLSL_MAT3,
		GLSL_MAT4,
		GLSL_IVEC2,
		GLSL_IVEC3,
		GLSL_IVEC4,
		SE3_STRUCT
	} type;

	string typeName() const;
	GLSLPrimitiveType(BuiltInType t);
	bool operator==(const GLSLPrimitiveType &other) const;
};

struct GLSLParameterType {
	GLSLPrimitiveType type;
	int arraySize;

	GLSLParameterType(GLSLPrimitiveType type, int arraySize = 1);
	GLSLParameterType(GLSLPrimitiveType::BuiltInType type, int arraySize = 1);

	string typeName() const;
	bool operator==(const GLSLParameterType &other) const;
};


string glslPrimitiveTypeName(GLSLPrimitiveType t);
string glslParameterTypeName(GLSLParameterType t);

// Forward declaration
struct GLSLStruct;

struct GLSLValidType
{
    variant<shared_ptr<GLSLStruct>, GLSLParameterType> type;

	GLSLValidType(shared_ptr<GLSLStruct> t);
	GLSLValidType(const GLSLParameterType &t);
	GLSLValidType(const GLSLPrimitiveType &t);
	GLSLValidType(GLSLPrimitiveType::BuiltInType t);

	string typeName() const;
	bool isStruct() const;
	bool operator==(const GLSLValidType &other) const;
	bool operator==(const GLSLStruct &other) const;
	bool operator==(const GLSLParameterType &other) const;
	bool operator==(const GLSLPrimitiveType &other) const;
	bool operator==(const GLSLPrimitiveType::BuiltInType &other) const;

	string declarationCode() const;
};





struct GLSLStruct
{
    string name;
    dict(string, GLSLValidType) members;
	GLSLStruct(const string &name, const dict(string, GLSLValidType) &members);

    string declarationCode() const;
    string typeName() const;
	bool operator==(const GLSLStruct &other) const;
	GLSLValidType operator[](const string &memberName) const;
};


class SDFStateStruct : public GLSLStruct {
public:
	SDFStateStruct(const GLSLValidType &base);
};

class GLSLFunctionSignature {
	string name;
	GLSLValidType returnType;
	vector<pair<GLSLValidType, string>> arguments;
public:
	GLSLFunctionSignature(const string &name, const GLSLValidType &returnType, const vector<pair<GLSLValidType, string>> &arguments);

	GLSLValidType getReturnType() const;
	string getName() const;
	vector<pair<GLSLValidType, string>> getArguments() const;
	vector<GLSLValidType> getArgumentTypes() const;
	vector<string> getArgumentNames() const;
	string generateDeclarationCode() const;

	string generateCall(const vector<string> &argNames) const;
};

class GLSLFunction {
public:
	GLSLFunctionSignature signature;
	string body;
	string returnLine;

	GLSLFunction(const GLSLFunctionSignature &signature, const string &body, const string &returnLine);
	GLSLFunction(const string &name, const GLSLValidType &returnType, const vector<pair<GLSLValidType, string>> &arguments, const string &body, const string &returnLine);

	string generateCode() const;
	string getName() const;
	GLSLValidType getReturnType() const;
	vector<pair<GLSLValidType, string>> getArguments() const;
	vector<GLSLValidType> getArgumentTypes() const;
	vector<string> getArgumentNames() const;
};


class SDFFunction : public GLSLFunction {
public:
	SDFFunction(const string &name, const string &body, const string &returnLine) ;
};

class SDFFunctionParameterised : public GLSLFunction {
	GLSLValidType _parameterType;
public:
	SDFFunctionParameterised(const string &name, const GLSLValidType &parameters, const string &body, const string &returnLine);

	SDFFunctionParameterised motionWrapper(const string &name) const;
	SDFFunction addIsometriesAndSubstitute(const string &name, const string &stateStructName) const;
	SDFFunction substituteParameters(const string &name, const string &paramStructName) const;
	GLSLValidType getParameterType() const;
};


class SDFObj {
	SDFFunction sdfPure;
	SE3 transform;

};







// class GLSLSmoothBinaryOperator : public GLSLFunction {
// 	GLSLValidType _parameterType;
// public:
// 	GLSLSmoothBinaryOperator(const string &name, const GLSLValidType &parameters, const string &body, const string &returnLine);
// };
