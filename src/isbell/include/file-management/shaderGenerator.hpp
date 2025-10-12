#pragma once
#include "filesUtils.hpp"
#include "macroParsing.hpp"
#include "SDFObjects2.hpp"



class ShaderFunctionModule : public CodeFileDescriptor {
public:
	string name;
	vector<string> functionNames;
	vector<string> structNames;

	using CodeFileDescriptor::getPath, CodeFileDescriptor::exists;
	ShaderFunctionModule(const Path &path, const string &name);

	void updateModuleContent();
	bool containsFunction(const string &functionName);
	bool containsStruct(const string &structName);
};

class ShaderCodeGenerator {
	CodeFileDescriptor baseFile, targetFile;
	string structCode, functionCode, uniformCode, moduleCode, macroKey;
	vector<string> functionNames, structNames, uniformNames;

public:
	ShaderCodeGenerator(const Path &baseFile, const Path &targetFile, const string &macro);

	void addModule(const ShaderFunctionModule &module);
	void addStruct(const string &code, const string &structName);
	void addFunction(const string &code, const string &functionName);
	void addUniform(const string &code, const string &uniformName);
	void generate() const;
};
