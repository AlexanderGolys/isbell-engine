#include "shaderGenerator.hpp"





bool startsWithStruct(const string &line) {
	size_t i = 0;
	while (i < line.size() && (line[i] == ' ' || line[i] == '\t'))
		i++;
	if (i + 6 <= line.size() && line.substr(i, 6) == "struct")
		return true;
	return false;
}


bool isPreprocessorLine(const string &line) {
	size_t i = 0;
	while (i < line.size() && (line[i] == ' ' || line[i] == '\t'))
		i++;
	if (i < line.size() && line[i] == '#')
		return true;
	return false;
}


optional<string> extractStructName(const string &line) {
	if (!startsWithStruct(line))
		return UNDEFINED;
	size_t i = 0;
	while (i < line.size() && (line[i] == ' ' || line[i] == '\t'))
		i++;
	i += 6;
	while (i < line.size() && (line[i] == ' ' || line[i] == '\t'))
		i++;
	if (i >= line.size() || !isIdentStart(line[i]))
		return UNDEFINED;
	size_t j = i + 1;
	while (j < line.size() && isIdentChar(line[j]))
		j++;
	string name = line.substr(i, j - i);
	return name;
}


optional<string> extractFunctionName(const string &line) {
	size_t searchFrom = 0;
	while (true) {
		size_t p = line.find('(', searchFrom);
		if (p == string::npos)
			return UNDEFINED;

		int parDepth = 1;
		size_t q = p + 1;
		while (q < line.size() && parDepth > 0) {
			if (line[q] == '(')
				parDepth++;
			else if (line[q] == ')')
				parDepth--;
			q++;
		}
		if (parDepth != 0)
			return UNDEFINED;

		// q points one past ')'
		size_t r = q;
		while (r < line.size() && (line[r] == ' ' || line[r] == '\t'))
			r++;
		// for function declaration/definition, next token should be ';' or end of line (after stripCurlyBlocks)
		if (r < line.size() && line[r] != ';') {
			// not a function declaration, move past this '(' and continue searching
			searchFrom = p + 1;
			continue;
		}

		if (p == 0)
			return UNDEFINED;
		size_t j = p;
		while (j > 0 && (line[j - 1] == ' ' || line[j - 1] == '\t'))
			j--;
		if (j == 0)
			return UNDEFINED;
		size_t i = j;
		while (i > 0 && isIdentChar(line[i - 1]))
			i--;
		if (i == j)
			return UNDEFINED;
		if (!isIdentStart(line[i]))
			return UNDEFINED;
		string name = line.substr(i, j - i);
		return name;
	}
}





ShaderFunctionModule::ShaderFunctionModule(const Path &path, const string &name)
: CodeFileDescriptor(path, true), name(name) {
	if (exists())
		updateModuleContent();
}


void ShaderFunctionModule::updateModuleContent() {
	if (!exists())
		THROW(FileNotFoundError, getFilename(), getDirectory().string());

	string content = readCode();
	string noComments = removeComments(content);
	string declsOnly = stripCurlyBlocks(noComments);

	functionNames.clear();
	structNames.clear();

	size_t start = 0;
	while (start <= declsOnly.size()) {
		size_t end = declsOnly.find('\n', start);
		if (end == string::npos)
			end = declsOnly.size();
		string line = trim(declsOnly.substr(start, end - start));
		if (!line.empty()) {
			if (isPreprocessorLine(line)) {
				// skip
			}
			else if (startsWithStruct(line)) {
				optional<string> sname = extractStructName(line);
				if (sname.has_value())
					structNames.push_back(*sname);
			}
			else {
				optional<string> fname = extractFunctionName(line);
				if (fname.has_value())
					functionNames.push_back(*fname);
			}
		}
		if (end == declsOnly.size())
			break;
		start = end + 1;
	}
}


bool ShaderFunctionModule::containsFunction(const string &functionName) {
	THROW_IF(!exists(), FileNotFoundError, getFilename(), getDirectory().string());

	if (functionNames.empty() and structNames.empty())
		updateModuleContent();

	return std::find(functionNames.begin(), functionNames.end(), functionName) != functionNames.end();
}


bool ShaderFunctionModule::containsStruct(const string &structName) {
	THROW_IF(!exists(), FileNotFoundError, getFilename(), getDirectory().string());

	if (functionNames.empty() and structNames.empty())
		updateModuleContent();

	return std::find(structNames.begin(), structNames.end(), structName) != structNames.end();
}

ShaderCodeGenerator::ShaderCodeGenerator(const Path &baseFile, const Path &targetFile, const string &macro)
: baseFile(baseFile), targetFile(targetFile), macroKey(macro) {}

void ShaderCodeGenerator::addModule(const ShaderFunctionModule &module) {
	if (!module.exists())
		THROW(FileNotFoundError, module.getPath().string());
	moduleCode += "\n" + module.readCode() + "\n";
	functionNames.insert(functionNames.end(), module.functionNames.begin(), module.functionNames.end());
	structNames.insert(structNames.end(), module.structNames.begin(), module.structNames.end());
}

void ShaderCodeGenerator::addStruct(const string &code, const string &structName) {
	if (contains(structNames, structName))
		return;
	structCode += "\n" + code + "\n";
	structNames.push_back(structName);
}

void ShaderCodeGenerator::addFunction(const string &code, const string &functionName) {
	if (contains(functionNames, functionName))
		THROW(ValueError, "Function " + functionName + " already exists in the generator");
	functionCode += "\n" + code + "\n";
	functionNames.push_back(functionName);
}

void ShaderCodeGenerator::addUniform(const string &code, const string &uniformName) {
	if (contains(uniformNames, uniformName))
		THROW(ValueError, "Uniform " + uniformName + " already exists in the generator");
	uniformCode += "\n" + code + "\n";
	uniformNames.push_back(uniformName);
}

void ShaderCodeGenerator::generate() const {
	if (!baseFile.exists())
		THROW(FileNotFoundError, baseFile.getPath().string());
	string insertCode = moduleCode + "\n" + structCode + "\n" + uniformCode + "\n" + functionCode + "\n";
	string baseCode = baseFile.readCode();
	CodeMacro macro = CodeMacro(macroKey, insertCode);
	string finalCode = macro.apply(baseCode);
	targetFile.writeCode(finalCode);
}
