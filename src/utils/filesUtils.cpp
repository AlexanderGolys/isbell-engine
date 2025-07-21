#include "filesUtils.hpp"

#include <fstream>
#include <io.h>
#include <iostream>
#include <sstream>

#include "exceptions.hpp"
#include "metaUtils.hpp"
#include "randomUtils.hpp"

using namespace glm;



CodeFileDescriptor::CodeFileDescriptor(const string &filename, const Path &directory, bool rootRelative)
: filename(filename),
  directory(rootRelative ? ConfigFile().getRoot() + directory : directory) {}

CodeFileDescriptor::CodeFileDescriptor(const Path &path)
: CodeFileDescriptor(path.filename(), path.pureDirectory(), true) {}

CodeFileDescriptor::CodeFileDescriptor(const string &filename, const string &directory, bool rootRelative)
: CodeFileDescriptor(filename, Path(directory), rootRelative) {}

CodeFileDescriptor::CodeFileDescriptor(const string &path, bool rootRelative) {
	Path p(path);
	if (rootRelative) p = p.makeAbsolute(ConfigFile().getRoot());
	directory = p.pureDirectory();
	filename = p.filename();
}


Path CodeFileDescriptor::getPath() const { return directory + filename; }
string CodeFileDescriptor::getFilename() const { return filename; }
Path CodeFileDescriptor::getDirectory() const { return directory; }



Path::Path() {
	char cwd[1024];
	if (getcwd(cwd, sizeof(cwd)) != nullptr)
		path = string(cwd);

	else throw SystemError("Current directory not found");
}

Path::Path(const string &path)
: path(path) {}

bool Path::unixStyle() const {
	int unix = 0, windows = 0;
	for (char c: path) {
		if (c == '/') unix++;
		if (c == '\\') windows++;
	}
	return unix > windows;
}

Path::Path(const Path &other)
: path(other.path) {}


Path::Path(Path &&other) noexcept
: path(std::move(other.path)) {}

Path &Path::operator=(const Path &other) {
	if (this == &other) return *this;
	path = other.path;
	return *this;
}

Path &Path::operator=(Path &&other) noexcept {
	if (this == &other) return *this;
	path = std::move(other.path);
	return *this;
}

string Path::slash() const {
	return unixStyle() ? "/" : "\\";
}

string Path::wrongSlash() const {
	return unixStyle() ? "\\" : "/";
}

Path::operator string() const {
	return path;
}

Path Path::operator+(const string &other) const {
	return Path(path + slash() + other);
}
Path Path::operator+(const Path &other) const {
	return Path(path + slash() + other.path);
}

void Path::switchStyle() {
	path.replace(path.begin(), path.end(), slash()[0], wrongSlash()[0]);
}

void Path::setStyle(bool newStyle) {
	if (newStyle != unixStyle())
		switchStyle();
}


bool Path::endsWithFile(const string &filename) const {
	if (!path.contains('.')) return false;
	int dotIndex = path.find_last_of('.');
	return dotIndex != 0 && path[dotIndex - 1] != slash()[0];
}




Path Path::pureDirectory() const {
	if (endsWithFile(path)) return Path(path.substr(0, path.find_last_of(slash())));
	return *this;
}




string Path::filename() const {
	if (endsWithFile(path)) return path.substr(path.find_last_of(slash()) + 1, path.size());
	throw ValueError("Path does not end with a file");
}




string Path::fileExtension() const {
	if (endsWithFile(path)) return path.substr(path.find_last_of('.') + 1, path.size());
	throw ValueError("Path does not end with a file");
}

Path Path::makeRelative(const Path &other) const {
	if (!path.contains(other.path)) throw ValueError("Path does not contain the other path");
	return Path(path.substr(other.path.size()));
}

Path Path::makeAbsolute(const Path &root) const {
	return root + *this;
}
Path Path::makeRelative(const string &other) const { return makeRelative(Path(other)); }

Path Path::makeAbsolute(const string &root) const { return makeAbsolute(Path(root)); }
Path Path::goUp() const { return Path(path.substr(0, path.find_last_of(slash()))); }

string Path::to_str() const { return path; }

// Path::operator string() const { return path; }


CodeMacro::CodeMacro(const string &key, const string &code)
: replacementKey(key),
  replacementCode(code) {}

CodeMacro::CodeMacro(const string &key, const CodeFileDescriptor &file)
: replacementKey(key) {
	CodeFileDescriptor f = file;
	replacementCode = f.getCode();
}

CodeMacro::CodeMacro(const CodeMacro &other)
: replacementKey(other.replacementKey),
  replacementCode(other.replacementCode) {}

CodeMacro::CodeMacro(CodeMacro &&other) noexcept
: replacementKey(std::move(other.replacementKey)),
  replacementCode(std::move(other.replacementCode)) {}

CodeMacro &CodeMacro::operator=(const CodeMacro &other) {
	if (this == &other) return *this;
	replacementKey = other.replacementKey;
	replacementCode = other.replacementCode;
	return *this;
}

CodeMacro &CodeMacro::operator=(CodeMacro &&other) noexcept {
	if (this == &other) return *this;
	replacementKey = std::move(other.replacementKey);
	replacementCode = std::move(other.replacementCode);
	return *this;
}

string CodeMacro::apply(const string &codeScheme) const {
	string result = codeScheme;
	size_t pos = 0;
	while ((pos = result.find(replacementKey, pos)) != std::string::npos) {
		result.replace(pos, replacementKey.length(), replacementCode);
		pos += replacementCode.length();
	}
	return result;
}

string CodeMacro::getKey() const { return replacementKey; }


void CodeFileDescriptor::changeLine(const string &line, int lineNumber) {
	std::ifstream file(getPath().to_str());
	if (!exists()) throw FileNotFoundError(getFilename());
	std::string code;
	int i = 0;
	while (std::getline(file, code)) {
		if (i == lineNumber) {
			code = line;
			break;
		}
		i++;
	}
	file.close();
	writeCode(code);
}

string CodeFileDescriptor::readLine(int lineNumber) const {
	std::ifstream file(getPath().to_str());
	if (!exists()) throw FileNotFoundError(getFilename());
	std::string code;
	int i = 0;
	while (std::getline(file, code)) {
		if (i == lineNumber) break;
		i++;
	}
	file.close();
	return code;
}

string ConfigFile::operator[](const string &key) const { return config.at(key); }
string ConfigFile::check(const string &key) const { return config.at(key); }
Path ConfigFile::getRoot() const { return Path(check("root")); }
Path ConfigFile::getMainShaderDirectory() const { return Path(check("shader-templates")); }
Path ConfigFile::getSDFMacroShader() const { return Path(check("sdfTools")); }
Path ConfigFile::getMathToolsShaderMacro() const { return Path(check("mathTools")); }
Path ConfigFile::getLightToolsShaderMacro() const { return Path(check("lightTools")); }
Path ConfigFile::getStructsShaderMacro() const { return Path(check("glslStructs")); }
Path ConfigFile::getShadersDir() const { return Path(check("shadersDir")); }
Path ConfigFile::getScreenshotsDir() const { return Path(check("screenshotsDir")); }



string CodeFileDescriptor::extension() const {
	return filename.substr(filename.find_last_of('.') + 1);
}

string CodeFileDescriptor::getCode() {
	std::ifstream shaderStream(getPath(), std::ios::in);
	if (shaderStream.is_open()) {
		string code;
		std::stringstream sstr;
		sstr << shaderStream.rdbuf();
		code = sstr.str();
		shaderStream.close();
		return code;
	}
	throw FileNotFoundError(filename);
}

string CodeFileDescriptor::readCode() const {
	std::ifstream shaderStream(getPath(), std::ios::in);
	if (shaderStream.is_open()) {
		string code;
		std::stringstream sstr;
		sstr << shaderStream.rdbuf();
		code = sstr.str();
		shaderStream.close();
		return code;
	}
	throw FileNotFoundError(filename);
}

bool CodeFileDescriptor::exists() const {
	std::ifstream shaderStream(getPath(), std::ios::in);
	return shaderStream.is_open();
}

void CodeFileDescriptor::writeCode(const string &code) const {
	if (exists()) std::cout << "Overwriting code in file " << getFilename() << std::endl;
	std::ofstream shaderStream(getPath(), std::ios::out);
	if (shaderStream.is_open()) {
		shaderStream << code;
		shaderStream.close();
		return;
	}
	throw FileNotFoundError(filename);
}

void CodeFileDescriptor::modifyCode(const string &code) const {
	if (!exists())
		throw FileNotFoundError(filename);
	writeCode(code);
}

void CodeFileDescriptor::saveNewCode(const string &code) const {
	if (exists())
		throw SystemError("File " + getFilename() + " already exists in " + getDirectory().to_str() + ".");
	writeCode(code);
}

bool CodeFileDescriptor::recogniseDirectoryNamingStyle() {
	return directory.to_str().find('/') != string::npos;
}




TemplateCodeFile::TemplateCodeFile(const CodeFileDescriptor &templateCode, const CodeFileDescriptor &target, const vector<CodeMacro> &macros)
: templateCode(templateCode),
  target(target),
  macros(macros) {}

string TemplateCodeFile::generateCodeText() const {
	string code = templateCode.readCode();
	for (const auto &macro: macros) code = macro.apply(code);
	return code;
}

void TemplateCodeFile::render() const {
	target.writeCode(generateCodeText()); }

CodeFileDescriptor TemplateCodeFile::generatedCodeFile() const {
	render();
	return CodeFileDescriptor(target);
}

void TemplateCodeFile::addMacro(const CodeMacro &macro) {
	for (auto &m: macros)
		if (m.getKey() == macro.getKey()) {
			m = macro;
			return;
		}
	macros.push_back(macro);
}

TemplateCodeFile TemplateCodeFile::operator+(const CodeMacro &other) const {
	TemplateCodeFile res = *this;
	res.addMacro(other);
	return res;
}

void TemplateCodeFile::operator+=(const CodeMacro &other) {
	addMacro(other); }

void TemplateCodeFile::operator+=(const TemplateCodeFile &other) {
	for (const auto &macro: other.macros) addMacro(macro);
	target = other.target;
}

TemplateCodeFile TemplateCodeFile::operator+(const TemplateCodeFile &other) const {
	TemplateCodeFile res = *this;
	res += other;
	return res;
}


ShaderMethodTemplate::ShaderMethodTemplate(const string &returnType, const string &name, const string &arguments, const string &body, const string &returnCode, const std::map<string, string> &templateParameters)
: returnType(returnType),
  name(name),
  arguments(arguments),
  body(body),
  returnCode(returnCode),
  templateParameters(templateParameters) {}

ShaderMethodTemplate::ShaderMethodTemplate(const string &returnType, const string &name, const string &arguments, const string &returnCode, const std::map<string, string> &templateParameters)
: ShaderMethodTemplate(returnType, name, arguments, "", returnCode, templateParameters) {}

string ShaderMethodTemplate::generateBody() const {
	string result = body;
	for (const auto &[fst, snd]: templateParameters)
		result = replaceAll(result, fst, snd);
	return result;
}




string ShaderMethodTemplate::generateCall(const string &args) const {
return name + "(" + args + ")"; }

string ShaderMethodTemplate::generateCall(const string &arg1, const string &arg2) const { return name + "(" + arg1 + ", " + arg2 + ")"; }

string ShaderMethodTemplate::generateCall(const string &arg1, const string &arg2, const string &arg3) const {
	return name + "(" + arg1 + ", " + arg2 + ", " + arg3 + ")";
}

string ShaderMethodTemplate::randomRename() {
	name = name + "_" + randomStringNumeric();
	return name;
}

void ShaderMethodTemplate::changeName(const string &newName) { name = newName; }

void ShaderMethodTemplate::replaceFunctionCalls(const string &oldName, const string &newName) {
	body = replaceAll(body, oldName + "(", newName + "(");
	returnCode = replaceAll(returnCode, oldName + "(", newName + "(");
}

ShaderMethodTemplate ShaderMethodTemplate::composeReturnWith(const ShaderMethodTemplate &f2) const {
	ShaderMethodTemplate res = *this;
	res.setReturnCode(f2.generateCall(returnCode));
	res.setReturnType(f2.returnType);
	return res;
}

std::map<string, string> ShaderMethodTemplateFromUniform::generateTemplateDict() {
	std::map<string, string> result = {};
	for (int i = 0; i < keys.size(); i++)
		result[keys[i]] = parameterInCodeStr(i);
	return result;
}

ShaderMethodTemplateFromUniform::ShaderMethodTemplateFromUniform(const string &returnType, const string &name, const string &arguments, const string &body, const string &returnCode, const vector<string> &keys, int firstIndex, const string &arrayName)
: ShaderMethodTemplate(returnType, name, arguments, body, returnCode),
  arrayName(arrayName),
  _firstIndex(firstIndex),
  keys(keys)
  {
	templateParameters = generateTemplateDict();
}

ShaderMethodTemplateFromUniform::ShaderMethodTemplateFromUniform(const string &returnType, const string &name, const string &arguments, const string &returnCode)
: ShaderMethodTemplateFromUniform(returnType, name, arguments, "", returnCode) {}

void ShaderMethodTemplateFromUniform::setFirstIndex(int newIndex) {
	_firstIndex = newIndex;
	templateParameters = generateTemplateDict();
}

void ShaderMethodTemplateFromUniform::increseFirstIndex(int n) { setFirstIndex(_firstIndex + n); }

int ShaderMethodTemplateFromUniform::firstFreeIndex() const { return _firstIndex + keys.size() ; }


ShaderMethodTemplateFromUniform ShaderMethodTemplateFromUniform::composeReturnWith(const ShaderMethodTemplate &f2) const { return ShaderMethodTemplateFromUniform(f2.getReturnType(), getName(), getArguments(), getBody(), f2.generateCall(getReturnCode()), keys, _firstIndex, arrayName); }

ShaderRealFunction::ShaderRealFunction(const string &name, const string &body, const string &returnCode)
: ShaderMethodTemplateFromUniform("float", name, "float x", body, returnCode) {}

ShaderRealFunction::ShaderRealFunction(const string &name, const string &returnCode)
: ShaderRealFunction(name, "", returnCode) {}

ShaderMethodTemplate operator&(const string &f1Name, const ShaderMethodTemplate &f2) {
	return ShaderMethodTemplate(f2.getReturnType(), f2.getName(), f2.getArguments(), f2.getBody(), f1Name + "(" + f2.getReturnCode() + ")", f2.templateParameters);
}

string ShaderMethodTemplate::generateCode() const
{
	return returnType + " " + name + "(" + arguments + ") {\n" + generateBody() + "\nreturn " + returnCode + ";\n}\n";
}

void ShaderMethodTemplate::applyRenaming(const string &oldName, const string &newName) {
	replaceFunctionCalls(oldName, newName);
	if (name == oldName) changeName(newName);
}

ShaderMethodTemplateFromUniform operator&(const string &f1Name, const ShaderMethodTemplateFromUniform &f2) {
	return ShaderMethodTemplateFromUniform(f2.getReturnType(), f2.getName(), f2.getArguments(), f2.getBody(), f1Name + "(" + f2.getReturnCode() + ")", f2.keys, f2._firstIndex, f2.arrayName);
}


int ShaderMethodTemplateFromUniform::parameterIndex(const string &key) const {
	return _firstIndex + indexOf(key, keys);
}

string ShaderMethodTemplateFromUniform::parameterInCodeStr(int i) const {
	return arrayName + "[" + str(i+_firstIndex) + "]";
}


int ShaderMethodTemplateFromUniform::firstIndex() const { return _firstIndex; }

vector<string> ShaderMethodTemplateFromUniform::getSortedKeys() const { return keys; }

string ShaderMethodTemplateFromUniform::getArrayName() const { return arrayName; }

void ShaderMethodTemplateFromUniform::addParameter(const string &key)
{
	keys.push_back(key);
	templateParameters = generateTemplateDict();
}

ShaderRealFunction operator&(const string &f1Name, const ShaderRealFunction &f2) {
	return ShaderRealFunction(f2.getName(), f2.getBody(), f1Name + "(" + f2.getReturnCode() + ")");
}

ShaderRealFunction ShaderRealFunction::composeReturnWith(const ShaderMethodTemplate &f2) const {
	return ShaderRealFunction(getName(), getBody(), f2.generateCall(getReturnCode()));
}

ShaderBinaryOperator::ShaderBinaryOperator(const string &name, const string &body, const string &returnCode)
: ShaderMethodTemplateFromUniform("float", name, "float x, float y", body, returnCode) {}

ShaderBinaryOperator::ShaderBinaryOperator(const string &name, const string &returnCode)
: ShaderBinaryOperator(name, "", returnCode) {}

ShaderRealFunctionR3::ShaderRealFunctionR3(const string &name, const string &body, const string &returnCode, const vector<string> &keys, int firstIndex, const string &arrayName)
: ShaderMethodTemplateFromUniform("float", name, "vec3 x", body, returnCode, keys, firstIndex, arrayName) {}

ShaderRealFunctionR3 operator&(const string &f1Name, const ShaderRealFunctionR3 &f2) {
	return ShaderRealFunctionR3(f2.getName(), f2.getBody(), f1Name + "(" + f2.getReturnCode() + ")");
}

ShaderRealFunctionR3 ShaderRealFunctionR3::composeReturnWith(const ShaderMethodTemplate &f2) const {
	return ShaderRealFunctionR3(getName(), getBody(), f2.generateCall(getReturnCode()));
}
