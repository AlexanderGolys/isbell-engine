#pragma once
#include "filesUtils.hpp"
#include "metaUtils.hpp"

#define DEFAULT_CONFIG_PATH "../../../../config/config.json"

class JSONParser : public CodeFileDescriptor {
	using CodeFileDescriptor::CodeFileDescriptor;

public:
	string operator[](const string &key) const;
};

class ConfigFile{
	JSONParser jsonParser;
public:
	explicit ConfigFile(const string &configPath = DEFAULT_CONFIG_PATH) : jsonParser(configPath) {}
	string operator[](const string &key) const { return jsonParser[key]; }


	string at(const string &key) const {
		return jsonParser[key];
	}

	Path getRoot() const;
	Path getMainShaderDirectory() const;
	Path getSDFMacroShader() const;
	Path getMathToolsShaderMacro() const;
	Path getLightToolsShaderMacro() const;
	Path getStructsShaderMacro() const;
	Path getShadersDir() const;
	Path getScreenshotsDir() const;

	static Path getProjectRoot();
};
