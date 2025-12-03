#include "configFiles.hpp"
#include "exceptions.hpp"


Path ConfigFile::getRoot() const {
	return Path(at("root_rel_build"));
}

Path ConfigFile::getMainShaderDirectory() const {
	return Path(at("shader-templates"));
}

Path ConfigFile::getSDFMacroShader() const {
	return Path(at("sdfTools"));
}

Path ConfigFile::getMathToolsShaderMacro() const {
	return Path(at("mathTools"));
}

Path ConfigFile::getLightToolsShaderMacro() const {
	return Path(at("lightTools"));
}

Path ConfigFile::getStructsShaderMacro() const {
	return Path(at("glslStructs"));
}

Path ConfigFile::getShadersDir() const {
	return Path(at("shadersDir"));
}

Path ConfigFile::getScreenshotsDir() const {
	return Path(at("screenshotsDir"));
}

Path ConfigFile::getProjectRoot() {
	static Path root = []() {
		Path current = filesystem::current_path();
		while (!filesystem::exists(current / "config" / "config.json")) {
			if (current.has_parent_path() && current != current.parent_path()) {
				current = current.parent_path();
			}
			else {
				THROW(FileSystemError, "Could not find project root (looked for config/config.json)");
			}
		}
		return current;
	}();
	return root;
}
