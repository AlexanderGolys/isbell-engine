#pragma once


#include <cstdio>
#include <cstdlib>
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "macroParsing.hpp"
#include "buffers.hpp"
#include "sceneRendering.hpp"

#include "indexedRendering.hpp"
#include "renderingUtils.hpp"
#include "filesUtils.hpp"
#include "logging.hpp"


enum Resolution {
	FHD,
	HD2K,
	UHD,
	UNKNOWN
};

int predefinedWidth(Resolution res);
int predefinedHeight(Resolution res);


class Window {
public:
	GLFWwindow* window;
	int width;
	int height;
	float aspectRatio;
	Window(int width, int height, const char* title);
	Window(Resolution resolution, const char* title);
	~Window();
	int destroy();

	void initViewport();

	void showCursor();
	void disableCursor();
	void hideCursorWithinWindow();
	void stickyKeys(bool sticky);
	void stickyMouseButtons(bool sticky);
	void setCallbacks(const GLFWkeyfun* keyCallback = nullptr, const GLFWcharfun* charCallback = nullptr, const GLFWmousebuttonfun* mouseButtonCallback = nullptr,
					  GLFWcursorposfun* cursorPosCallback = nullptr, GLFWcursorenterfun* cursorEnterCallback = nullptr, GLFWscrollfun* scrollCallback = nullptr,
					  GLFWdropfun* dropCallback = nullptr);
	bool isOpen();
	void renderFramebufferToScreen();
};


enum ShaderType {
	CLASSIC,
	GEOMETRY1,
};


class Shader {
protected:
	GLenum shaderType;
	GLuint shaderID;
	string code = "";

public:
	explicit Shader(CodeFileDescriptor& file);
	explicit Shader(CodeFileDescriptor&& file);
	explicit Shader(const Path& file) : Shader(CodeFileDescriptor(file)) {}
	explicit Shader(const TemplateCodeFile& file) : Shader(file.generatedCodeFile()) {}
	Shader(const string& code, GLenum shaderType);
	Shader(const Shader& other);
	Shader(Shader&& other) noexcept;
	Shader& operator=(const Shader& other);
	Shader& operator=(Shader&& other) noexcept;
	virtual ~Shader() = default;

	GLuint getID() const;
	void compile();
	string getCode() const;

	static GLenum getTypeFromExtension(const string& extension);
};


class ShaderProgram {
protected:
	Shader vertexShader;
	Shader fragmentShader;
	optional<Shader> geometryShader;
	ShaderType shaderType;

public:
	unordered_map<string, GLuint> uniformLocations;
	unordered_map<string, GLSLType> uniformTypes;
	GLuint programID;

	ShaderProgram(const Shader& vertexShader, const Shader& fragmentShader);
	ShaderProgram(const Shader& vertexShader, const Shader& fragmentShader, const Shader& geometryShader);
	ShaderProgram(const string& vertexPath, const string& fragPath);
	explicit ShaderProgram(const string& standard_file_path);

	void linkShaders();
	~ShaderProgram();

	void use();
	void initUniforms(const unordered_map<string, GLSLType>& uniforms);
	void initTextureSampler(const Texture* texture);

	void setTextureSampler(const Texture* texture) const;
	void setUniforms(const unordered_map<string, const GLfloat*>& uniformValues);
	void setUniform(const string& uniformName, const GLfloat* uniformValue);
	void setUniform(const string& uniformName, float uniformValue);
	void setUniform(const string& uniformName, int uniformValue);
	void setUniform(const string& uniformName, vec2 uniformValue);
	void setUniform(const string& uniformName, vec3 uniformValue);
	void setUniform(const string& uniformName, vec4 uniformValue);
	void setUniform(const string& uniformName, mat2 uniformValue);
	void setUniform(const string& uniformName, mat3 uniformValue);
	void setUniform(const string& uniformName, mat4 uniformValue);
	void setUniform(const string& uniformName, float x, float y);
	void setUniform(const string& uniformName, float x, float y, float z);
	void setUniform(const string& uniformName, float x, float y, float z, float w);
};


class RenderingStep {
protected:
	shared_ptr<ShaderProgram> shader;
	vector<shared_ptr<AttributeBuffer>> attributes;
	shared_ptr<ElementBuffer> elementBuffer;
	shared_ptr<VertexArray> vao;

	shared_ptr<IndexedMesh> mesh = nullptr;
	shared_ptr<MaterialPhong> material = nullptr;

	std::map<string, GLSLType> uniforms;
	std::map<string, shared_ptr<std::function<void(float, shared_ptr<ShaderProgram>)>>> uniformSetters;
	HOM(float, void) customStep;
public:
	RenderingStep(const shared_ptr<ShaderProgram>& shader, const shared_ptr<MaterialPhong>& material, const shared_ptr<IndexedMesh>& mesh);
	virtual ~RenderingStep();

	void initAttributes();

	void initMaterialTextures() const;
	void bindTextures() const;

	void addCustomAction(const HOM(float, void)& action);
	virtual void init(const shared_ptr<Camera>& cam, const vector<Light>& lights);

	void addUniforms(const std::map<string, GLSLType>& uniforms, std::map<string, shared_ptr<std::function<void(float, shared_ptr<ShaderProgram>)>>> setters);
	void addUniform(string uniformName, GLSLType uniformType, shared_ptr<std::function<void(float, shared_ptr<ShaderProgram>)>> setter);
	void addConstFloats(const std::map<string, float>& uniforms);
	void addConstVec4(const string& name, vec4 value);
	void addConstColor(const string& name, vec4 value);
	void addMaterialUniforms();
	void setUniforms(float t);

	virtual void addCameraUniforms(const shared_ptr<Camera>& camera);
	void addLightUniform(const Light& pointLight, int lightIndex = 1);
	void addLightsUniforms(const vector<Light>& lights);

	virtual void renderStep(float t);
};


struct RenderSettings {
	vec4 bgColor;
	bool alphaBlending;
	bool depthTest;
	bool timeUniform;
	float speed;
	int maxFPS;
	Resolution resolution;
	string windowTitle;

	bool takeScreenshots;
	float screenshotFrequency;
	Path screenshotDirectory;

	RenderSettings(vec4 bgColor, bool alphaBlending, bool depthTest, bool timeUniform, float speed, int maxFPS, bool takeScreenshots, Resolution resolution, float screenshotFrequency = 0.f, const string& windowTitle = "window");
};


class Renderer {
protected:
	float last_time_capture = 0;

	unique_ptr<Window> window;
	vector<shared_ptr<RenderingStep>> renderingSteps;

	shared_ptr<Camera> camera;
	vector<Light> lights;

	float time = 0;
	float dt = 0;

	END(float) animSpeed;
	std::function<void(float, float)> perFrameFunction;
	RenderSettings settings;

public:
	explicit Renderer(float animSpeed = 1.f, vec4 bgColor = BLACK, const string& screenshotDirectory = "screenshots/", float screenshotFrequency = -1);
	explicit Renderer(const RenderSettings& settings);
	virtual ~Renderer();

	void initMainWindow(int width, int height, const char* title);
	void initMainWindow(Resolution resolution, const char* title);
	void initMainWindow();
	void resetTimer();

	void setCamera(const shared_ptr<Camera>& camera);
	void setLights(const vector<Light>& lights);

	void setLightWithMesh(const Light& light, const shared_ptr<MaterialPhong>& material, const ShaderProgram& shader, float radius);
	void setLightsWithMesh(const vector<Light>& lights, const shared_ptr<MaterialPhong>& material, const ShaderProgram& shader, float radius);
	void setLightWithMesh(const Light& light, float ambient, float diff, float spec, float shine, const ShaderProgram& shader, float radius);
	void setLightsWithMesh(const vector<Light>& lights_, float ambient, float diff, float spec, float shine, const ShaderProgram& shader, float radius);

	void addRenderingStep(shared_ptr<RenderingStep> renderingStep);
	void addMeshStep(const ShaderProgram& shader, const shared_ptr<IndexedMesh>& model, const shared_ptr<MaterialPhong>& material);

	float initFrame();
	float lastDeltaTime() const;

	void addPerFrameUniforms(const std::map<string, GLSLType>& uniforms, const std::map<string, shared_ptr<std::function<void(float, shared_ptr<ShaderProgram>)>>>& setters) const;
	void addPerFrameUniform(const string& uniformName, GLSLType uniformType, const shared_ptr<std::function<void(float, shared_ptr<ShaderProgram>)>>& setter) const;
	void addConstUniforms(const std::map<string, GLSLType>& uniforms, std::map<string, shared_ptr<std::function<void(shared_ptr<ShaderProgram>)>>> setters) const;
	void addConstUniform(const string& uniformName, GLSLType uniformType, shared_ptr<std::function<void(shared_ptr<ShaderProgram>)>> setter) const;
	void addTimeUniform() const;
	void addConstFloats(const std::map<string, float>& uniforms) const;
	void addCustomAction(HOM(float, void) action);
	void addCustomAction(std::function<void(float, float)> action);
	void nonlinearSpeed(const END(float)& speed);
	void addSurfaceFamilyDeformer(SurfaceParametricPencil& pencil, IndexedMesh& surface);

	void addFloorWorkingArea(vec3 corner1, vec3 corner2, vec3 corner3, vec3 corner4, const MaterialPhong& material, float height, float pyramid_height);

	virtual void initRendering();
	virtual void renderAllSteps();
	virtual int mainLoop();
};
