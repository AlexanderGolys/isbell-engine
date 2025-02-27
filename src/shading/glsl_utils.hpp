#pragma once
#include "indexedRendering.hpp"

#include <cstdio>
#include <cstdlib>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>


// #include "renderingUtils.hpp"

#include <fstream>
#include <io.h>
#include <map>
#include <string>
#include <memory>
#include <sstream>

class Texture;

void error_callback(int error, const char* description);
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);

GLuint bindVAO();
void disableAttributeArrays(int how_many=4);
mat4 generateMVP(vec3 camPosition, vec3 camLookAt, vec3 upVector, float fov, float aspectRatio, float clippingRangeMin, float clippingRangeMax, mat4 modelTransform);
GLuint LoadShaders(const char* vertex_file_path, const char* fragment_file_path);
void setUniformTextureSampler(GLuint programID, Texture* texture, int textureSlot);


enum Resolution {
	FHD,
	UHD
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
	void setCallbacks(const GLFWkeyfun* keyCallback = nullptr, const GLFWcharfun* charCallback = nullptr, const GLFWmousebuttonfun* mouseButtonCallback = nullptr, GLFWcursorposfun* cursorPosCallback = nullptr, GLFWcursorenterfun* cursorEnterCallback = nullptr, GLFWscrollfun* scrollCallback = nullptr, GLFWdropfun* dropCallback = nullptr);
	bool isOpen();

	void renderFramebufferToScreen();
};



size_t sizeOfGLSLType(GLSLType type);
int lengthOfGLSLType(GLSLType type);
GLenum primitiveGLSLType(GLSLType type);

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
	explicit Shader(CodeFileDescriptor &file);
	explicit Shader(CodeFileDescriptor &&file);
	explicit Shader(const Path &file) : Shader(CodeFileDescriptor(file)) {}
	explicit Shader(const TemplateCodeFile &file) : Shader(file.generatedCodeFile()) {}

	Shader(const string &code, GLenum shaderType);

	Shader(const Shader &other);
	Shader(Shader &&other) noexcept;
	Shader & operator=(const Shader &other);
	Shader & operator=(Shader &&other) noexcept;
	virtual ~Shader() = default;


	GLuint getID() const;
	void compile();
	string getCode() const;

	static GLenum getTypeFromExtension(const string &extension);
};


class ShaderProgram {
protected:
	Shader vertexShader;
	Shader fragmentShader;
	std::optional<Shader> geometryShader;

	ShaderType shaderType;

public:
	std::map<string, GLuint> uniformLocations;
	std::map<string, GLSLType> uniformTypes;

	GLuint programID;

	ShaderProgram(const Shader &vertexShader, const Shader &fragmentShader);
	ShaderProgram(const Shader &vertexShader, const Shader &fragmentShader, const Shader &geometryShader);
	ShaderProgram(const string &vertexPath, const string &fragPath);


	void linkShaders();

//	ShaderProgram(const string &vertexShaderCode, const string &fragmentShaderCode, const string &geometryShaderCode="");
	explicit ShaderProgram(const std::string &standard_file_path);
	~ShaderProgram();
	void use();
	void initUniforms(const std::map<std::string, GLSLType> &uniforms);

	void setTextureSampler(const Texture* texture, int textureSlot) const;
	void setUniforms(const std::map<std::string, const GLfloat*> &uniformValues);
	void setUniform(const std::string &uniformName, const GLfloat* uniformValue);
	void setUniform(const std::string &uniformName, float uniformValue);
	void setUniform(const std::string &uniformName, int uniformValue);
	void setUniform(const std::string &uniformName, vec2 uniformValue);
	void setUniform(const std::string &uniformName, vec3 uniformValue);
	void setUniform(const std::string &uniformName, vec4 uniformValue);
	void setUniform(const std::string &uniformName, mat2 uniformValue);
	void setUniform(const std::string &uniformName, mat3 uniformValue);
	void setUniform(const std::string& uniformName, mat4 uniformValue);
	void setUniform(const std::string& uniformName, float x, float y);
	void setUniform(const std::string &uniformName, float x, float y, float z);
	void setUniform(const std::string &uniformName, float x, float y, float z, float w);
};

class Camera {
public:
	float fov_x;
	float aspectRatio;
	float clippingRangeMin;
	float clippingRangeMax;
	bool moving;
	std::shared_ptr<SmoothParametricCurve> trajectory;
    std::shared_ptr<SmoothParametricCurve> lookAtFunc;
    std::function<vec3(float)> up;
	mat4 projectionMatrix;

	Camera();

	Camera(vec3 position, vec3 lookAtPos, vec3 upVector=vec3(0, 0, 1), float fov_x=PI/4,
		   float aspectRatio=16/9.f, float clippingRangeMin=.01f, float clippingRangeMax=100.f);

	Camera(float radius, float speed, float height, vec3 lookAtPos=vec3(0), vec3 upVector=vec3(0, 0, 1), float fov_x=PI/4,
		   float aspectRatio=16/9.f, float clippingRangeMin=.01f, float clippingRangeMax=100.f);

	Camera(const std::shared_ptr<SmoothParametricCurve> &trajectory, vec3 lookAtPos, vec3 upVector, float fov_x = PI / 4,
	       float aspectRatio = 16 / 9.f, float clippingRangeMin=.01f, float clippingRangeMax=100.f);

    Camera(const std::shared_ptr<SmoothParametricCurve> &trajectory, const std::shared_ptr<SmoothParametricCurve> &lookAtPos, vec3 upVector, float fov_x = PI / 4,
       float aspectRatio = 16 / 9.f, float clippingRangeMin=.01f, float clippingRangeMax=100.f);

    Camera(const std::shared_ptr<SmoothParametricCurve> &trajectory, const std::shared_ptr<SmoothParametricCurve> &lookAtPos, const std::function<vec3(float)> &upVector, float fov_x = PI / 4,
   float aspectRatio = 16 / 9.f, float clippingRangeMin=.01f, float clippingRangeMax=100.f);

	vec3 position(float t) { return trajectory->operator()(t); }
    vec3 lookAtPoint(float t) { return lookAtFunc->operator()(t); }
    vec3 upVector(float t) { return up(t); }
	mat4 mvp(float t, const mat4 &modelTransform);
	mat4 viewMatrix(float t);
	mat4 vp(float t);
};

class Attribute {
public:
	std::string name;
	GLuint bufferAddress;
	size_t size;
	GLSLType type;
	int inputNumber;
	bool enabled;
	bool bufferInitialized;

	Attribute(std::string name, GLSLType type, int inputNumber);
	virtual ~Attribute();

	virtual void initBuffer();
	virtual void enable();
	virtual void disable();
	virtual void load(const void* firstElementAdress, int bufferLength);
	virtual void freeBuffer();
};






class RenderingStep {
    void weakMeshRenderStep(float t);
    void superMeshRenderStep(float t);
    void modelRenderStep(float t);

public:
	explicit RenderingStep(const std::shared_ptr<ShaderProgram> &shader);
	RenderingStep(const RenderingStep &other);

	virtual ~RenderingStep();
	std::shared_ptr<ShaderProgram> shader;
	std::vector<std::shared_ptr<Attribute>> attributes;
	std::shared_ptr<Model3D> model;
	std::shared_ptr<SuperMesh> super = nullptr;
    std::shared_ptr<WeakSuperMesh> weak_super = nullptr;
    GLuint elementBufferLoc = 0;

	std::map<std::string, GLSLType> uniforms;
	std::map<std::string, std::shared_ptr<std::function<void(float, std::shared_ptr<ShaderProgram>)>>> uniformSetters;
	std::function<void(float)> customStep;

	void setModel(const std::shared_ptr<Model3D> &model);
	void setSuperMesh(const std::shared_ptr<SuperMesh> &super);
    void setWeakSuperMesh(const std::shared_ptr<WeakSuperMesh> &super);

	int findAttributeByName(const std::string &name);
	void initStdAttributes();
	void initMaterialAttributes();
    void initElementBuffer();
	void resetAttributeBuffers();
	void initUnusualAttributes(const std::vector<std::shared_ptr<Attribute>>& attributes);
	void initExtraAttribute(int i);
	void loadMeshAttributes(); // 0:position, 1:normal, 2:color, 3:uv
	void initWeakMeshAttributes();
    void loadElementBuffer();
	void enableAttributes();
	void disableAttributes();

	void addCustomAction(const std::function<void(float)> &action);
	bool superLoaded() const { return super != nullptr; }
    bool weakSuperLoaded() const { return weak_super != nullptr; }
    virtual void init(const std::shared_ptr<Camera> &cam, const std::vector<Light> &lights);
    void initTextures();
    void bindTextures();

	void addUniforms(const std::map<std::string, GLSLType> &uniforms, std::map<std::string, std::shared_ptr<std::function<void(float, std::shared_ptr<ShaderProgram>)>>> setters);
	void addUniform(std::string uniformName, GLSLType uniformType, std::shared_ptr<std::function<void(float, std::shared_ptr<ShaderProgram>)>> setter);
	void addConstFloats(const std::map<std::string, float>& uniforms);
	void addConstVec4(const std::string& name, vec4 value);
	void addConstColor(const std::string &name, vec4 value) { addConstVec4(name, value); }
	void setUniforms(float t);

    virtual void addCameraUniforms(const std::shared_ptr<Camera>& camera);
	void addLightUniform(const Light &pointLight, int lightIndex = 1);
	void addLightsUniforms(const std::vector<Light> &lights);
    void addTexturedMaterialUniforms();

	virtual void renderStep(float t);

};








class Renderer {
	float frameOlderTimeThanThePublicOne = 0;
	float since_last_scr = 0;

public:
	std::unique_ptr<Window> window;
	GLuint vao;
	std::vector<std::shared_ptr<RenderingStep>> renderingSteps;
	std::shared_ptr<Camera> camera;
	std::vector<Light> lights;
	float time = 0;
    float dt = 0;
	Fooo animSpeed;
	vec4 bgColor;
	std::unique_ptr<std::function<void(float, float)>> perFrameFunction;
	std::string screenshotDirectory;
	float screenshotPeriod = 10;
	bool takeScreenshots = false;

	explicit Renderer(float animSpeed=1.f, vec4 bgColor=BLACK, const std::string &screenshotDirectory="screenshots/", float screenshotFrequency=-1);

	virtual ~Renderer();
	
	void initMainWindow(int width, int height, const char* title);
	void initMainWindow(Resolution resolution, const char* title);





	void setCamera(const std::shared_ptr<Camera> &camera);
	void setLights(const std::vector<Light> &lights);
	void setLightWithMesh(const Light &light, const MaterialPhong& material, const ShaderProgram &shader, float radius);
	void setLightsWithMesh(const std::vector<Light> &lights, const MaterialPhong& material, const ShaderProgram &shader, float radius);
	void setLightWithMesh(const Light &light, float ambient, float diff, float spec, float shine, const ShaderProgram &shader, float radius);
	void setLightsWithMesh(const std::vector<Light> &lights, float ambient, float diff, float spec, float shine, const ShaderProgram &shader, float radius);

	void addRenderingStep(std::shared_ptr<RenderingStep> renderingStep);
	void addMeshStep(const ShaderProgram& shader, const std::shared_ptr<WeakSuperMesh> &model, const MaterialPhong& material);

	float initFrame();
	float lastDeltaTime() const;

	void addPerFrameUniforms(const std::map<std::string, GLSLType> &uniforms, std::map<std::string, std::shared_ptr<std::function<void(float, std::shared_ptr<ShaderProgram>)>>> setters);
	void addPerFrameUniform(const std::string &uniformName, GLSLType uniformType, std::shared_ptr<std::function<void(float, std::shared_ptr<ShaderProgram>)>> setter);
	void addConstUniforms(const std::map<std::string, GLSLType>& uniforms, std::map<std::string, std::shared_ptr<std::function<void(std::shared_ptr<ShaderProgram>)>>> setters);
	void addConstUniform(const std::string &uniformName, GLSLType uniformType, std::shared_ptr<std::function<void(std::shared_ptr<ShaderProgram>)>> setter);
	void addTimeUniform();
	void addConstFloats(const std::map<std::string, float> &uniforms);
	void addCustomAction(std::function<void(float)> action);
    void addCustomAction(std::function<void(float, float)> action);
	void nonlinearSpeed(const Fooo &speed) { animSpeed = speed; }
	void screenshot(const std::string &filename) const;
	void screenshot() const;
	void addSurfaceFamilyDeformer(SurfaceParametricPencil &pencil, WeakSuperMesh &surface);

	virtual void initRendering();
	virtual void renderAllSteps();
	virtual int mainLoop();
};
