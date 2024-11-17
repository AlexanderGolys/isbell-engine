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

#include <map>
#include <string>
#include <memory>

class Texture;

void error_callback(int error, const char* description);
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);

GLuint bindVAO();
void disableAttributeArrays(int how_many=4);
glm::mat4 generateMVP(glm::vec3 camPosition, glm::vec3 camLookAt, glm::vec3 upVector, float fov, float aspectRatio, float clippingRangeMin, float clippingRangeMax, glm::mat4 modelTransform);
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
	void setCallbacks(GLFWkeyfun* keyCallback = nullptr, GLFWcharfun* charCallback = nullptr, GLFWmousebuttonfun* mouseButtonCallback = nullptr, GLFWcursorposfun* cursorPosCallback = nullptr, GLFWcursorenterfun* cursorEnterCallback = nullptr, GLFWscrollfun* scrollCallback = nullptr, GLFWdropfun* dropCallback = nullptr);
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
public:
	GLuint programID;
	const char* vertex_file_path;
	const char* fragment_file_path;
	const char* geometry_file_path = nullptr;
	ShaderType shaderType;
	std::map<std::string, GLuint> uniformLocations;
	std::map<std::string, GLSLType> uniformTypes;

	Shader(const char* vertex_file_path, const char* fragment_file_path);
	Shader(const char* vertex_file_path, const char* fragment_file_path, const char* geometry_file_path);
	explicit Shader(std::string standard_file_path);
	~Shader();
	void use();
	void initUniforms(std::map<std::string, GLSLType> uniforms);

	void setTextureSampler(Texture* texture, int textureSlot);
	void setUniforms(std::map<std::string, const GLfloat*> uniformValues);
	void setUniform(std::string uniformName, const GLfloat* uniformValue);
	void setUniform(std::string uniformName, float uniformValue);
	void setUniform(std::string uniformName, int uniformValue);
	void setUniform(std::string uniformName, glm::vec2 uniformValue);
	void setUniform(std::string uniformName, glm::vec3 uniformValue);
	void setUniform(std::string uniformName, glm::vec4 uniformValue);
	void setUniform(std::string uniformName, glm::mat2 uniformValue);
	void setUniform(std::string uniformName, glm::mat3 uniformValue);
	void setUniform(std::string uniformName, glm::mat4 uniformValue);
	void setUniform(std::string uniformName, float x, float y);
	void setUniform(std::string uniformName, float x, float y, float z);
	void setUniform(std::string uniformName, float x, float y, float z, float w);
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
    std::function<glm::vec3(float)> up;
	glm::mat4 projectionMatrix;

	Camera();

	Camera(glm::vec3 position, glm::vec3 lookAtPos, glm::vec3 upVector, float fov_x=PI/4, 
		   float aspectRatio=16/9.f, float clippingRangeMin=.01f, float clippingRangeMax=100.f);

	Camera(const std::shared_ptr<SmoothParametricCurve> &trajectory, glm::vec3 lookAtPos, glm::vec3 upVector, float fov_x = PI / 4,
	       float aspectRatio = 16 / 9.f, float clippingRangeMin=.01f, float clippingRangeMax=100.f);

    Camera(const std::shared_ptr<SmoothParametricCurve> &trajectory, const std::shared_ptr<SmoothParametricCurve> &lookAtPos, glm::vec3 upVector, float fov_x = PI / 4,
       float aspectRatio = 16 / 9.f, float clippingRangeMin=.01f, float clippingRangeMax=100.f);

    Camera(const std::shared_ptr<SmoothParametricCurve> &trajectory, const std::shared_ptr<SmoothParametricCurve> &lookAtPos, const std::function<glm::vec3(float)> &upVector, float fov_x = PI / 4,
   float aspectRatio = 16 / 9.f, float clippingRangeMin=.01f, float clippingRangeMax=100.f);

	glm::vec3 position(float t) { return trajectory->operator()(t); }
    glm::vec3 lookAtPoint(float t) { return lookAtFunc->operator()(t); }
    glm::vec3 upVector(float t) { return up(t); }
	glm::mat4 mvp(float t, const glm::mat4 &modelTransform);
	glm::mat4 viewMatrix(float t);
	glm::mat4 vp(float t);
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
	std::shared_ptr<Shader> shader;
	std::vector<std::shared_ptr<Attribute>> attributes;
	std::shared_ptr<Model3D> model;
	std::shared_ptr<SuperMesh> super = nullptr;
    std::shared_ptr<WeakSuperMesh> weak_super = nullptr;
    GLuint elementBufferLoc = 0;

	std::map<std::string, GLSLType> uniforms;
	std::map<std::string, std::shared_ptr<std::function<void(float, std::shared_ptr<Shader>)>>> uniformSetters;
	std::function<void(float)> customStep;

	explicit RenderingStep(const std::shared_ptr<Shader> &shader);
	RenderingStep(const RenderingStep &other);
	~RenderingStep();

	void setModel(const std::shared_ptr<Model3D> &model);
	void setSuperMesh(const std::shared_ptr<SuperMesh> &super);
    void setWeakSuperMesh(const std::shared_ptr<WeakSuperMesh> &super);
	void initStdAttributes();
	void initMaterialAttributes();
    void initElementBuffer();
	void resetAttributeBuffers();
	void initUnusualAttributes(std::vector<std::shared_ptr<Attribute>> attributes);
	void loadStandardAttributes(); // 0:position, 1:normal, 2:color, 3:uv
    void loadElementBuffer();
	void enableAttributes();
	void disableAttributes();
	void addCustomAction(const std::function<void(float)> &action);
	bool superLoaded() const { return super != nullptr; }
    bool weakSuperLoaded() const { return weak_super != nullptr; }
    void init(const std::shared_ptr<Camera> &cam, const std::vector<std::shared_ptr<PointLight>> &lights);
    void initTextures();
    void bindTextures();

	void addUniforms(const std::map<std::string, GLSLType> &uniforms, std::map<std::string, std::shared_ptr<std::function<void(float, std::shared_ptr<Shader>)>>> setters);
	void addUniform(std::string uniformName, GLSLType uniformType, std::shared_ptr<std::function<void(float, std::shared_ptr<Shader>)>> setter);
	void addConstFloats(const std::map<std::string, float>& uniforms);
	void addConstVec4(const std::string& name, glm::vec4 value);
	void addConstColor(std::string name, glm::vec4 value) { addConstVec4(name, value); }
	void setUniforms(float t);

    void addCameraUniforms(const std::shared_ptr<Camera>& camera);
	void addLightUniform(const std::shared_ptr<PointLight>& pointLight, int lightIndex=1);
	void addLightsUniforms(const std::vector<std::shared_ptr<PointLight>> &lights);
	void addMaterialUniform();
    void addTexturedMaterialUniforms();

	void renderStep(float t);

};

class Renderer {
private:
	float frameOlderTimeThanThePublicOne = 0;

public:
	std::unique_ptr<Window> window;
	GLuint vao;
	std::vector<std::shared_ptr<RenderingStep>> renderingSteps;
	std::shared_ptr<Camera> camera;
	std::vector<std::shared_ptr<PointLight>> lights;
	float time = 0;
    float dt = 0;
	float animSpeed;
	glm::vec4 bgColor;
	std::unique_ptr<std::function<void(float, float)>> perFrameFunction;
	

	explicit Renderer(float animSpeed=1.f, glm::vec4 bgColor=BLACK);

	Renderer(	int width, int height, 
				const char* title,
				const std::shared_ptr<Camera> &camera,
				const std::vector<std::shared_ptr<PointLight>> &lights,
				const std::vector<std::shared_ptr<RenderingStep>>& renderingSteps,
				float animSpeed=1.f, 
				glm::vec4 bgColor=BLACK);

	Renderer(	Resolution resolution, 
				const char* title,
				const std::shared_ptr<Camera> &camera,
				const std::vector<std::shared_ptr<PointLight>> &lights,
				const std::vector<std::shared_ptr<RenderingStep>>& renderingSteps,
				float animSpeed=1.f, 
				glm::vec4 bgColor=BLACK);

	~Renderer();
	
	void initMainWindow(int width, int height, const char* title);
	void initMainWindow(Resolution resolution, const char* title);

	void setCamera(const std::shared_ptr<Camera> &camera);
	void setLights(const std::vector<std::shared_ptr<PointLight>> &lights);

	void addRenderingStep(std::shared_ptr<RenderingStep> renderingStep);
	
	float initFrame();
	float lastDeltaTime() const;

	void addPerFrameUniforms(std::map<std::string, GLSLType> uniforms, std::map<std::string, std::shared_ptr<std::function<void(float, std::shared_ptr<Shader>)>>> setters);
	void addPerFrameUniform(std::string uniformName, GLSLType uniformType, std::shared_ptr<std::function<void(float, std::shared_ptr<Shader>)>> setter);
	void addConstUniforms(std::map<std::string, GLSLType> uniforms, std::map<std::string, std::shared_ptr<std::function<void(std::shared_ptr<Shader>)>>> setters);
	void addConstUniform(std::string uniformName, GLSLType uniformType, std::shared_ptr<std::function<void(std::shared_ptr<Shader>)>> setter);
	void addTimeUniform();
	void addConstFloats(std::map<std::string, float> uniforms);
	void addCustomAction(std::function<void(float)> action);
    void addCustomAction(std::function<void(float, float)> action);

	void initRendering();
	void renderAllSteps();
	int mainLoop();
};
