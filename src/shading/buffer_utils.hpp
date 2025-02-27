
#pragma once

#include "glsl_utils.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>



class VBO {
  GLuint id;
public:
  VBO() { glGenBuffers(1, &id); }
  VBO(const void* data, size_t size, GLenum usage) : VBO() { load(data, size, usage); }
  void destroy() { glDeleteBuffers(1, &id); }
  ~VBO() { destroy(); }
  void bind() { glBindBuffer(GL_ARRAY_BUFFER, id); }
  static void unbind() { glBindBuffer(GL_ARRAY_BUFFER, 0); }
  void load(const void* data, size_t size, GLenum usage);
  GLuint getID() { return id; }
};


class EBO {
  GLuint id;
public:
  EBO() { glGenBuffers(1, &id); }
  EBO(const void* data, size_t size, GLenum usage) : EBO() { load(data, size, usage); }
  void destroy() { glDeleteBuffers(1, &id); }
  ~EBO() { destroy(); }  void bind() { glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, id); }
  static void unbind() { glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0); }
  void load(const void* data, size_t size, GLenum usage);
  void load(const std::vector<glm::ivec3> &data, GLenum usage);

};

class VAO {
  GLuint id;
public:
  VAO() { glGenVertexArrays(1, &id); }
  void destroy() { glDeleteVertexArrays(1, &id); }
  ~VAO() { destroy(); }
  void bind() { glBindVertexArray(id); }
  static void unbind() { glBindVertexArray(0); }
  void linkAttribute(VBO vbo, const Attribute &attr, const void* offset);
  void linkAttribute(VBO vbo, int layout, int len, GLenum type, int blockSize, const void* offset);

};

enum LOGGING_INTENCITY {
  NOPE,
  TIMID,
  MEH___FINE,
  TRADITIONAL_MEDIA_WITHOUT_ADBLOCK
};

class LOGGING {
public:
  LOGGING_INTENCITY intencity;
  LOGGING(LOGGING_INTENCITY intencity) : intencity(intencity) {}
  bool operator>=(int level) { return intencity >= level; }
  bool operator<=(int level) { return intencity <= level; }
    bool operator==(int level) { return intencity == level; }
  void log (const std::string &m1, const std::string &m2, const std::string &m3) {
    if (intencity >= TIMID)
      std::cout << m1 << std::endl;
    if (intencity >= MEH___FINE)
      std::cout << m2 << std::endl;
    if (intencity >= TRADITIONAL_MEDIA_WITHOUT_ADBLOCK)
      std::cout << m3 << std::endl;
  }
  void log (const std::string &m2, const std::string &m3) {
    if (intencity >= MEH___FINE)
      std::cout << m2 << std::endl;
    if (intencity >= TRADITIONAL_MEDIA_WITHOUT_ADBLOCK)
      std::cout << m3 << std::endl;
  }
  void log (const std::string &m3) {
    if (intencity >= TRADITIONAL_MEDIA_WITHOUT_ADBLOCK)
      std::cout << m3 << std::endl;
  }
  void log2 (const std::string &m2, const std::string &m3) {
    if (intencity >= TIMID)
      std::cout << m2 << std::endl;
    if (intencity >= MEH___FINE)
      std::cout << m3 << std::endl;
  }
  void log2 (const std::string &m3) {
    if (intencity >= MEH___FINE)
      std::cout << m3 << std::endl;
  }
  void log3 (const std::string &m3) {
    if (intencity >= TIMID)
      std::cout << m3 << std::endl;
  }


};

class VAORenderingObject {
  Shader *shader;
  VAO vao;
  WeakSuperMesh *mesh;
  VBO stdVBO, matVBO, ex0VBO, exVBO;
  EBO ebo;
  std::vector<std::function<void(float, Shader&)>> setters = {};

public:
  mutable LOGGING logger = NOPE;
  VAORenderingObject(Shader *shader, WeakSuperMesh *mesh);

  void loadStdVBO(GLenum usage);
  void loadMatVBO(GLenum usage);
  void loadEx0VBO(GLenum usage);
  void loadExVBO(GLenum usage);
  void loadEBO(GLenum usage);
  void loadActiveAttrs(GLenum usage);

  void addUniform(std::string uniformName, GLSLType uniformType, const std::function<void(float, Shader&)> &setter);
  void setUniforms(float t);

  void linkStdVBO();
  void linkMatVBO();
  void linkEx0VBO();
  void linkExVBO();

  void init();
  void render(float t);
};



class RENDER_SETTINGS {
public:
  int width, height;
  const char* title;
  GLenum shading;
  bool depthTest;
  glm::vec3 bgColor;
  int maxFPS;
  LOGGING logging;

  RENDER_SETTINGS(int width, int height, const char* title, GLenum shading, bool depthTest, glm::vec3 bgColor, int maxFPS, LOGGING logging) : width(width), height(height), title(title), shading(shading), depthTest(depthTest), bgColor(bgColor), maxFPS(maxFPS), logging(logging) {}
  RENDER_SETTINGS(Resolution res, const char* title, GLenum shading, bool depthTest, glm::vec3 bgColor, int maxFPS, LOGGING logging) : RENDER_SETTINGS(predefinedWidth(res), predefinedHeight(res), title, shading, depthTest, bgColor, maxFPS, logging) {}
  float min_dt() { return 1.f / maxFPS; }
};

class VAORenderer {
  Window window;
  std::vector<VAORenderingObject> steps = {};
  Camera *camera;
  std::vector<PointLight> lights;
  float time=0;
  std::function<void(float)> perFrameFunction = [](float t) {};
  RENDER_SETTINGS settings;
  LOGGING logger;

public:
  VAORenderer(RENDER_SETTINGS settings, Camera *cam, std::vector<PointLight> lights);
  void addStep(const VAORenderingObject &step);
  void initCameraUniform();
  void initLightsUniforms();
  void initTimeUniform();
  void addPerFrameFunction(std::function<void(float)> function);
  void init();
  void renderFrame(float dt);
  int mainLoop();

};
