
#include "buffer_utils.hpp"
#include <stdio.h>
#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>

#include <stdlib.h>
#include <string.h>

#include <GL/glew.h>

using namespace glm;
using std::vector, std::string, std::map, std::shared_ptr, std::unique_ptr, std::make_shared, std::make_unique;



void VBO::load(const void *data, size_t size, GLenum usage) {
    bind();
    glBufferData(GL_ARRAY_BUFFER, size, data, usage);
}

void EBO::load(const void *data, size_t size, GLenum usage) {
    bind();
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, size, data, usage);
}
void EBO::load(const std::vector<ivec3> &data, GLenum usage) {
    load(data.data(), data.size() * sizeof(ivec3), usage);
}

void VAO::linkAttribute(VBO vbo, const Attribute &attr, const void *offset) {
    bind();
    vbo.bind();
    glVertexAttribPointer(attr.inputNumber, lengthOfGLSLType(attr.type), primitiveGLSLType(attr.type), GL_FALSE, sizeOfGLSLType(attr.type), offset);
    glEnableVertexAttribArray(attr.inputNumber);
    vbo.unbind();
}
void VAO::linkAttribute(VBO vbo, int layout, int len, GLenum type, int blockSize, const void *offset) {
    bind();
    vbo.bind();
    glVertexAttribPointer(layout, len, type, GL_FALSE, blockSize, offset);
    glEnableVertexAttribArray(layout);
    vbo.unbind();
}
VAORenderingObject::VAORenderingObject(Shader *shader, WeakSuperMesh *mesh) {
    this->shader = shader;
    this->mesh = mesh;
    vao = VAO();
    stdVBO = VBO();
    matVBO = VBO();
    ex0VBO = VBO();
    exVBO = VBO();
    ebo = EBO();
}

void VAORenderingObject::loadStdVBO(GLenum usage) {
    vao.bind();
    stdVBO.bind();
    stdVBO.load(mesh->getBufferLocation(POSITION), mesh->getBufferLength(POSITION)*sizeof(Stds), usage);
}
void VAORenderingObject::loadMatVBO(GLenum usage) {
    matVBO.load(mesh->getBufferLocation(MATERIAL1), mesh->getBufferLength(MATERIAL1)*sizeof(mat4), usage);

}
void VAORenderingObject::loadEx0VBO(GLenum usage) {
    ex0VBO.load(mesh->getBufferLocation(EXTRA0), mesh->getBufferLength(EXTRA0)*sizeof(vec4), GL_STATIC_DRAW);
}
void VAORenderingObject::loadExVBO(GLenum usage) {
    exVBO.load(mesh->getBufferLocation(EXTRA1), mesh->getBufferLength(EXTRA1)*sizeof(mat4), GL_STATIC_DRAW);
}
void VAORenderingObject::loadEBO(GLenum usage) {
    ebo.load(mesh->bufferIndexLocation(), mesh->bufferIndexSize(), GL_STATIC_DRAW);
}
void VAORenderingObject::loadActiveAttrs(GLenum usage) {
    loadStdVBO(usage);
    loadEBO(usage);
    if (mesh->isActive(MATERIAL1))
        loadMatVBO(usage);
    if (mesh->isActive(EXTRA0))
        loadEx0VBO(usage);
    if (mesh->isActive(EXTRA1) || mesh->isActive(EXTRA2) || mesh->isActive(EXTRA3) || mesh->isActive(EXTRA4))
        loadExVBO(usage);
}
void VAORenderingObject::addUniform(std::string uniformName, GLSLType uniformType, const std::function<void(float, Shader &)> &setter) {
    shader->initUniforms({{uniformName, uniformType}});
    setters.push_back(setter);
}
void VAORenderingObject::setUniforms(float t) {
    shader->use();
    for (auto &setter : setters)
        setter(t, *shader);
}

void VAORenderingObject::linkStdVBO() {

    vao.linkAttribute(stdVBO, 0, 3, GL_FLOAT, sizeof(Stds), (void*)0);
    vao.linkAttribute(stdVBO, 1, 3, GL_FLOAT, sizeof(Stds), (void*)(sizeof(vec3)));
    vao.linkAttribute(stdVBO, 2, 2, GL_FLOAT, sizeof(Stds), (void*)(2*sizeof(vec3)));
    vao.linkAttribute(stdVBO, 3, 4, GL_FLOAT, sizeof(Stds), (void*)(2*sizeof(vec3)+sizeof(vec2)));
}
void VAORenderingObject::linkMatVBO() {
    vao.linkAttribute(matVBO, 4, 4, GL_FLOAT, sizeof(vec4x4), (void *)mesh->offset(MATERIAL1));
    vao.linkAttribute(matVBO, 5, 4, GL_FLOAT, sizeof(vec4x4), (void*)mesh->offset(MATERIAL2));
    vao.linkAttribute(matVBO, 6, 4, GL_FLOAT, sizeof(vec4x4), (void*)mesh->offset(MATERIAL3));
    vao.linkAttribute(matVBO, 7, 4, GL_FLOAT, sizeof(vec4x4), (void*)mesh->offset(MATERIAL4));
}
void VAORenderingObject::linkEx0VBO() {
    vao.linkAttribute(ex0VBO, 8, 4, GL_FLOAT, sizeof(vec4), (void *)0);
}
void VAORenderingObject::linkExVBO() {
    vao.linkAttribute(exVBO, 9, 4, GL_FLOAT, 0, (void*)mesh->offset(EXTRA1));
    vao.linkAttribute(exVBO, 10, 4, GL_FLOAT,0, (void*)mesh->offset(EXTRA2));
    vao.linkAttribute(exVBO, 11, 4, GL_FLOAT,0, (void*)mesh->offset(EXTRA3));
    vao.linkAttribute(exVBO, 12, 4, GL_FLOAT,0, (void*)mesh->offset(EXTRA4));
}

void VAORenderingObject::render(float t) {

    vao.bind();
    stdVBO.bind();
    ebo.bind();

    glDrawElements(GL_TRIANGLES, mesh->bufferIndexLength()*3, GL_UNSIGNED_INT, 0);

}

void VAORenderingObject::init() {
    shader->use();
    vao.bind();
    stdVBO.bind();
    glBufferData(GL_ARRAY_BUFFER, mesh->getBufferLength(POSITION)*sizeof(vec3), mesh->getBufferLocation(POSITION), GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Stds), (void*)0);
    glEnableVertexAttribArray(0);

    glBufferData(GL_ARRAY_BUFFER, mesh->getBufferLength(POSITION)*sizeof(vec3), mesh->getBufferLocation(NORMAL), GL_STATIC_DRAW);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Stds), (void*)sizeof(vec3));
    glEnableVertexAttribArray(1);

    glBufferData(GL_ARRAY_BUFFER, mesh->getBufferLength(POSITION)*sizeof(vec2), mesh->getBufferLocation(UV), GL_STATIC_DRAW);
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(Stds), (void*)(2*sizeof(vec3)));
    glEnableVertexAttribArray(2);

    glBufferData(GL_ARRAY_BUFFER, mesh->getBufferLength(POSITION)*sizeof(vec4), mesh->getBufferLocation(COLOR), GL_STATIC_DRAW);
    glVertexAttribPointer(3, 4, GL_FLOAT, GL_FALSE, sizeof(Stds), (void*)(2*sizeof(vec3) + sizeof(vec2)));
    glEnableVertexAttribArray(3);

    glBindBuffer(GL_ARRAY_BUFFER, 0);

    ebo.bind();
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, mesh->bufferIndexSize()*sizeof(GLuint)*3, mesh->bufferIndexLocation(), GL_STATIC_DRAW);

    vao.unbind();
    ebo.unbind();
    stdVBO.unbind();


    // loadActiveAttrs(GL_STATIC_DRAW);
    // linkStdVBO();

}

VAORenderer::VAORenderer(RENDER_SETTINGS settings, Camera *cam, std::vector<PointLight> lights) :
    window(settings.width, settings.height, settings.title), settings(settings), logger(settings.logging) {
    glewExperimental = true;
    if (glewInit() != GLEW_OK)
        exit(2138);
    camera = cam;
    this->lights = lights;
    time = 0;
    steps = std::vector<VAORenderingObject>();

}


void VAORenderer::addStep(const VAORenderingObject &step) {
    step.logger = settings.logging;
    steps.push_back(step);
    logger.log2("Added rendering step");
}
void VAORenderer::initCameraUniform() {
    auto setter = [cam=this->camera](float t, Shader &shader) {
        shader.setUniform("mvp", cam->mvp(t, mat4(1)));
    };
    for (auto &step : steps)
        step.addUniform("mvp", MAT4, setter);
    logger.log2("mvp uniform setted in " + std::to_string(steps.size()) + " steps.");


    auto setter2 = [cam=this->camera](float t, Shader &shader) {
        shader.setUniform("camPosition", cam->position(t));
    };
    for (auto &step : steps)
        step.addUniform("camPosition", VEC3, setter2);
    logger.log2("CamPosition uniform setted in " + std::to_string(steps.size()) + " steps.");

}

void VAORenderer::initLightsUniforms() {
    for (int i = 0; i < lights.size(); i++) {
        auto setter = [i, m=lights[i].compressToMatrix()](float t, Shader &shader) {
            shader.setUniform("light" + std::to_string(i+1), m);
        };
        for (auto &step : steps) {
            step.addUniform("light" + std::to_string(i+1), MAT4, setter);
        }
        logger.log2("light" + std::to_string(i+1) + " uniform setted in " + std::to_string(steps.size()) + " steps.");
    }
}
void VAORenderer::initTimeUniform() {
    auto setter = [](float t, Shader &shader) {
        shader.setUniform("t", t);
    };
    for (auto &step : steps)
        step.addUniform("t", FLOAT, setter);
    logger.log2("t uniform setted in " + std::to_string(steps.size()) + " steps.");

}

void VAORenderer::addPerFrameFunction(std::function<void(float)> function) {
    perFrameFunction = [f=perFrameFunction, g=function](float t) {
        f(t);
        g(t);
    };
    logger.log3("Custom per frame function updated.");
}

void VAORenderer::init() {
    // glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
    window.initViewport();
    if (settings.depthTest)
        glEnable(GL_DEPTH_TEST);
    logger.log2("Depth test " + std::string(settings.depthTest ? "enabled" : "disabled"));
    glShadeModel(settings.shading);
    logger.log2("Shading model setted to " + std::string(settings.shading == GL_FLAT ? "FLAT" : "SMOOTH"));
    initCameraUniform();
    initLightsUniforms();
    initTimeUniform();
    logger.log("Initialised camera, lights and time uniforms in each step.");
    for (auto &step : steps) {
        step.init();
        logger.log("Step initialised.");
    }
}

void VAORenderer::renderFrame(float dt) {
    time += dt;
    logger.log("New frame started.", "t = " + std::to_string(time), "(+  " + std::to_string(dt) + ")");

    glClearColor(settings.bgColor.x, settings.bgColor.y, settings.bgColor.z, 1.0f);

    if (settings.depthTest)
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    else
        glClear(GL_COLOR_BUFFER_BIT);


    perFrameFunction(time);
    for (auto &step : steps) {
        step.setUniforms(time);
        logger.log2("Uniforms setted.");

        step.render(time);
        logger.log2("Step rendered.");
    }
    glfwSwapBuffers(window.window);
    glfwPollEvents();
    logger.log("Buffers swapped.");
}

int VAORenderer::mainLoop() {
    init();
    logger.log3("Global init completed.");

    while (window.isOpen()) {
        time = glfwGetTime();
        renderFrame(time);
        logger.log2("Frame rendered.");
        }

    return window.destroy();
}
