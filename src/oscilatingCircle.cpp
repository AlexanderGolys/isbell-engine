#include "common/glsl_utils.hpp"
#include "common/specific.hpp"
#include <memory>
#include <cmath>
#include <stdlib.h>

using namespace glm;
using std::vector, std::string, std::map, std::shared_ptr, std::unique_ptr, std::pair, std::make_unique, std::make_shared;

float debugSin(float t)
{
	return -cos(t)/2 + 0.5;
}

float width1(float t)
{
    return .01f + .00f*sin(4*t)*sin(4*t);
}

float width2(float t)
{
    return .01f + .00f*sin(4*t)*sin(4*t);
}

int main(void)
{
    Renderer renderer = Renderer(.08f, vec4(.07f, .0409f, 0.05585f, 1.0f));
    renderer.initMainWindow(FHD, "flows");
    float camSpeed = 1.f;

    PolyGroupID curveID = PolyGroupID(420);
    PolyGroupID curvatureId = PolyGroupID(69);

    shared_ptr<SmoothParametricCurve> curve = make_shared<SmoothParametricCurve>([](float t) { return vec3(cos(t+.1), sin(t+.1), cos(5*t+1.5)/5.f); },

        curveID, 0, TAU, true, .01);
    SmoothParametricCurve curvatureCenter = SmoothParametricCurve([curve](float t) { return (*curve)(t) - curve->normal(t)*(min(2.f, curve->curvature_radius(t))); },
        curvatureId, curve->getT0(), curve->getT1(), curve->isPeriodic(), curve->getEps());

	shared_ptr<SmoothParametricCurve> camCurve = make_shared<SmoothParametricCurve>([curve](float t) { return vec3(cos(t)*1, sin(t)*1, .5)*3.f; }, 0, TAU, .1, true, .01);
    shared_ptr<SmoothParametricCurve> lookCamCurve = make_shared<SmoothParametricCurve>([curve](float t) { return vec3(cos(t), sin(t), 0); }, 0, TAU, .1, true, .01);

    shared_ptr<Camera> camera = make_shared<Camera>(camCurve, lookCamCurve, [curve](float t){return vec3(0, 0, 1);}, PI/4);
    auto lights = vector({std::make_shared<PointLight>(vec3(-0.4, -3, 2), vec4(.0952795, .785579, .4197285861, 1), 15.0f),
                          std::make_shared<PointLight>(vec3(0, -2, 2.1), vec4(.898769, .75369864, .03, 1), 15.0f),
                          std::make_shared<PointLight>(vec3(2, -1, -2), vec4(.698, .292598, .39785938, 1), 15.0f)});

    auto tex1 = make_shared<Texture>("C:\\Users\\PC\\Desktop\\ogl-master\\src\\textures\\texture1.bmp", 0, "texture_ambient", true);
    auto tex2 = make_shared<Texture>("C:\\Users\\PC\\Desktop\\ogl-master\\src\\textures\\texture1.bmp", 1, "texture_diffuse", true);
    auto tex3 = make_shared<Texture>("C:\\Users\\PC\\Desktop\\ogl-master\\src\\textures\\texture1.bmp", 2, "texture_specular", true);

    auto blue1 = make_shared<Texture>(vec3(BLUE_PALLETTE[2]), 0, "texture_ambient");
    auto blue2 = make_shared<Texture>(vec3(BLUE_PALLETTE[2]), 1, "texture_diffuse");
    auto blue3 = make_shared<Texture>(vec3(BLUE_PALLETTE[2]), 2, "texture_specular");

    auto red1 = make_shared<Texture>("C:\\Users\\PC\\Desktop\\ogl-master\\src\\textures\\texture_red.bmp", 0, "texture_ambient");
    auto red2 = make_shared<Texture>("C:\\Users\\PC\\Desktop\\ogl-master\\src\\textures\\texture_red.bmp", 1, "texture_diffuse");
    auto red3 = make_shared<Texture>("C:\\Users\\PC\\Desktop\\ogl-master\\src\\textures\\texture_red.bmp", 2, "texture_specular");

    auto matcurva =  MaterialPhong(std::move(red1), std::move(red2),std::move(red3), .282257031423, .666483956641656 , .5219131160145739731, 90.0);
    auto matsp =  MaterialPhong(std::move(tex1), std::move(tex2),std::move(tex3), 0.051031423, .3962453956641656 , .0931160145739731, 60.0);

    auto bluemat =  MaterialPhong(std::move(blue1), std::move(blue2),std::move(blue3), 0.051031423, .3962453956641656 , .0931160145739731, 60.0);

    auto step = make_shared<RenderingStep>(make_shared<Shader>(
            "C:\\Users\\PC\\Desktop\\ogl-master\\src\\shaders\\curvaCircle.vert",
            "C:\\Users\\PC\\Desktop\\ogl-master\\src\\shaders\\curvaCircle.frag"));
    auto step2 = make_shared<RenderingStep>(make_shared<Shader>(
            "C:\\Users\\PC\\Desktop\\ogl-master\\src\\shaders\\curvaCircle.vert",
            "C:\\Users\\PC\\Desktop\\ogl-master\\src\\shaders\\curvaCircle.frag"));

    PolyGroupID sphereID = PolyGroupID(421);

    shared_ptr<WeakSuperMesh> c1 = make_shared<WeakSuperMesh>(SmoothParametricSurface(*curve, width1), 150, 10, curveID);

    auto oscCirc = make_shared<Disk3D>("C:\\Users\\PC\\Desktop\\meshes\\wheel1.obj", vec3(0), vec3(1, 0, 0) ,vec3(0, 0, 1), sphereID);
    //
    // oscCirc.addUniformSurface(SmoothParametricSurface(curvatureCenter, width2), 150, 10, curvatureId);

    c1->addGlobalMaterial(matcurva);
    oscCirc->addGlobalMaterial(matsp);

    // supeer->addUniformSurface(curva1(0), 150, 10, curveID);

    step->setWeakSuperMesh(c1);
    step2->setWeakSuperMesh(oscCirc);
    // step->addUniform("intencities", MAT4, spec_uni([camera](float t, const shared_ptr<Shader> &shader) {shader->setUniform("mvp", camera->mvp(t));}));
    renderer.addRenderingStep(step2);
    renderer.addRenderingStep(step);


    renderer.setCamera(camera);
    renderer.setLights(lights);
    // renderer.addConstUniform("intencities", VEC4

    std::function<void(float, float)> def = [&oscCirc, sphereID, &curvatureCenter, &curve](float t, float dt){
        oscCirc->moveRotate(curvatureCenter(t), curve->tangent(t), curve->normal(t));
        oscCirc->setR(min(curve->curvature_radius(t), 2.f));
        // v.applyFunction(SpaceAutomorphism::deltaRotation(curve->normal(t-dt), curve->normal(t), curvatureCenter(t)));
        // v.applyFunction(SpaceAutomorphism::scaling(curve->curvature_radius(t)/curve->curvature_radius(t-dt), curvatureCenter(t)));
    };

	renderer.addCustomAction(def);

    return renderer.mainLoop();
}
