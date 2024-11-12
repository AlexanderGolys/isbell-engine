#include "common/glsl_utils.hpp"
#include "common/specific.hpp"
#include <memory>
#include <cmath>
#include <stdlib.h>


using namespace glm;
using std::vector, std::string, std::map, std::shared_ptr, std::unique_ptr, std::pair, std::make_unique, std::make_shared;

#define spec_uni make_shared<std::function<void(float, shared_ptr<Shader>)>>
#define setter_uni shared_ptr<std::function<void(float, shared_ptr<Shader>)>>
#define shd shared_ptr<Shader>



float debugSin(float t)
{
	return -cos(t)/2 + 0.5;
}

float w(float t)
{
	return .007f + .00f*sin(1.5*t)*sin(1.5*t) + .01f/(1.5f + abs(.5*t));;
}



MaterialPhong mater(float t)
{
	auto c = lerp(REDPINK_PALLETTE[6], REDPINK_PALLETTE[6], sin(1.5*t)*sin(1.5*t));
	return MaterialPhong(c, c, REDPINK_PALLETTE[1], .2031423, .756641656 , .131160145739731, 10.0);
}

float w2(float t)
{
    return .005f + .00f*sin(1.5*t)*sin(1.5*t) + .005f/(1.f + abs(.5*t));;
}
float w1(float t)
{
    return .005f + .01f/(1.f + abs(.5*t));;
}


MaterialPhong mater2(float t)
{
    auto c = lerp(BLUE_PALLETTE[6], BLUE_PALLETTE[7], cos(1.5*t)*cos(1.5*t));
    return MaterialPhong(c, c, BLUE_PALLETTE[1], .2031423, .756641656 , .131160145739731, 10.0);
}




int main(void)
{
    Renderer renderer = Renderer(.05f, vec4(.08f, .0809f, 0.1385f, 1.0f));
    renderer.initMainWindow(FHD, "flows");

	shared_ptr<SmoothParametricCurve> camCurve = make_shared<SmoothParametricCurve>([](float t) { return vec3(cos(t*1)*sqrt(8), sin(t*1)*2, sin(t*1)*2); }, 0, TAU, .1, true, .01);
    shared_ptr<Camera> camera = make_shared<Camera>(camCurve, vec3(0.0, 0.0, 0), vec3(0, 0, 1), PI/4);
    auto lights = vector({std::make_shared<PointLight>(vec3(-0.4, -1, 3), vec4(.95, .79, .861, 1), 15.0f),
                          std::make_shared<PointLight>(vec3(0, 3, -0.1), vec4(.9498769, .69864, .694764, 1), 15.0f),
                          std::make_shared<PointLight>(vec3(3, -3, 0.5), vec4(.98, .98, .938, 1), 20.0f)});

    auto m1 =  MaterialPhong(REDPINK_PALLETTE[3], REDPINK_PALLETTE[3],REDPINK_PALLETTE[1], .2031423, .756641656 , .131160145739731, 10.0);
    auto m2 =  MaterialPhong(REDPINK_PALLETTE[7], REDPINK_PALLETTE[7], REDPINK_PALLETTE[1], .1031423, .956641656 , .131160145739731, 10.0);
    auto m3 =  MaterialPhong(BLUE_PALLETTE[8], BLUE_PALLETTE[9], BLUE_PALLETTE[1], .2031423, .756641656 , .131160145739731, 10.0);


    auto step = make_shared<RenderingStep>(make_shared<Shader>(
            "C:\\Users\\PC\\Desktop\\ogl-master\\src\\shaders\\hyperbolicAut.vert",
            "C:\\Users\\PC\\Desktop\\ogl-master\\src\\shaders\\hyperbolicAut.frag"));
    auto step2 = make_shared<RenderingStep>(make_shared<Shader>(
                "C:\\Users\\PC\\Desktop\\ogl-master\\src\\shaders\\hyperbolicAut.vert",
                "C:\\Users\\PC\\Desktop\\ogl-master\\src\\shaders\\hyperbolicAut.frag"));

    // WeakSuperMesh sph = icosphere(.9, 3, vec3(0, 0, 0), m1, PolyGroupID(222));

    PolyGroupID curveID = PolyGroupID(2139);
    auto supeer = make_shared<WeakSuperMesh>(icosphere(1, 3, vec3(0, 0, 0), m3, PolyGroupID(222)));
    auto curva1 = [&curveID](float t) { return SmoothParametricSurface(sphericalSpiral(.1 + .05*sin(4*t), 1.0f, TAU*2, curveID, .01), w2); };
    auto curva2 = [&curveID](float t) { return SmoothParametricSurface(sphericalSpiral(.2 + .07*sin(4.3*t+.41), 1.0f, TAU*2.2, curveID, .01), w1); };

    supeer->addUniformSurface(curva1(0), 200, 10, curveID, m1);
    auto supeer2 = make_shared<WeakSuperMesh>(curva2(0), 200, 10, curveID, m2);

    step->setWeakSuperMesh(supeer);
    renderer.addRenderingStep(step);
    step2->setWeakSuperMesh(supeer2);
    renderer.addRenderingStep(step2);
    renderer.setCamera(camera);
    renderer.setLights(lights);

    auto surfaceDeformer = [](float t, auto &curva) {
        return [t, &curva](BufferedVertex &v) {v.setPosition(curva(t).parametersNormalised(v.getUV()));};
    };


    auto deformer = [&surfaceDeformer, &supeer, &supeer2, &curveID, &curva1, &curva2](float t) {
        supeer->deformPerVertex(curveID, surfaceDeformer(t, curva1));
        supeer2->deformPerVertex(curveID, surfaceDeformer(t, curva2));
    };

	renderer.addCustomAction(deformer);

    return renderer.mainLoop();
}