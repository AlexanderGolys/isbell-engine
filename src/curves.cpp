#include "common/glsl_utils.hpp"
#include "common/specific.hpp"
#include <memory>
#include <cmath>

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
	return .01f + .005f*sin(11*t)*sin(11*t);
}

MaterialPhong mat(float t)
{
	auto c = lerp(BLUE_PALLETTE[4], BLUE_PALLETTE[8], sin(11*t)*sin(11*t));
	return MaterialPhong(c, c, WHITE, .2531423, .956 , .45739731, 50.0);
}

float w2(float t)
{
	return .025f + .015f*sin(9*t);
}

MaterialPhong mater2(float t)
{
	auto c = lerp(BLUE_PALLETTE[5], BLUE_PALLETTE[9], sin(9*t)*sin(9*t));
	return MaterialPhong(c, c, WHITE, .08031423, .8656 , .0145739731, 80.0);
}

float w3(float t)
{
	return .02f + .01f*sin(2*t);
}

MaterialPhong mater3(float t)
{
	auto c = lerp(BLUE_PALLETTE[2], BLUE_PALLETTE[2], sin(27*t)*sin(27*t));
	return MaterialPhong(c, c, WHITE, .14531423, .8656 , 0.45739731, 200.0);
}


int main(void)
{
	auto M = Matrix<Complex, 2>(-I*.7+1/2.f, I*1.1-2.14, I/0.8+.07, I*0.8f+1.8);

    Renderer renderer = Renderer(.05f, vec4(.01f, .0109f, 0.0285f, 1.0f));
    renderer.initMainWindow(UHD, "flows");
    Shader geoShader1 = Shader("C:\\Users\\PC\\Desktop\\ogl-master\\src\\shaders\\hyperbolicAut.vert",
                    "C:\\Users\\PC\\Desktop\\ogl-master\\src\\shaders\\hyperbolicAut.frag");



	shared_ptr<SmoothParametricCurve> curve = make_shared<SmoothParametricCurve>([](float t) { return vec3(-1.5, -2.55, 2); }, .01f);
    shared_ptr<Camera> camera = make_shared<Camera>(curve, vec3(0.0, 0.0, 0), vec3(0, 0, 1), PI/4);
    auto lights = vector<shared_ptr<PointLight>>({std::make_shared<PointLight>(vec3(1, 3, -1), vec4(.95, .79, .861, 1), 15.0f),
                                                 std::make_shared<PointLight>(vec3(-3, -1, 3), vec4(.869, .9864, .864, 1), 18.0f),
                                                 std::make_shared<PointLight>(vec3(2, 3, 1.5), vec4(.9, 1, .9, 1), 15.0f)});


    MaterialPhong materialmid = MaterialPhong(BLUE_PALLETTE[7], BLUE_PALLETTE[7], BLUE_PALLETTE[2], .06031423, .222156, .025739731, 2.0); // mid
    MaterialPhong black = MaterialPhong(PALETTE2[3], PALETTE2[3], WHITE, .1545, .38, .275, 100.0); // right, reddish
	MaterialPhong yelue = MaterialPhong(PALETTE2[4], PALETTE2[4], WHITE, .22155, .95, .997, 5.0); // floor
	MaterialPhong yelue2 = MaterialPhong(PALETTE2[4]*1.05f, PALETTE2[4], WHITE, .29155, .95, .997, 30.0);

	vector<pair<Complex, Complex>> endpoints = {{I+.01, I*2+.02}, {I*.02+1, I*.02+3},
												{I*1.5f+.02, I*2.5f+.02}, {I*.02+1.5, I*.02+3.5},
												{I*2.f+.02, I*3.f+.02}, {I*.02+2.f, I*.02+4.f},
	{I+.01, I*.02+1}, {I*3.f+.025, I*.02+4.f}, {I/2.f+.01, I*.02+5.5f}};

	SuperPencilCurve circ = circlePencil(vec3(1, 0, 0), vec3(0, 1, 0), vec3(0, 0, 0), .4f, w, mater2, 500);
	float speed = 5.0f;

	circ.addDeformationAlongVectorField(VectorFieldR3([](vec3 v){return vec3(sin(v.y), cos(v.x), sin(v.x*v.y+1))/5.f;}));
	auto mesh = circ.associateMesh(100);

	for (int i = 0; i < 5; i++)
	{
		SuperCurve circ1 = circle(vec3(1, 0, 0), vec3(0, 1, 0), vec3(0, 0, 0), .5f+0.1f*i, w, mater2, 400);
		auto mesh1 = *circ1.associateMesh(60);
		mesh->merge(mesh1);
	}



	mesh->addPolyGroup(PlanarUnitDisk(100, 10).embeddInR3(-.015)*mat3(15), materialmid);
	mesh->precomputeBuffers();
	mesh->precomputeExtraBuffer("curvePoint");
	// circ.transformMeshByAmbientMap(SpaceEndomorphism::scaling(1, 1, 1));





	auto step = make_shared<RenderingStep>(make_shared<Shader>(geoShader1));
	step->setSuperMesh(mesh);

	auto deform = [&circ](float t){circ.transformMesh(t);};
	renderer.addCustomAction(deform);
	renderer.addRenderingStep(step);

    renderer.setCamera(camera);
    renderer.setLights(lights);
    return renderer.mainLoop();
}