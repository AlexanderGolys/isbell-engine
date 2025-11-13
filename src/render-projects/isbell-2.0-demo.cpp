#include "../core/include/controllers/controllers.hpp"
#include "specific.hpp"
#include "planarMesh.hpp"
#include "sceneRendering.hpp"


WindowSettings windowSettings = WindowSettings(Resolution::FHD, "window");
RenderSettings settings = RenderSettings(windowSettings, vec4(0.05f,0.05f,0.05,1), true, true, true, true, 2.f, 1.0f);
CameraSettings cameraSettings = CameraSettings();

Path shaderDir = Path("C:\\Users\\shitstem\\Desktop\\cppProjects\\isbell-engine\\src\\core\\shaders\\shaders");

int main() {
	Renderer renderer = Renderer(settings);

	sptr<Camera> camera = make_shared<Camera>(cameraSettings, vec3(-3, -3, 5), vec3(0, 0, 0), vec3(0, 0, 1));
	sptr<ShaderProgram> shader3D = ShaderProgram::standardShaderProgram(shaderDir / "basic.frag", shaderDir / "basic.vert");
	sptr<ShaderProgram> shader2D = ShaderProgram::standardShaderProgram(shaderDir / "ui.frag", shaderDir / "ui.vert");

	auto materialModel = make_shared<MaterialModelPhong>(
	vec4(0.2f, 0.2f, 0.8f, 1.0f),
	vec4(0.2f, 0.2f, 0.8f, 1.0f),
	vec4(0.2f, 0.5f, 0.8f, 1.0f),
	0.2f, 0.8f, 1.0f, 2.0f);

	auto flatsurf = SmoothParametricSurface([](float t, float u){
		return vec3(t, u, t*u/4-sin(2*t)*cos(3*u));
	}, vec2(-3, 2), vec2(-2, 2), false, false, 0.01f);

	auto ring = PlanarSurface([](vec2 uv){
		return vec2(uv.x*cos(uv.y), uv.x*sin(uv.y));
	}, vec2(.5, 1), vec2(0, TAU), false, true);

	auto ringMesh = make_shared<BasicPlanarMesh>(ring, 15, 50);
	auto mesh = make_shared<IndexedMesh3D>(flatsurf, 100, 300);

	auto light1 = make_shared<SimplePointLight>(vec3(5, -2, 13), WHITE, vec3(1, 0.01f, 0.001f));
	auto light2 = make_shared<SimplePointLight>(vec3(-5, 5, 15), WHITE, vec3(1, 0.01f, 0.001f));
	auto light3 = make_shared<SimplePointLight>(vec3(0, 3, 20), WHITE, vec3(1, 0.01f, 0.001f));

	auto meshLayer = make_shared<Mesh3DLayer>(shader3D, mesh, camera, materialModel, vector{light1, light2, light3});
	auto ringLayer = make_shared<GenericMeshLayer>(shader2D, ringMesh);

	auto tex = make_shared<Texture2D>(make_shared<ConstColorTextureData>(Color(0, 1, 1)), 0, "u_texture");
	ringLayer->addComponent(tex);

	auto camController = make_shared<CameraPositionController>(camera, [](vec3 pos, float t, float dt){
		return vec3(3, 0, 0)*cos(t) + vec3(0, 3, 0)*sin(t) + vec3(0, 0, 5);
	});
	meshLayer->addComponent(camController);

	renderer.addLayer(meshLayer);
	renderer.addLayer(ringLayer);
	renderer.run();
}
