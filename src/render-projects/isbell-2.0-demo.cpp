#include "controllers.hpp"
#include "specific.hpp"
#include "planarMesh.hpp"
#include "sceneRendering.hpp"


WindowSettings windowSettings = WindowSettings(Resolution::HD2K, "window", false, false);
RenderSettings settings = RenderSettings(windowSettings, Color(0.05f,0.05f,0.05), true, true, true, true, 2.f, 1.0f);
CameraSettings cameraSettings = CameraSettings();

Path shaderDir = Path("C:\\Users\\shitstem\\Desktop\\cppProjects\\isbell-engine\\src\\core\\shaders\\shaders");

int main() {
	Renderer renderer = Renderer(settings);

	sptr<Camera> camera = make_shared<Camera>(cameraSettings, vec3(-6, -6, 5), vec3(0, 0, 0), vec3(0, 0, 1));
	sptr<ShaderProgram> shader3D = ShaderProgram::standardShaderProgram(shaderDir / "basic.frag", shaderDir / "basic.vert");
	sptr<ShaderProgram> shader2D = ShaderProgram::standardShaderProgram(shaderDir / "ui.frag", shaderDir / "ui.vert");

	auto materialModel1 = make_shared<MaterialModelPhong>(
	Color(0.2f, 0.2f, 0.8f),
	Color(0.2f, 0.2f, 0.8f),
	Color(0.2f, 0.5f, 0.8f),
	0.1f, 0.4f, 1.0f, 12.0f);

	auto materialModel2 = make_shared<MaterialModelPhong>(
	Color(0.8f, 0.2f, 0.2f),
	Color(0.95f, 0.2f, 0.2f),
	Color(1, 0.9f, 0.9f),
	0.08f, 0.6f, 2.0f, 22.0f);

	auto surf1 = DupinCyclide(3, 4, 1.3f, 0.4f);
	auto ellips = ellipsoid(2, 1, 5.3f, vec3(0, 0, 1));
	auto surf2 = SmoothParametricSurface([ellips](float u, float v) {
		vec3 p = ellips(u, v);
		return p + vec3((.5+p.z)*sin(2*p.z)/3, 0, 0);
	}, ellips.get_rangeT(), ellips.get_rangeU(), ellips.get_periodicT(), ellips.get_periodicU(), ellips.get_epsilon());

	auto mesh1 = make_shared<Mesh3D>(surf1, 501, 501);
	auto mesh2 = make_shared<Mesh3D>(surf2, 1500, 1500);
	// mesh2->translate(vec3(0, 12, 0));

	mesh1->transformVertices([&](const Vertex3D& v) {
		float z = v.position.z;
		if (z < 0 and dot(v.normal, vec3(0, 0, 1)) > 0 or z >= 0 and dot(v.normal, vec3(0, 0, 1)) <= 0)
			return Vertex3D(v.position, -v.normal, v.uv, v.color);
		return v;
		});


	auto light1 = SimplePointLight(vec3(5, -2, 14), WHITE, vec3(1, 0.01f, 0.001f));
	auto light2 = SimplePointLight(vec3(-5, 5, 15), WHITE, vec3(1, 0.01f, 0.001f));
	auto light3 = SimplePointLight(vec3(1, 1, 20), WHITE, vec3(1, 0.01f, 0.001f));
	auto lights = make_shared<arrayStruct<SimplePointLight>>(vector{light1, light2, light3});

	auto scene = make_shared<Scene3DLayer>(camera, lights);
	scene->addMeshLayer(shader3D, mesh1, materialModel1);
	scene->addMeshLayer(shader3D, mesh2, materialModel2);

	auto dist = *Distribution::inverseGamma(3, 1, 3);
	auto camShiftHor = make_shared<CameraHorisontalShiftOnKey>(camera, dist, GLFW_KEY_RIGHT, GLFW_KEY_LEFT, 0.1f);
	auto camShiftVer = make_shared<CameraVerticalShiftOnKey>(camera, dist, GLFW_KEY_UP, GLFW_KEY_DOWN, 0.1f);
	auto camRotationHor = make_shared<CameraHorisontalRotationOnKey>(camera, dist, GLFW_KEY_D, GLFW_KEY_A, 0.1f);
	auto camRotationVer = make_shared<CameraVerticalRotationOnKey>(camera, dist/5, GLFW_KEY_W, GLFW_KEY_S, 0.1f);
	auto camZoom = make_shared<CameraZoomOnScroll>(camera, dist, 0.5f);
	auto camRotate = make_shared<CameraRotateOnScroll>(camera, dist, 0.1f);

	scene->addUpdateComponent(camShiftHor);
	scene->addUpdateComponent(camShiftVer);
	scene->addUpdateComponent(camRotationHor);
	scene->addUpdateComponent(camRotationVer);
	scene->addUpdateComponent(camZoom);
	scene->addUpdateComponent(camRotate);

	scene->addEventListener(camShiftHor);
	scene->addEventListener(camShiftVer);
	scene->addEventListener(camRotationHor);
	scene->addEventListener(camRotationVer);
	scene->addEventListener(camZoom);
	scene->addEventListener(camRotate);

	renderer.addLayer(scene);
	renderer.run();
}
