#include"Bases.h"
#include"Pipeline.h"
#include<Windows.h>

int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE HPrevInstance,
	LPSTR lpCmdLine, int nCmdShow) {
	IJWorld world;
	world.camera.position = IJVector(0, 2, 0, 0);
	world.camera.upwards = IJAuxVector(0, 0, 1);
	world.camera.direction = IJAuxVector(0, -1, 0);
	world.camera.type = IJ_ORTHOGRAPHIC;
	IJShape cube;
	cube.data[0] = IJVector(1.0, 0.0, 0.0, 0.0);
	cube.data[1] = IJVector(0.0, 1.0, 1.0, 0.0);
	cube.type = IJ_CUBE;
	world.shapes.push_back(cube);
	IJPatch *patch = VertexShaderStage1(world);
	patch = VertexShaderStage2(world, patch);
	MessageBox(NULL, "Hello World!", "Note", MB_OK);
	return 0;
}