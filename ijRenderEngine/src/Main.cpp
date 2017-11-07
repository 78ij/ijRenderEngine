/*
*                 Creates a window,handle the messages.
*
*                       Created by 78ij in 2017.11
*/

#include"Bases.h"
#include"Pipeline.h"

extern BYTE buffer[WIDTH * HEIGHT * 3];
LRESULT CALLBACK WndProc(HWND, UINT, WPARAM, LPARAM);
IJWorld world;
int originx = 0, originy = 0;
bool isbuttondown;

static HDC screen_hdc;
static HWND screen_hwnd;
static HDC hCompatibleDC; 
static HBITMAP hCompatibleBitmap; 
static HBITMAP hOldBitmap;                 
static BITMAPINFO binfo;

int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE HPrevInstance,
	LPSTR lpCmdLine, int nCmdShow) {
	for (int i = 0; i < WIDTH * HEIGHT * 3; i++) {
		buffer[i] = 255;
	}
	world.camera.position = IJVector(1, 1, 0,1);
	world.camera.upwards = IJAuxVector(0, 0, 1);
	world.camera.direction = IJAuxVector(-1, -1, 0.4);
	world.camera.type = IJ_ORTHOGRAPHIC;
	IJShape sphere;
    sphere.data[0] = IJVector(0.0, 0.0, 0.0, 1.0);
	sphere.radius = 0.7;
	sphere.type = IJ_SPHERE;
	sphere.step[0] = 30;
	sphere.step[1] = 30;
	world.shapes.push_back(sphere);
	/*IJShape cube;
	cube.data[0] = IJVector(1.0, 0.0, 0.0, 1.0);
	cube.data[1] = IJVector(0.0, 1.0, 1.0, 1.0);
	cube.type = IJ_CUBE;
	world.shapes.push_back(cube);*/
	IJPatch *patch = VertexShaderStage1(world);
	patch = VertexShaderStage2(world, patch);
	RasterizationStage1(world, patch);
	FreePatch(patch, world);
	patch = NULL;
	//Line(IJVector(1,-1,0,0), IJVector(-1,1,0,0));
	static TCHAR szAppName[] = TEXT("BitBlt");
	HWND         hwnd;
	MSG          msg;
	WNDCLASS     wndclass;

	wndclass.style = CS_HREDRAW | CS_VREDRAW;
	wndclass.lpfnWndProc = WndProc;
	wndclass.cbClsExtra = 0;
	wndclass.cbWndExtra = 0;
	wndclass.hInstance = hInstance;
	wndclass.hIcon = LoadIcon(NULL, IDI_INFORMATION);
	wndclass.hCursor = LoadCursor(NULL, IDC_ARROW);
	wndclass.hbrBackground = (HBRUSH)GetStockObject(WHITE_BRUSH);
	wndclass.lpszMenuName = NULL;
	wndclass.lpszClassName = szAppName;

	if (!RegisterClass(&wndclass))
	{
		MessageBox(NULL, TEXT("This program requires Windows NT!"),
			szAppName, MB_ICONERROR);
		return 0;
	}

	hwnd = CreateWindow(szAppName, TEXT("Rendering Demo"),
		WS_OVERLAPPEDWINDOW^WS_THICKFRAME,
		CW_USEDEFAULT, CW_USEDEFAULT,
		WIDTH + 20, HEIGHT + 40,
		NULL, NULL, hInstance, NULL);

	ShowWindow(hwnd, nCmdShow);
	UpdateWindow(hwnd);

	ZeroMemory(&binfo, sizeof(BITMAPINFO));
	binfo.bmiHeader.biBitCount = 24; 
	binfo.bmiHeader.biCompression = BI_RGB;
	binfo.bmiHeader.biHeight = HEIGHT;
	binfo.bmiHeader.biPlanes = 1;
	binfo.bmiHeader.biSizeImage = 0;
	binfo.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
	binfo.bmiHeader.biWidth = WIDTH;

	screen_hwnd = hwnd;
	screen_hdc = GetDC(screen_hwnd);
	hCompatibleDC = CreateCompatibleDC(screen_hdc);
	hCompatibleBitmap = CreateCompatibleBitmap(screen_hdc, WIDTH, HEIGHT);
	hOldBitmap = (HBITMAP)SelectObject(hCompatibleDC, hCompatibleBitmap);

	//These two functions must be called every time we refresh the output.
	SetDIBits(screen_hdc, hCompatibleBitmap, 0, HEIGHT, buffer, (BITMAPINFO*)&binfo, DIB_RGB_COLORS);
	BitBlt(screen_hdc, -1, -1, WIDTH, HEIGHT, hCompatibleDC, 0, 0, SRCCOPY);

	while (GetMessage(&msg, NULL, 0, 0))
	{
		TranslateMessage(&msg);
		DispatchMessage(&msg);
	}
	return msg.wParam;
	return 0;
}

LRESULT CALLBACK WndProc(HWND hwnd, UINT message, WPARAM wParam, LPARAM lParam)
{
	static int  cxClient, cyClient, cxSource, cySource;
	HDC         hdcClient, hdcWindow;
	int         x, y;
	PAINTSTRUCT ps;

	switch (message)
	{
	case WM_CREATE:
		cxSource = GetSystemMetrics(SM_CXSIZEFRAME) +
			GetSystemMetrics(SM_CXSIZE);

		cySource = GetSystemMetrics(SM_CYSIZEFRAME) +
			GetSystemMetrics(SM_CYCAPTION);
		return 0;

	case WM_SIZE:
		cxClient = LOWORD(lParam);
		cyClient = HIWORD(lParam);
		return 0;
	
	case WM_LBUTTONDOWN:
		isbuttondown = true;
		originx = LOWORD(lParam);
		originy = HIWORD(lParam);
		break;
	case WM_MOUSEMOVE: {
		int x = LOWORD(lParam);
		int y = HIWORD(lParam);
		if (isbuttondown == true) {
			HDC hdc = GetDC(hwnd);
			double offsetx = (double)(x - originx) / 399.5;
			double offsety = (double)(y - originy) / 399.5;
			for (int i = 0; i < WIDTH * HEIGHT * 3; i++) {
				buffer[i] = 255;
			}
			world.camera.position = IJVector(
				world.camera.position[0] + offsetx,
				world.camera.position[1] - offsety,
				world.camera.position[2],
				world.camera.position[3]
			);
		IJPatch *patch = VertexShaderStage1(world);
		patch = VertexShaderStage2(world, patch);
		RasterizationStage1(world, patch);
		FreePatch(patch, world);
		patch = NULL;
		SetDIBits(screen_hdc, hCompatibleBitmap, 0, HEIGHT, buffer, (BITMAPINFO*)&binfo, DIB_RGB_COLORS);
		BitBlt(screen_hdc, -1, -1, WIDTH, HEIGHT, hCompatibleDC, 0, 0, SRCCOPY);
		}
		originx = x;
		originy = y;
		break;
	}
	case WM_LBUTTONUP:
		isbuttondown = false;
		break;

	case WM_DESTROY:
		PostQuitMessage(0);
		return 0;
	}
	return DefWindowProc(hwnd, message, wParam, lParam);
}