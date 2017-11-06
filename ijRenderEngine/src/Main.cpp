#include"Bases.h"
#include"Pipeline.h"
extern BYTE buffer[WIDTH * HEIGHT * 3];
LRESULT CALLBACK WndProc(HWND, UINT, WPARAM, LPARAM);

int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE HPrevInstance,
	LPSTR lpCmdLine, int nCmdShow) {
	for (int i = 0; i < WIDTH * HEIGHT * 3; i++) {
		buffer[i] = 255;
	}
	IJWorld world;
	world.camera.position = IJVector(0, 0, 1, 0);
	world.camera.upwards = IJAuxVector(0, 1, 0);
	world.camera.direction = IJAuxVector(0, 0, -1);
	world.camera.type = IJ_ORTHOGRAPHIC;
	IJShape cube;
	cube.data[0] = IJVector(1.0, 0.0, 0.0, 0.0);
	cube.data[1] = IJVector(0.0, 1.0, 1.0, 0.0);
	cube.type = IJ_CUBE;
	world.shapes.push_back(cube);
	IJPatch *patch = VertexShaderStage1(world);
	patch = VertexShaderStage2(world, patch);
	//RasterizationStage1(world, patch);

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

	hwnd = CreateWindow(szAppName, TEXT("BitBlt Demo"),
		WS_OVERLAPPEDWINDOW,
		CW_USEDEFAULT, CW_USEDEFAULT,
		WIDTH + 20, HEIGHT + 40,
		NULL, NULL, hInstance, NULL);

	ShowWindow(hwnd, nCmdShow);
	UpdateWindow(hwnd);

	static HDC screen_hdc;
	static HWND screen_hwnd;
	static HDC hCompatibleDC; //兼容HDC  
	static HBITMAP hCompatibleBitmap; //兼容BITMAP  
	static HBITMAP hOldBitmap; //旧的BITMAP                   
	static BITMAPINFO binfo; //BITMAPINFO结构体  
    Line(IJVector(-0.5,-1, 0, 0), IJVector(0.7,0.3, 0, 0));
	ZeroMemory(&binfo, sizeof(BITMAPINFO));
	binfo.bmiHeader.biBitCount = 24;      //每个像素多少位，也可直接写24(RGB)或者32(RGBA)  
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
	//将颜色数据打印到屏幕上，这下面两个函数每帧都得调用  
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

	case WM_DESTROY:
		PostQuitMessage(0);
		return 0;
	}
	return DefWindowProc(hwnd, message, wParam, lParam);
}