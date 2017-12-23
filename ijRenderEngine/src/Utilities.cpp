/*
*                  The main pipeline functions of the engine.
*
*                       Created by 78ij in 2017.11
*/

#include "Utilities.h"

extern BYTE buffer[WIDTH * HEIGHT * 3];
extern float depthbuffer[WIDTH * HEIGHT];

bool compare(IJVector a, IJVector b) {
	return a[1] < b[1];
}

void SortColor(IJTriangle *triangle) {
	if (compare(triangle->data[0], triangle->data[1])) {
		swap((triangle->color[0]), (triangle->color[1]));
		swap(triangle->data[0], triangle->data[1]);
	}
	if (compare(triangle->data[1], triangle->data[2])) {
		swap(triangle->color[1], triangle->color[2]);
		swap(triangle->data[1], triangle->data[2]);
	}
	if (compare(triangle->data[0], triangle->data[1])) {
		swap(triangle->color[0], triangle->color[1]);
		swap(triangle->data[0], triangle->data[1]);
	}
	if (compare(triangle->data[1], triangle->data[2])) {
		swap(triangle->color[1], triangle->color[2]);
		swap(triangle->data[1], triangle->data[2]);
	}
}

void DrawOneLine(int x1, int y1, int x2, int y2, IJColor color1, IJColor color2, float zbuffer) {

	if (x2 < x1) {
		swap(x2, x1);
		swap(y2, y1);
		swap(color2, color1);
	}
	int x, y;
	x = x1;
	y = y1;
	// pallel to the axis
	if (y1 == y2) {
		if (x1 == x2) return;
		//pallel to x
		while (x <= x2) {
			x++;
			if (!isexceed(x, y)) {
			if (depthbuffer[y * WIDTH + x] < zbuffer) {
					DrawPoint(x, y, interpolate<IJColor>(color1, color2, x1, x2, x));
					depthbuffer[y * WIDTH + x] = zbuffer;
				}
			}
		}
		return;
	}
	if (x1 == x2) {

		//pallel to y
		if (y1 > y2) {
			swap(x2, x1);
			swap(y2, y1);
			swap(color2, color1);
			y = y1;
		}
		while (y < y2) {
			y++;
			if (!isexceed(x, y)) {
			if (depthbuffer[y * WIDTH + x] < zbuffer) {
					DrawPoint(x, y, interpolate<IJColor>(color1, color2, y1, y2, y));
					depthbuffer[y * WIDTH + x] = zbuffer;
				}
			}
		}
		return;
	}
	int dx = x2 - x1, dy = y2 - y1;
	int p;
	int twoDy = 2 * dy, twoMinusDx = 2 * (dy - dx), twoDx = 2 * dx, twoMinusDy = 2 * (dx - dy);
	int twoSum = 2 * (dy + dx);
	double k = (double)dy / (double)dx;

	//0<k<1
	if (k<1.0&&k>0.0) {
		p = 2 * dy - dx;
		while (x < x2) {
			x++;
			if (p < 0)
				p += twoDy;
			else {
				y++;
				p += twoMinusDx;
			}
			if (!isexceed(x, y)) {
			if (depthbuffer[y * WIDTH + x] < zbuffer) {
					DrawPoint(x, y, interpolate<IJColor>(color1, color2, x1, x2, x));
					depthbuffer[y * WIDTH + x] = zbuffer;
				}
			}
		}
	}
	//k>=1 
	if (k >= 1.0) {
		p = dy;
		while (y < y2) {
			y++;
			if (p < 0)
				p += twoDx;
			else {
				x++;
				p += twoMinusDy;
			}
			if (!isexceed(x, y)) {
			if (depthbuffer[y * WIDTH + x] < zbuffer) {
				
					DrawPoint(x, y, interpolate<IJColor>(color1, color2, y1, y2, y));
					depthbuffer[y * WIDTH + x] = zbuffer;
				}
			}
		}
	}
	//0>k>-1
	if (k > -1 && k < 0) {
		p = 2 * dy + dx;
		while (x < x2) {
			x++;
			if (p >= 0)
				p += twoDy;
			else {
				y--;
				p += twoSum;
			}
			if (!isexceed(x, y)) {
			if (depthbuffer[y * WIDTH + x] < zbuffer) {
					DrawPoint(x, y, interpolate<IJColor>(color1, color2, x1, x2, x));
					depthbuffer[y * WIDTH + x] = zbuffer;
				}
			}
		}
	}
	//k<-1
	if (k <= -1) {
		p = 2 * dx - dy;
		while (y > y2) {
			y--;
			if (p >= 0)
				p -= twoDx;
			else {
				x++;
				p -= twoSum;
			}
			if (!isexceed(x, y)) {
			if (depthbuffer[y * WIDTH + x] < zbuffer) {
					DrawPoint(x, y, interpolate<IJColor>(color1, color2, y1, y2, y));
					depthbuffer[y * WIDTH + x] = zbuffer;
				}
			}
		}
	}
}

template<class T>
T &interpolate(T start, T end, int v1, int v2,int mid) {
	double portion = (mid - v1) / (v2 - v1);
	T ret = start + (end - start) * portion;
	return ret;
}

void DrawPoint(int x, int y, IJColor color) {
	buffer[(y * WIDTH + x) * 3] = color.b;
	buffer[(y * WIDTH + x) * 3 + 1] = color.g;
	buffer[(y * WIDTH + x) * 3 + 2] = color.r;
}

bool isexceed(int x, int y) {
	if (((y * WIDTH + x) * 3 + 2) > (WIDTH * HEIGHT * 3 - 1))  return true;
	else if ((y * WIDTH + x) * 3 + 2 < 0) return true;
	else if (y >= WIDTH - 1 || x >= WIDTH - 1)  return true;
	else if (y < 0 || x < 0)  return true;
	return false;
}

void DrawFlatBottomTriangle(IJVector a, IJVector b, IJVector c, IJColor color1, IJColor color2, IJColor color3, float zbuffer)
{
	a = a * 399.5 + IJVector(399.5, 399.5, 0, 0);
	b = b * 399.5 + IJVector(399.5, 399.5, 0, 0);
	c = c * 399.5 + IJVector(399.5, 399.5, 0, 0);
	int x1 = a[0];
	int y1 = a[1];
	int x2 = b[0];
	int y2 = b[1];
	int x3 = c[0];
	int y3 = c[1];
	DrawOneLine(x1, y1, x2, y2, color2, color1, zbuffer);
	DrawOneLine(x1, y1, x3, y3, color1, color3, zbuffer);
	DrawOneLine(x2, y2, x3, y3, color2, color3, zbuffer);
	for (int y = y1; y > y2; --y)
	{
		int xs, xe;
		IJColor colors = interpolate<IJColor>(color1, color2, y1, y2, y);
		IJColor colore = interpolate<IJColor>(color1, color3, y1, y3, y);
		xs = (y1 - y) * (x2 - x1) / (y1 - y2) + x1;
		xe = (y1 - y) * (x3 - x1) / (y1 - y3) + x1;
		DrawOneLine(xs, y, xe, y, colors, colore, zbuffer);
	}
}

void DrawFlatTopTriangle(IJVector a, IJVector b, IJVector c, IJColor color1, IJColor color2, IJColor color3, float zbuffer)
{
	a = a * 399.5 + IJVector(399.5, 399.5, 0, 0);
	b = b * 399.5 + IJVector(399.5, 399.5, 0, 0);
	c = c * 399.5 + IJVector(399.5, 399.5, 0, 0);
	int x1 = a[0];
	int y1 = a[1];
	int x2 = b[0];
	int y2 = b[1];
	int x3 = c[0];
	int y3 = c[1];
	DrawOneLine(x1, y1, x3, y3, color1, color3, zbuffer);
	DrawOneLine(x2, y2, x3, y3, color3, color2, zbuffer);
	DrawOneLine(x2, y2, x3, y3, color2, color3, zbuffer);
	for (int y = y1; y > y3; --y)
	{
		IJColor colors = interpolate<IJColor>(color1, color3, y1, y3, y);
		IJColor colore = interpolate<IJColor>(color2, color3, y2, y3, y);
		int xs, xe;
		xs = (y1 - y) * (x3 - x1) / (y1 - y3) + x1;
		xe = (y2 - y) * (x3 - x2) / (y2 - y3) + x2;
		DrawOneLine(xs, y, xe, y, colors, colore, zbuffer);
	}
}

void TriangleRasterization(IJTriangle *triangle) {
	SortColor(triangle);
	double xmiddle, deltax = (triangle->data[2][0]) - (triangle->data[0][0]);
	double portion = ((triangle->data[0][1]) - (triangle->data[1][1]))
		/ ((triangle->data[0][1]) - (triangle->data[2][1]));
	double portionx = portion * deltax;
	xmiddle = triangle->data[0][0] + portionx;
	IJColor middlecolor = (triangle->color[2] - triangle->color[0]) * portion + triangle->color[0];
	//middlecolor.normalize();
	if (xmiddle >= triangle->data[1][0]) {
		DrawFlatBottomTriangle(
			IJVector(triangle->data[0][0], triangle->data[0][1], 0, 1),
			IJVector(triangle->data[1][0], triangle->data[1][1], 0, 1),
			IJVector(xmiddle, triangle->data[1][1], 0, 1),
			triangle->color[0],triangle->color[1], middlecolor, triangle->zbuffer);
		DrawFlatTopTriangle(
			IJVector(triangle->data[1][0], triangle->data[1][1], 0, 1),
			IJVector(xmiddle, triangle->data[1][1], 0, 1),
			IJVector(triangle->data[2][0], triangle->data[2][1], 0, 1),
			triangle->color[1], middlecolor, triangle->color[2], triangle->zbuffer);
	}
	else {
		DrawFlatBottomTriangle(
			IJVector(triangle->data[0][0], triangle->data[0][1], 0, 1),
			IJVector(xmiddle, triangle->data[1][1], 0, 1),
			IJVector(triangle->data[1][0], triangle->data[1][1], 0, 1),
			triangle->color[0], middlecolor, triangle->color[1], triangle->zbuffer);
		DrawFlatTopTriangle(IJVector(xmiddle, triangle->data[1][1], 0, 1),
			IJVector(triangle->data[1][0], triangle->data[1][1], 0, 1),
			IJVector(triangle->data[2][0], triangle->data[2][1], 0, 1),
			middlecolor, triangle->color[1], triangle->color[2], triangle->zbuffer);
	}
}

void FreePatch(IJPatch *data, IJWorld world) {
	IJuint size = world.shapes.size();
	for (int shapecount = 0; shapecount < size; shapecount++) {
		IJuint patchsize = (data + shapecount)->size;
		free((data + shapecount)->data);
	}
	free(data);
}

IJVector getPoint(double u, double v, IJVector center, double radius) {

	double x = radius*sin(PI*v)*cos(PI2*u);
	double y = radius*sin(PI*v)*sin(PI2*u);
	double z = radius*cos(PI*v);
	return IJVector(center[0] + x,
		center[1] + y,
		center[2] + z,
		1.0);
}

void DividePolygon(IJTriangle *triangle, IJVector a, IJVector b, IJVector c, IJVector d) {
	triangle->data[0] = a;
	triangle->data[1] = b;
	triangle->data[2] = c;
	(triangle + 1)->data[0] = a;
	(triangle + 1)->data[1] = c;
	(triangle + 1)->data[2] = d;
}

void CalculateCubeVertices(IJShape cube, IJVector *vertices) {
	IJVector vertex1 = cube.data[0];
	IJVector vertex2 = cube.data[1];
	IJVector tempvertices[8] = {
		IJVector(vertex1[0],vertex1[1],vertex1[2],1),
		IJVector(vertex1[0],vertex2[1],vertex1[2],1),
		IJVector(vertex2[0],vertex2[1],vertex1[2],1),
		IJVector(vertex2[0],vertex1[1],vertex1[2],1),
		IJVector(vertex1[0],vertex1[1],vertex2[2],1),
		IJVector(vertex1[0],vertex2[1],vertex2[2],1),
		IJVector(vertex2[0],vertex2[1],vertex2[2],1),
		IJVector(vertex2[0],vertex1[1],vertex2[2],1),
	};
	for (int i = 0; i < 8; i++) {
		vertices[i] = tempvertices[i];
	}
}


void BlinnPhong(IJWorld world, IJVector point, IJVector normal, IJColor *color) {
	double brightness;
	IJAuxVector norm(normal[0], normal[1], normal[2]);
	norm.normalize();
	IJVector direction = world.light.position - point;
	IJAuxVector dir(direction[0], direction[1], direction[2]);
	dir.normalize();
	IJVector viewposition = world.camera.position - point;
	IJAuxVector view(viewposition[0], viewposition[1], viewposition[2]);
	view.normalize();
	IJAuxVector half = (dir + view);
	half.normalize();
	double temp = half.dot(norm);
	brightness = 0.5 * max(dir.dot(norm), 0.0) + 0.5 * std::pow(max(temp, 0.0), 8);
	BYTE coloroffset = brightness * 255;
	*color = *color + coloroffset;
}
 
IJVector Calculatenormal(IJTriangle triangle) {
	IJAuxVector vector1(
		triangle.data[0][0] - triangle.data[1][0],
		triangle.data[0][1] - triangle.data[1][1],
		triangle.data[0][2] - triangle.data[1][2]);
	IJAuxVector vector2(
		triangle.data[1][0] - triangle.data[2][0],
		triangle.data[1][1] - triangle.data[2][1],
		triangle.data[1][2] - triangle.data[2][2]);
	IJAuxVector cross = vector1.cross(vector2);
	return IJVector(
		cross[0],
		cross[1],
		cross[2],
		0);
}

void Processply(IJObject *object,IJWorld world) {
	vector<string> data;
	string current;
	fstream plystream;
	int vertexcount;
	int facecount;
	IJVector *points;
	IJTriangle *faces;
	plystream.open(object->path, fstream::in);
	getline(plystream, current);
	if (current != "ply")
		throw std::runtime_error("the file format must be ply!");
	getline(plystream, current);
	split(current, data, " ");
	if(data[1] != "ascii")
		throw std::runtime_error("the file encoding must be ascii!");
	while (getline(plystream, current)) {
		if (current == "end_header") break;
		else {
			split(current, data, " ");
			string token = data[0];
			if (token == "comment") continue;
			if (token == "element") {
				if (data[1] == "vertex"){
					string number = data[2];
					vertexcount = static_cast<int>(strtol(number.c_str(), nullptr, 10));
				}
				if (data[1] == "face") {
					string number = data[2];
					facecount = static_cast<int>(strtol(number.c_str(), nullptr, 10));
				}
			}
		}
	}
	points = (IJVector *)malloc(sizeof(IJVector) * vertexcount);
	for (int i = 0; i < vertexcount; i++) {
		vector<string> data;
		string current;
		IJVector *point1 = points + i;
		getline(plystream, current);
		split(current, data, " ");
		(*point1)[0] = static_cast<double>(strtof(data[0].c_str(), nullptr)) * 10;
		(*point1)[1] = static_cast<double>(strtof(data[1].c_str(), nullptr)) * 10;
		(*point1)[2] = static_cast<double>(strtof(data[2].c_str(), nullptr)) * 10;
		(*point1)[3] = 1;
	}
	faces = (IJTriangle *)malloc(sizeof(IJTriangle) * facecount);
#pragma omp parallel for
	for (int i = 0; i < facecount; i++) {
		vector<string> data;
		string current;
		IJTriangle *face1 = faces + i;
		getline(plystream, current);
		split(current, data, " ");
		int first = static_cast<int>(strtol(data[1].c_str(), nullptr, 10));
		int second = static_cast<int>(strtol(data[2].c_str(), nullptr, 10));
		int third = static_cast<int>(strtol(data[3].c_str(), nullptr, 10));
		(*face1).data[0] = points[first];
		(*face1).data[1] = points[second];
		(*face1).data[2] = points[third];
		for (int i = 0; i < 3; i++) {
			(*face1).color[i] = IJColor(0,0,0);
		}
		BlinnPhong(world, points[first], Calculatenormal(*face1),  &(*face1).color[0]);
		BlinnPhong(world, points[second], Calculatenormal(*face1), &(*face1).color[1]);
		BlinnPhong(world, points[second], Calculatenormal(*face1), &(*face1).color[2]);
	}
	plystream.close();
	free(points);
	object->size = facecount;
	object->triangles = faces;
}