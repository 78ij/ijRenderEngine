/*
*                  The main pipeline functions of the engine.
*
*                       Created by 78ij in 2017.11
*/

#include "Pipeline.h"

//global buffer
BYTE buffer[WIDTH * HEIGHT * 3];
IJint depthbuffer[WIDTH * HEIGHT];
using std::swap;
using std::max;
//---------------------------------------------------
//Auxiliry functions used to calculate vertices,etc
//---------------------------------------------------
bool compare(IJVector a, IJVector b) {
	return a[1] < b[1];
}
void SortColor(IJTriangle *triangle){
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

void DrawOneLine(int x1, int y1, int x2, int y2, IJColor *color1,IJColor *color2,IJint zbuffer)
{
	if (x2<x1)
	{
		swap(x2, x1);
		swap(y2, y1);
	}
	int x, y;
	x = x1;
	y = y1;
	// pallel to the axis
	if (y1 == y2)
	{
		//pallel to x
		while (x<x2)
		{
			x++;
			if (((y * WIDTH + x) * 3 + 2) >(WIDTH * HEIGHT * 3 - 1)) continue;
			if ((y * WIDTH + x) * 3 + 2 < 0) continue;
			if (y >= WIDTH - 1 || x >= WIDTH - 1) continue;
			if (y < 0 || x < 0) continue;
			if(depthbuffer[y * WIDTH + x] < zbuffer){
				buffer[(y * WIDTH + x) * 3] = color1[0] + (color2[0] - color1[0]) * (x - x1) / (x2 - x1);
				buffer[(y * WIDTH + x) * 3 + 1] = color1[1] + (color2[1] - color1[1]) * (x - x1) / (x2 - x1);
				buffer[(y * WIDTH + x) * 3 + 2] = color1[2] + (color2[2] - color1[2]) * (x - x1) / (x2 - x1);
				depthbuffer[y * WIDTH + x] = zbuffer;
			}
		}
		return;
	}
	if (x1 == x2)
	{
		if (y1 = y2)
			return;
		//pallel to y
		if (y1 > y2) {
			y = y2;
			y2 = y1;
		}
		while (y<y2)
		{
			y++;
			if (((y * WIDTH + x) * 3 + 2) >(WIDTH * HEIGHT * 3 - 1)) continue;
			if ((y * WIDTH + x) * 3 + 2 < 0) continue;
			if (y >= WIDTH - 1 || x >= WIDTH - 1) continue;
			if (y < 0 || x < 0) continue;
			if (depthbuffer[y * WIDTH] < zbuffer) {
				buffer[(y * WIDTH + x) * 3] = color1[0] + (color2[0] - color1[0]) * (y - y1) / (y2 - y1);
				buffer[(y * WIDTH + x) * 3 + 1] = color1[1] + (color2[1] - color1[1]) * (y - y1) / (y2 - y1);
				buffer[(y * WIDTH + x) * 3 + 2] = color1[2] + (color2[2] - color1[2]) * (y - y1) / (y2 - y1);
				depthbuffer[y * WIDTH + x] = zbuffer;
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
	if (k<1.0&&k>0.0)
	{
		p = 2 * dy - dx;
		while (x<x2)
		{
			x++;
			if (p<0)
				p += twoDy;
			else
			{
				y++;
				p += twoMinusDx;
			}
			if (((y * WIDTH + x) * 3 + 2) >(WIDTH * HEIGHT * 3 - 1)) continue;
			if ((y * WIDTH + x) * 3 + 2 < 0) continue;
			if (y >= WIDTH - 1 || x >= WIDTH - 1) continue;
			if (y < 0 || x < 0) continue;
			if (depthbuffer[y * WIDTH + x] < zbuffer) {
				buffer[(y * WIDTH + x) * 3] = color1[0] + (color2[0] - color1[0]) * (x - x1) / (x2 - x1);
				buffer[(y * WIDTH + x) * 3 + 1] = color1[1] + (color2[1] - color1[1]) * (x - x1) / (x2 - x1);
				buffer[(y * WIDTH + x) * 3 + 2] = color1[2] + (color2[2] - color1[2]) * (x - x1) / (x2 - x1);
				depthbuffer[y * WIDTH + x] = zbuffer;
			}

		}
	}
	//k>=1 
	if (k >= 1.0)
	{
		p = dy;
		while (y<y2)
		{
			y++;
			if (p<0)
				p += twoDx;
			else
			{
				x++;
				p += twoMinusDy;
			}
			if (((y * WIDTH + x) * 3 + 2) >(WIDTH * HEIGHT * 3 - 1)) continue;
			if ((y * WIDTH + x) * 3 + 2 < 0) continue;
			if (y >= WIDTH - 1 || x >= WIDTH - 1) continue;
			if (y < 0 || x < 0) continue;
			if (depthbuffer[y * WIDTH + x] < zbuffer) {
				buffer[(y * WIDTH + x) * 3] = color1[0] + (color2[0] - color1[0]) * (y - y1) / (y2 - y1);
				buffer[(y * WIDTH + x) * 3 + 1] = color1[1] + (color2[1] - color1[1]) * (y - y1) / (y2 - y1);
				buffer[(y * WIDTH + x) * 3 + 2] = color1[2] + (color2[2] - color1[2]) * (y - y1) / (y2 - y1);
				depthbuffer[y * WIDTH + x] = zbuffer;
			}
		}
	}
	//0>k>-1
	if (k>-1 && k<0)
	{
		p = 2 * dy + dx;
		while (x<x2)
		{
			x++;
			if (p >= 0)
				p += twoDy;
			else
			{
				y--;
				p += twoSum;
			}
			if (((y * WIDTH + x) * 3 + 2) >(WIDTH * HEIGHT * 3 - 1)) continue;
			if ((y * WIDTH + x) * 3 + 2 < 0) continue;
			if (y >= WIDTH - 1 || x >= WIDTH - 1) continue;
			if (y < 0 || x < 0) continue;
			if (depthbuffer[y * WIDTH + x] < zbuffer) {
				buffer[(y * WIDTH + x) * 3] = color1[0] + (color2[0] - color1[0]) * (x - x1) / (x2 - x1);
				buffer[(y * WIDTH + x) * 3 + 1] = color1[1] + (color2[1] - color1[1]) * (x - x1) / (x2 - x1);
				buffer[(y * WIDTH + x) * 3 + 2] = color1[2] + (color2[2] - color1[2]) * (x - x1) / (x2 - x1);
				depthbuffer[y * WIDTH + x] = zbuffer;
			}
		}
	}
	//k<-1
	if (k <= -1)
	{
		p = 2 * dx - dy;
		while (y>y2)
		{
			y--;
			if (p >= 0)
				p -= twoDx;
			else
			{
				x++;
				p -= twoSum;
			}
			if (((y * WIDTH + x) * 3 + 2) >(WIDTH * HEIGHT * 3 - 1)) continue;
			if ((y * WIDTH + x) * 3 + 2 < 0) continue;
			if (y >= WIDTH - 1 || x >= WIDTH - 1) continue;
			if (y < 0 || x < 0) continue;
			if (depthbuffer[y * WIDTH + x] < zbuffer) {
				buffer[(y * WIDTH + x) * 3] = color1[0] + (color2[0] - color1[0]) * (y - y1) / (y2 - y1);
				buffer[(y * WIDTH + x) * 3 + 1] = color1[1] + (color2[1] - color1[1]) * (y - y1) / (y2 - y1);
				buffer[(y * WIDTH + x) * 3 + 2] = color1[2] + (color2[2] - color1[2]) * (y - y1) / (y2 - y1);
				depthbuffer[y * WIDTH + x] = zbuffer;
			}
		}
	}
}

void DrawFlatBottomTriangle(IJVector a,IJVector b,IJVector c,IJColor *color1,IJColor *color2,IJColor * color3,IJint zbuffer)
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
	//DrawOneLine(x1, y1, x2, y2, color2, color1, zbuffer);
    //DrawOneLine(x1, y1, x3, y3, color1, color3, zbuffer);
	for (int y = y1; y > y2; --y)
	{
		int xs, xe;
		IJColor colors[3] = {
			(color1[0] - color2[0]) *  (y2 - y) / (y2 - y1) + color2[0],
			(color1[1] - color2[1]) *  (y2 - y) / (y2 - y1) + color2[1],
			(color1[2] - color2[2]) *  (y2 - y) / (y2 - y1) + color2[2]
		};
		IJColor colore[3] = {
			(color1[0] - color3[0]) *  (y2 - y) / (y2 - y1) + color3[0],
			(color1[1] - color3[1]) *  (y2 - y) / (y2 - y1) + color3[1],
			(color1[2] - color3[2]) *  (y2 - y) / (y2 - y1) + color3[2]
		};
		xs = (y1 - y) * (x2 - x1) / (y1 - y2) + x1;
		xe = (y1 - y) * (x3 - x1) / (y1 - y3) + x1;
		DrawOneLine(xs, y, xe, y, colors,colore,zbuffer);
	}
}
void DrawFlatTopTriangle(IJVector a, IJVector b, IJVector c, IJColor *color1, IJColor *color2, IJColor * color3,IJint zbuffer)
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
	//DrawOneLine(x1, y1, x3, y3, color1, color3, zbuffer);
	//DrawOneLine(x2, y2, x3, y3, color3, color2, zbuffer);
	for (int y = y1; y > y3; --y)
	{
		IJColor colors[3] = {
			(color1[0] - color3[0]) *  (y3 - y) / (y3 - y1) + color3[0],
			(color1[1] - color3[1]) *  (y3 - y) / (y3 - y1) + color3[1],
			(color1[2] - color3[2]) *  (y3 - y) / (y3 - y1) + color3[2]
		};
		IJColor colore[3] = {
			(color2[0] - color3[0]) *  (y3 - y) / (y3 - y1) + color3[0],
			(color2[1] - color3[1]) *  (y3 - y) / (y3 - y1) + color3[1],
			(color2[2] - color3[2]) *  (y3 - y) / (y3 - y1) + color3[2]
		};
		int xs, xe;
		xs = (y1 - y) * (x3 - x1) / (y1 - y3) + x1;
		xe = (y2 - y) * (x3 - x2) / (y2 - y3) + x2;
		DrawOneLine(xs, y, xe, y, colors,colore,zbuffer);
	}
}
void TriangleRasterization(IJTriangle *triangle) {
	SortColor(triangle);
	double xmiddle, deltax = (triangle->data[2][0]) - (triangle->data[0][0]);
	double portion = ((triangle->data[0][1]) - (triangle->data[1][1]))
		/ ((triangle->data[0][1]) - (triangle->data[2][1]));
	double portionx = portion * deltax;
	xmiddle = triangle->data[0][0] + portionx;
	IJColor middlecolor[3] = {
		(triangle->color[2][0] - triangle->color[0][0]) * portion + triangle->color[0][0],
		(triangle->color[2][1] - triangle->color[0][1]) * portion + triangle->color[0][1],
		(triangle->color[2][2] - triangle->color[0][2]) * portion + triangle->color[0][2]
	};

	if (xmiddle >= triangle->data[1][0]) {
		DrawFlatBottomTriangle(
			IJVector(triangle->data[0][0], triangle->data[0][1], 0, 1), 		
			IJVector(triangle->data[1][0], triangle->data[1][1], 0, 1), 
			IJVector(xmiddle, triangle->data[1][1], 0, 1),
			triangle->color[0], triangle->color[1], middlecolor,  triangle->zbuffer);
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
			middlecolor, triangle->color[1], triangle->color[2] , triangle->zbuffer);
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

IJVector getPoint(double u, double v, IJVector center,double radius){

	double x = radius*sin(PI*v)*cos(PI2*u);
	double y = radius*sin(PI*v)*sin(PI2*u);
	double z = radius*cos(PI*v);
	return IJVector(center[0] + x,
		center[1] + y,
		center[2] + z,
		1.0);
}

void DividePolygon(IJTriangle *triangle,IJVector a,IJVector b,IJVector c,IJVector d) {
	triangle->data[0] = a;
	triangle->data[1] = b;
	triangle->data[2] = c;
	(triangle + 1)->data[0] = a;
	(triangle + 1)->data[1] = c;
	(triangle + 1)->data[2] = d;
}

void CalculateCubeVertices(IJShape cube,IJVector *vertices) {
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
void Line(IJVector a, IJVector b,IJColor *color) {
	a = a * 399.5 + IJVector(399.5, 399.5, 0, 0);
	b = b * 399.5 + IJVector(399.5, 399.5, 0, 0);
	int x1 = a[0];
	int x2 = b[0];
	int y1 = a[1];
	int y2 = b[1];
	//Bresenhamline(x1, y1, x2, y2, color);
	//DrawOneLine(x1, y1, x2, y2, color,0);
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
	brightness = 0.5 * max(dir.dot(norm),0.0) + 0.5 * std::pow(max(temp, 0.0), 8);
	BYTE coloroffset = brightness * 255;
	color[0] = ((color[0] + coloroffset ) > 255) ? 255 : color[0] + coloroffset; //* 0.1;
	color[1] = ((color[1] + coloroffset ) > 255) ? 255 : color[1] + coloroffset; //* 0.6;
	color[2] = ((color[2] + coloroffset ) > 255) ? 255 : color[2] + coloroffset; //* 0.3;
}
//---------------------------------------------------
//Core pipeline functions
//---------------------------------------------------
IJPatch *VertexShaderStage1(IJWorld world){
	IJShape tempshape;
	IJuint size = world.shapes.size();
	IJPatch *ret = (IJPatch *)calloc(size, sizeof(IJPatch));
	if (!size) return NULL;
	for (int shapecount = 0; shapecount < size; shapecount++) {
		tempshape = world.shapes[shapecount];
		switch (tempshape.type) {
		case IJ_CUBE: {
			IJVector vertices[8];
			IJTriangle *triangles = (IJTriangle *)calloc(12, sizeof(IJTriangle));
			for (int i = 0; i < 12; i++) {
				IJTriangle *triangle = new(triangles + i) IJTriangle;
				for (int j = 0; j < 3; j++) {
					triangle->color[j][0] = tempshape.color[0];
					triangle->color[j][1] = tempshape.color[1];
					triangle->color[j][2] = tempshape.color[2];
				}
			}
			CalculateCubeVertices(tempshape, vertices);
			DividePolygon(triangles,vertices[0], vertices[3], vertices[2], vertices[1]);
			DividePolygon(triangles + 2, vertices[0], vertices[4], vertices[5], vertices[1]);
			DividePolygon(triangles + 4, vertices[0], vertices[4], vertices[7], vertices[3]);
			DividePolygon(triangles + 6, vertices[3], vertices[7], vertices[6], vertices[2]);
			DividePolygon(triangles + 8, vertices[1], vertices[5], vertices[6], vertices[2]);
			DividePolygon(triangles + 10, vertices[4], vertices[7], vertices[6], vertices[5]);
			IJPatch *patch = new(ret + shapecount) IJPatch;
			patch->data = triangles;
			patch->size = 12;
			break;
		}
		case IJ_SPHERE: {
			IJuint ustep = tempshape.step[0];
			IJuint vstep = tempshape.step[1];
			IJTriangle *triangles = (IJTriangle *)calloc(((vstep - 1) * ustep) * 2, sizeof(IJTriangle));
			for (int i = 0; i < ((vstep - 1) * ustep) * 2; i++) {
				IJTriangle *triangle2 = new(triangles + i) IJTriangle;
				for (int j = 0; j < 3; j++) {
					triangle2->color[j][0] = tempshape.color[0];
					triangle2->color[j][1] = tempshape.color[1];
					triangle2->color[j][2] = tempshape.color[2];
				}
			}
			double ustepinterval = 1 / (double)ustep;
			double vstepinterval = 1 / (double)(vstep -1 );
			// the triangles at bottom
			double u = 0, v = 0;
			for (int i = 0; i < ustep; i++) {
				IJTriangle *triangle = triangles + i;
				triangle->data[0] = getPoint(u + ustepinterval, vstepinterval, tempshape.data[0], tempshape.radius);
				triangle->data[1] = getPoint(0, 0, tempshape.data[0], tempshape.radius);
				triangle->data[2] = getPoint(u, vstepinterval, tempshape.data[0], tempshape.radius);
				BlinnPhong(world, triangle->data[0], triangle->data[0] - tempshape.data[0], triangle->color[0]);
				BlinnPhong(world, triangle->data[1], triangle->data[1] - tempshape.data[0], triangle->color[1]);
				BlinnPhong(world, triangle->data[2], triangle->data[2] - tempshape.data[0], triangle->color[2]);
				u += ustepinterval;
			}
			// the polygons of in the middle
			u = 0, v = vstep;
			IJTriangle * triangle1 = triangles + ustep;
			for (int i = 1; i < vstep - 1; i++) {
				for (int j = 0; j < ustep; j++) {
					DividePolygon(triangle1, getPoint(u, v, tempshape.data[0], tempshape.radius),
						getPoint(u, v + vstepinterval, tempshape.data[0], tempshape.radius),
						getPoint(u + ustepinterval, v + vstepinterval, tempshape.data[0], tempshape.radius),
						getPoint(u + ustepinterval, v, tempshape.data[0], tempshape.radius));
					BlinnPhong(world, triangle1->data[0], triangle1->data[0] - tempshape.data[0], triangle1->color[0]);
					BlinnPhong(world, triangle1->data[1], triangle1->data[1] - tempshape.data[0], triangle1->color[1]);
					BlinnPhong(world, triangle1->data[2], triangle1->data[2] - tempshape.data[0], triangle1->color[2]);
					BlinnPhong(world, (triangle1 + 1)->data[0], (triangle1 + 1)->data[0] - tempshape.data[0], (triangle1 + 1)->color[0]);
					BlinnPhong(world, (triangle1 + 1)->data[1], (triangle1 + 1)->data[1] - tempshape.data[0], (triangle1 + 1)->color[1]);
					BlinnPhong(world, (triangle1 + 1)->data[2], (triangle1 + 1)->data[2] - tempshape.data[0], (triangle1 + 1)->color[2]);
					u += ustepinterval;
					triangle1 += 2;
				}
				v += vstepinterval;
		    }
			// the triangles on the top
			u = 0;
			for (int i = 0; i < ustep; i++) {
				IJTriangle *triangle = triangles + ((vstep - 1) * ustep) * 2 - ustep + i;
				triangle->data[0] = getPoint(0, 1, tempshape.data[0], tempshape.radius);
				triangle->data[1] = getPoint(u, 1 - vstepinterval, tempshape.data[0], tempshape.radius);
				triangle->data[2] = getPoint(u + ustepinterval, 1 - vstepinterval, tempshape.data[0], tempshape.radius);
				BlinnPhong(world, triangle->data[0], triangle->data[0] - tempshape.data[0], triangle->color[0]);
				BlinnPhong(world, triangle->data[1], triangle->data[1] - tempshape.data[0], triangle->color[1]);
				BlinnPhong(world, triangle->data[2], triangle->data[2] - tempshape.data[0], triangle->color[2]);
			    u += ustepinterval;
			}
			IJPatch *patch = new(ret + shapecount) IJPatch;
			patch->data = triangles;
			patch->size = ((vstep - 1) * ustep) * 2;
			break;
		}
		}
	}
	return ret;
}

IJPatch *VertexShaderStage2(IJWorld world, IJPatch *data) {
	IJuint size = world.shapes.size();
	if (world.camera.type == IJ_ORTHOGRAPHIC) {
		IJtransform transform;
		IJtransform auxtransform;
		IJVector offset = world.camera.position;
		IJAuxVector w = -1 * world.camera.direction;
		w = w / w.norm();
		IJAuxVector u = world.camera.upwards.cross(w);
		u = u / u.norm();
		IJAuxVector v = w.cross(u);
		transform << u[0], u[1], u[2], 0,
			v[0], v[1], v[2], 0,
			w[0], w[1], w[2], 0,
			0, 0, 0, 1;
		auxtransform << 1, 0, 0, -offset[0],
			0, 1, 0, -offset[1],
			0, 0, 1, -offset[2],
			0, 0, 0, 1;
		transform = transform * auxtransform;
		for (int shapecount = 0; shapecount < size; shapecount++) {
			IJuint patchsize = (data + shapecount)->size;
			for (int primitivecount = 0; primitivecount < patchsize; primitivecount++) {
				((data + shapecount)->data + primitivecount)->data[0] = transform *
					((data + shapecount)->data + primitivecount)->data[0];
				((data + shapecount)->data + primitivecount)->data[1] = transform *
					((data + shapecount)->data + primitivecount)->data[1];
				((data + shapecount)->data + primitivecount)->data[2] = transform *
					((data + shapecount)->data + primitivecount)->data[2];
				((data + shapecount)->data + primitivecount)->zbuffer =
					((data + shapecount)->data + primitivecount)->data[0][2] +
					((data + shapecount)->data + primitivecount)->data[1][2] +
					((data + shapecount)->data + primitivecount)->data[2][2] * 200;
			}
		}
	}
	return data;
}

IJPatch *RasterizationStage1(IJWorld world, IJPatch *data) {
	IJuint size = world.shapes.size();
	for (int shapecount = 0; shapecount < size; shapecount++) {
		IJuint patchsize = (data + shapecount)->size;
		for (int primitivecount = 0; primitivecount < patchsize; primitivecount++) {
			TriangleRasterization((data + shapecount)->data + primitivecount);
		}
		
	}
	return NULL;
}

