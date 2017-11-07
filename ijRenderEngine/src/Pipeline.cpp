/*
*                  The main pipeline functions of the engine.
*
*                       Created by 78ij in 2017.11
*/

#include"Pipeline.h"

//global buffer
BYTE buffer[WIDTH * HEIGHT * 3];

//---------------------------------------------------
//Auxiliry functions used to calculate vertices,etc
//---------------------------------------------------
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

bool Linearithmatic(IJVector a, IJVector b, double c, double d) {
	if ((a[1] - b[1]) * c +
		(b[0] - a[0]) * d +
		a[0] * b[1] -
		b[0] * a[1] > 0) return true;
	else return false;
}
void Bresenhamline(int x1, int y1, int x2, int y2, IJColor *color) {
	int x, y, dx, dy, s1, s2, p, temp, interchange, i; 
	x = x1; 
	y = y1;   
	dx = abs(x2 - x1); 
	dy = abs(y2 - y1);  
	if (x2>x1)   s1 = 1; 
	else   s1 = -1;  
	if (y2>y1)   s2 = 1;
	else   s2 = -1;   
	if (dy>dx) { 
		temp = dx;  
		dx = dy;
		dy = temp;   
		interchange = 1; 
	}
	else  
		interchange = 0;  
	p = (dy << 1) - dx; 
	for (i = 1; i <= dx; i++) {
		if (((y * WIDTH + x) * 3 + 2) >(WIDTH * HEIGHT * 3 - 1)) continue;
		if ((y * WIDTH + x) * 3 + 2 < 0) continue;
		if (y >= WIDTH - 1 || x >= WIDTH - 1) continue;
		if (y < 0 || x < 0) continue;
		buffer[(y * WIDTH + x) * 3] = color[0];
		buffer[(y * WIDTH + x) * 3 + 1] = color[1];
		buffer[(y * WIDTH + x) * 3 + 2] = color[2];
		if (p >= 0) {
			if (interchange == 0) 
				y = y + s2; 
			else    
				x = x + s1;
			p = p -(dx << 1);
		}
		else { 
			if (interchange == 0) 
				x = x + s1; 
			else   
				y = y + s2;  
			p = p + (dy << 1);  
		}
	}
}
void Line(IJVector a, IJVector b,IJColor *color) {
	a = a * 399.5 + IJVector(399.5, 399.5, 0, 0);
	b = b * 399.5 + IJVector(399.5, 399.5, 0, 0);
	int x1 = a[0];
	int x2 = b[0];
	int y1 = a[1];
	int y2 = b[1];
	Bresenhamline(x1, y1, x2, y2, color);
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
		case IJ_TRIANGLE: {
			IJTriangle *triangle = (IJTriangle *)calloc(1,sizeof(IJTriangle));
			triangle = new(triangle) IJTriangle;
			for (int i = 0; i < 3; i++) {
				triangle->data[i] = tempshape.data[i];
				triangle->color[i] = tempshape.color[i];
			}
			IJPatch *patch = new(ret + shapecount) IJPatch;
			patch->data = triangle;
			patch->size = 1;
			break;
		}
		case IJ_CUBE: {
			IJVector vertices[8];
			IJTriangle *triangles = (IJTriangle *)calloc(12, sizeof(IJTriangle));
			for (int i = 0; i < 12; i++) {
				IJTriangle *triangle = new(triangles + i) IJTriangle;
				triangle->color[0] = tempshape.color[0];
				triangle->color[0] = tempshape.color[1];
				triangle->color[0] = tempshape.color[2];
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
				triangle2->color[0] = tempshape.color[0];
				triangle2->color[0] = tempshape.color[1];
				triangle2->color[0] = tempshape.color[2];
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
			    u += ustepinterval;
			}
			IJPatch *patch = new(ret + shapecount) IJPatch;
			patch->data = triangles;
			patch->size = ((vstep - 1) * ustep) * 2;
			break;
		}
		}
		return ret;
	}
}

IJPatch *RasterizationStage1(IJWorld world, IJPatch *data) {
	IJuint size = world.shapes.size();
	for (int shapecount = 0; shapecount < size; shapecount++) {
		IJuint patchsize = (data + shapecount)->size;
		for (int primitivecount = 0; primitivecount < patchsize; primitivecount++) {
			Line(((data + shapecount)->data + primitivecount)->data[0],
				((data + shapecount)->data + primitivecount)->data[1],
				((data + shapecount)->data + primitivecount)->color);
			Line(((data + shapecount)->data + primitivecount)->data[1],
				((data + shapecount)->data + primitivecount)->data[2],
				((data + shapecount)->data + primitivecount)->color);
			Line(((data + shapecount)->data + primitivecount)->data[2],
				((data + shapecount)->data + primitivecount)->data[0],
				((data + shapecount)->data + primitivecount)->color);
		}
	}
	return NULL;
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
					((data + shapecount)->data + primitivecount)->data[2][2];
			}
		}
	}
	return data;
}