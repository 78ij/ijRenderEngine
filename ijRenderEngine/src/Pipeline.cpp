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
		center[3]);
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

void Line(IJVector a, IJVector b) {
	a = a * 399.5 + IJVector(399.5, 399.5, 0, 0);
	b = b * 399.5 + IJVector(399.5, 399.5, 0, 0);

	double m = (b[1] - a[1]) / (b[0] - a[0]);
	IJVector origin = a[0] > b[0] ? b : a;
	IJVector destination = a[0] < b[0] ? b : a;

	if (abs(m) <= 1) {
		double start = min(a[0], b[0]);
		double finish = max(a[0], b[0]);
		double y0 = a[0] > b[0] ? b[1] : a[1];
		double y1 = a[0] <= b[0] ? b[1] : a[1];
		bool isbigger = ((y1 - y0) >= 0);
		int st = start;
		int y = y0;
		for (; start <= finish; start++) {
			st = start;
			y = y0;
			if (isbigger) {
				if (!Linearithmatic(origin, destination, start + 1, y0 + 0.5)) {
					y0++;
				}
			}
			if (!isbigger) {
				if (Linearithmatic(origin, destination, start + 1, y0 + 0.5)) {
					y0--;

				}
			}
			if (((y * WIDTH + st) * 3 + 2) >(WIDTH * HEIGHT * 3 - 1)) continue;
			if ((y * WIDTH + st) * 3 + 2 < 0) continue;
			if (y >= WIDTH - 1 || st >= WIDTH - 1) continue;
			if (y < 0 || st < 0) continue;
			buffer[(y * WIDTH + st) * 3] = 0;
			buffer[(y * WIDTH + st) * 3 + 1] = 0;
			buffer[(y * WIDTH + st) * 3 + 2] = 0;
		}
	}
	else {
		double y0 = min(a[1], b[1]);
		double y1 = max(a[1], b[1]);
		double start = a[1] > b[1] ? b[0] : a[0];
		double finish = a[1] <= b[1] ? b[0] : a[0];
		bool isbigger = ((finish - start) >= 0);
		for (; y0 <= y1; y0++) {
			int st = start;
			int y = y0;

			if (isbigger) {
				if (Linearithmatic(origin, destination, start + 1, y0 + 0.5)) {
					start++;
				}
			}
			if (!isbigger) {
				if (!Linearithmatic(origin, destination, start + 1, y0 + 0.5)) {
					start--;
				}
			}
			if ((y * WIDTH + st) * 3 + 2 > WIDTH * HEIGHT * 3 - 1) continue;
			if ((y * WIDTH + st) * 3 + 2 < 0)continue;
			if (y > WIDTH - 1 || st > WIDTH - 1) continue;
			if (y < 0 || st < 0) continue;
			buffer[(y * WIDTH + st) * 3] = 0;
			buffer[(y * WIDTH + st) * 3 + 1] = 0;
			buffer[(y * WIDTH + st) * 3 + 2] = 0;
		}
	}
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
			double ustepinterval = 1 / (double)ustep;
			double vstepinterval = 1 / (double)vstep;
			// the triangles at bottom
			double u = 0, v = 0;
			for (int i = 0; i < ustep; i++) {
				IJTriangle *triangle = new(triangles + i) IJTriangle;
				triangle->data[0] = getPoint(0, 0, tempshape.data[0], tempshape.radius);
				triangle->data[1] = getPoint(u, vstepinterval, tempshape.data[0], tempshape.radius);
				triangle->data[2] = getPoint(u + ustepinterval, vstepinterval, tempshape.data[0], tempshape.radius);
			}
			// the polygons of in the middle
			u = 0, v = vstep;
			IJTriangle * triangle1 = new(triangles + ustep) IJTriangle;
			IJTriangle * triangle2 = new(triangles + ustep + 1) IJTriangle;
			for (int i = 1; i < vstep - 1; i++) {
				for (int j = 0; j < ustep; j++) {
					DividePolygon(triangle1, getPoint(u, v, tempshape.data[0], tempshape.radius),
						getPoint(u + ustepinterval, v, tempshape.data[0], tempshape.radius),
						getPoint(u + ustepinterval, v + vstepinterval, tempshape.data[0], tempshape.radius),
						getPoint(u, v + vstepinterval, tempshape.data[0], tempshape.radius));
					u += ustepinterval;
					triangle1 = new(triangles + i + j) IJTriangle;
					triangle2 = new(triangles + i + j + 1) IJTriangle;
				}
				v += vstepinterval;
		    }
			// the polygons on the top
			u = 0;
			for (int i = 0; i < ustep; i++) {
				IJTriangle *triangle = new(triangles + ((vstep - 1) * ustep) * 2 - ustep + i) IJTriangle;
				triangle->data[0] = getPoint(0, 1, tempshape.data[0], tempshape.radius);
				triangle->data[1] = getPoint(u, 1 - vstepinterval, tempshape.data[0], tempshape.radius);
				triangle->data[2] = getPoint(u + ustepinterval, 1 - vstepinterval, tempshape.data[0], tempshape.radius);
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


IJPatch *RasterizationStage1(IJWorld world, IJPatch *data) {
	IJuint size = world.shapes.size();
	for (int shapecount = 0; shapecount < size; shapecount++) {
		IJuint patchsize = (data + shapecount)->size;
		for (int primitivecount = 0; primitivecount < patchsize; primitivecount++) {
			Line(((data + shapecount)->data + primitivecount)->data[0],
				((data + shapecount)->data + primitivecount)->data[1]);
			Line(((data + shapecount)->data + primitivecount)->data[1],
				((data + shapecount)->data + primitivecount)->data[2]);
			Line(((data + shapecount)->data + primitivecount)->data[2],
				((data + shapecount)->data + primitivecount)->data[0]);
		}
	}
	return NULL;
}
