/*
*                  The main pipeline functions of the engine.
*
*                       Created by 78ij in 2017.11
*/

#include "Pipeline.h"
#include "Utilities.h"

//global buffer
BYTE buffer[WIDTH * HEIGHT * 3];
float depthbuffer[WIDTH * HEIGHT];

//---------------------------------------------------
//Core pipeline functions
//---------------------------------------------------
IJPatch *VertexShaderStage1(IJWorld world) {
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
					triangle->color[j] = tempshape.color;
				}
			}
			CalculateCubeVertices(tempshape, vertices);
			DividePolygon(triangles, vertices[0], vertices[3], vertices[2], vertices[1]);
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
					triangle2->color[j] = tempshape.color;
				}
			}
			double ustepinterval = 1 / (double)ustep;
			double vstepinterval = 1 / (double)(vstep);
			// the triangles at the bottom
#pragma omp parallel for
			for (int i = 0; i < ustep; i++) {
				double u = ustepinterval * i, v = 0;
				IJTriangle *triangle = triangles + i;
				triangle->data[0] = getPoint(u + ustepinterval, vstepinterval, tempshape.data[0], tempshape.radius);
				triangle->data[1] = getPoint(0, 0, tempshape.data[0], tempshape.radius);
				triangle->data[2] = getPoint(u, vstepinterval, tempshape.data[0], tempshape.radius);
				BlinnPhong(world, triangle->data[0], triangle->data[0] - tempshape.data[0], &triangle->color[0]);
				BlinnPhong(world, triangle->data[1], triangle->data[1] - tempshape.data[0], &triangle->color[1]);
				BlinnPhong(world, triangle->data[2], triangle->data[2] - tempshape.data[0], &triangle->color[2]);
			}
			// the polygons in the middle
#pragma omp parallel for
			for (int i = 1; i < vstep - 1; i++) {
				for (int j = 0; j < ustep; j++) {
					double u = ustepinterval * j;
					double v = vstepinterval * i;
					IJTriangle * triangle1 = triangles + ustep + 2 * (j + (i - 1) * ustep);
					DividePolygon(triangle1, getPoint(u, v, tempshape.data[0], tempshape.radius),
						getPoint(u, v + vstepinterval, tempshape.data[0], tempshape.radius),
						getPoint(u + ustepinterval, v + vstepinterval, tempshape.data[0], tempshape.radius),
						getPoint(u + ustepinterval, v, tempshape.data[0], tempshape.radius));
					BlinnPhong(world, triangle1->data[0], triangle1->data[0] - tempshape.data[0], &triangle1->color[0]);
					BlinnPhong(world, triangle1->data[1], triangle1->data[1] - tempshape.data[0], &triangle1->color[1]);
					BlinnPhong(world, triangle1->data[2], triangle1->data[2] - tempshape.data[0], &triangle1->color[2]);
					BlinnPhong(world, (triangle1 + 1)->data[0], (triangle1 + 1)->data[0] - tempshape.data[0], &(triangle1 + 1)->color[0]);
					BlinnPhong(world, (triangle1 + 1)->data[1], (triangle1 + 1)->data[1] - tempshape.data[0], &(triangle1 + 1)->color[1]);
					BlinnPhong(world, (triangle1 + 1)->data[2], (triangle1 + 1)->data[2] - tempshape.data[0], &(triangle1 + 1)->color[2]);
				}
			}
			// the triangles on the top
#pragma omp parallel for
			for (int i = 0; i < ustep; i++) {
				double u = ustepinterval * i;
				IJTriangle *triangle = triangles + ((vstep - 1) * ustep) * 2 - ustep + i;
				triangle->data[0] = getPoint(0, 1, tempshape.data[0], tempshape.radius);
				triangle->data[1] = getPoint(u, 1 - vstepinterval, tempshape.data[0], tempshape.radius);
				triangle->data[2] = getPoint(u + ustepinterval, 1 - vstepinterval, tempshape.data[0], tempshape.radius);
				BlinnPhong(world, triangle->data[0], triangle->data[0] - tempshape.data[0], &triangle->color[0]);
				BlinnPhong(world, triangle->data[1], triangle->data[1] - tempshape.data[0], &triangle->color[1]);
				BlinnPhong(world, triangle->data[2], triangle->data[2] - tempshape.data[0], &triangle->color[2]);
			}
			IJPatch *patch = new(ret + shapecount) IJPatch;
			patch->data = triangles;
			patch->size = ((vstep - 1) * ustep) * 2;
			break;
		}
		case IJ_OBJECT: {
			IJPatch *patch = new(ret + shapecount) IJPatch;
			Processply(&tempshape.object, world);
			patch->data = tempshape.object.triangles;
			patch->size = tempshape.object.size;
		}
		}
	}
	return ret;
}

IJPatch *VertexShaderStage2(IJWorld world, IJPatch *data) {
	IJuint size = world.shapes.size();
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
	if (world.camera.type == IJ_PERSPECTIVE) {
		double n = world.camera.znear;
		double f = world.camera.zfar;
		double r = world.camera.fov;
		IJtransform perstransform;
		perstransform << (n / r), 0, 0, 0,
			0, (n / r), 0, 0,
			0, 0, -(f + n) / (f - n), -(2 * f * n) / (f - n),
			0, 0, -1, 0;
			transform = perstransform *  transform * auxtransform ;
	}
	else {
		transform = transform * auxtransform;

	}
	for (int shapecount = 0; shapecount < size; shapecount++) {
		IJuint patchsize = (data + shapecount)->size;
#pragma omp parallel for
		for (int primitivecount = 0; primitivecount < patchsize; primitivecount++) {
			((data + shapecount)->data + primitivecount)->data[0] = transform *
				((data + shapecount)->data + primitivecount)->data[0];
			((data + shapecount)->data + primitivecount)->data[1] = transform *
				((data + shapecount)->data + primitivecount)->data[1];
			((data + shapecount)->data + primitivecount)->data[2] = transform *
				((data + shapecount)->data + primitivecount)->data[2];

			double w =((data + shapecount)->data + primitivecount)->data[0][3];
			((data + shapecount)->data + primitivecount)->data[0] /= w;
			w = ((data + shapecount)->data + primitivecount)->data[1][3];
			((data + shapecount)->data + primitivecount)->data[1] /= w;
            w = ((data + shapecount)->data + primitivecount)->data[2][3];
			((data + shapecount)->data + primitivecount)->data[2] /= w;
			((data + shapecount)->data + primitivecount)->zbuffer =
				(((data + shapecount)->data + primitivecount)->data[0][2] +
				((data + shapecount)->data + primitivecount)->data[1][2] +
					((data + shapecount)->data + primitivecount)->data[2][2]);
		}
	}

	return data;
}

IJPatch *RasterizationStage1(IJWorld world, IJPatch *data) {
	IJuint size = world.shapes.size();
	for (int shapecount = 0; shapecount < size; shapecount++) {
		IJuint patchsize = (data + shapecount)->size;
#pragma omp parallel for
		for (int primitivecount = 0; primitivecount < patchsize; primitivecount++) {
			TriangleRasterization((data + shapecount)->data + primitivecount);
		}
		
	}
	return NULL;
}

