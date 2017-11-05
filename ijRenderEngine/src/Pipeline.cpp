#include"Pipeline.h"

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
			triangles = new(triangles) IJTriangle;
			CalculateCubeVertices(tempshape, vertices);
			DividePolygon(triangles,vertices[0], vertices[3], vertices[2], vertices[1]);
			DividePolygon(triangles + 2, vertices[0], vertices[4], vertices[5], vertices[1]);
			DividePolygon(triangles + 8, vertices[3], vertices[7], vertices[4], vertices[0]);
			DividePolygon(triangles + 6, vertices[3], vertices[7], vertices[6], vertices[2]);
			DividePolygon(triangles + 4, vertices[1], vertices[5], vertices[6], vertices[2]);
			DividePolygon(triangles + 10, vertices[4], vertices[7], vertices[6], vertices[5]);
			IJPatch *patch = new(ret + shapecount) IJPatch;
			patch->data = triangles;
			patch->size = 12;
			break;
		}
		case IJ_SPHERE: {
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