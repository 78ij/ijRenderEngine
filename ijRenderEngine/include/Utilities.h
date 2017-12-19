/*
*				Auxiliry functions used to calculate vertices,etc
*
*                       Created by 78ij in 2017.12
*/


#ifndef UTILITIES_H
#define UTILITIES_H

#include "Bases.h"

bool compare(IJVector a, IJVector b);
void SortColor(IJTriangle * triangles);
void DrawOneLine(int x1, int y1, int x2, int y2, IJColor *color1, IJColor *color2, IJint zbuffer);
void DrawFlatBottomTriangle(IJVector a, IJVector b, IJVector c, 
	IJColor *color1, IJColor *color2, IJColor * color3, IJint zbuffer);
void DrawFlatTopTriangle(IJVector a, IJVector b, IJVector c,
	IJColor *color1, IJColor *color2, IJColor * color3, IJint zbuffer);
void TriangleRasterization(IJTriangle *triangle);
void FreePatch(IJPatch *data, IJWorld world);
IJVector getPoint(double u, double v, IJVector center, double radius);
void DividePolygon(IJTriangle *triangle, IJVector a, IJVector b, IJVector c, IJVector d);
void CalculateCubeVertices(IJShape cube, IJVector *vertices);
void Line(IJVector a, IJVector b, IJColor *color);
void BlinnPhong(IJWorld world, IJVector point, IJVector normal, IJColor *color);
void Processply(IJObject *object,IJWorld world);


#endif