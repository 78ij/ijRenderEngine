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
void DrawOneLine(int x1, int y1, int x2, int y2, IJColor color1, IJColor color2, float zbuffer);
void DrawFlatBottomTriangle(IJVector a, IJVector b, IJVector c, 
	IJColor color1, IJColor color2, IJColor  color3, float zbuffer);
void DrawFlatTopTriangle(IJVector a, IJVector b, IJVector c,
	IJColor color1, IJColor color2, IJColor color3, float zbuffer);
void TriangleRasterization(IJTriangle *triangle);
void FreePatch(IJPatch *data, IJWorld world);
IJVector getPoint(double u, double v, IJVector center, double radius);
void DividePolygon(IJTriangle *triangle, IJVector a, IJVector b, IJVector c, IJVector d);
void CalculateCubeVertices(IJShape cube, IJVector *vertices);
void BlinnPhong(IJWorld world, IJVector point, IJVector normal, IJColor *color);
void DrawPoint(int x, int y, IJColor color);
bool isexceed(int x, int y);
void Processply(IJObject *object,IJWorld world);
template<class T> T &interpolate(T start, T end, int v1, int v2, int mid);


#endif