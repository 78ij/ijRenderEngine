
#ifndef BASE_H
#define BASE_H

#include "Dense"
#include <cmath>
#include <vector>
#include <new>
#include<Windows.h>
#include "StdVector"

//types of shape
#define IJ_CUBE                                             0
#define IJ_SPHERE                                           1
#define IJ_TRIANGLE                                         2
#define IJ_POLYGEN                                          3  //4 vertices

//types of camera
#define IJ_PERSPECTIVE                                      4
#define IJ_ORTHOGRAPHIC                                     5

//types of data
#define IJuint                                              unsigned long
#define IJint                                               long

//essential math data types
typedef Eigen::Matrix<double,4,1,Eigen::DontAlign>          IJVector;
typedef Eigen::Matrix<double, 3, 1, Eigen::DontAlign>       IJAuxVector;
typedef Eigen::Matrix <double,4,4,Eigen::RowMajor>          IJtransform;
typedef BYTE                                                IJColor;

// Returns 0 when the result is below 0
inline double ClampedCos(double angle) { 
	return std::cos(angle) < 0 ? 0 :
		std::cos(angle);
}  

struct IJShape {
	IJuint type;    // cube,sphere,rectangle
    IJVector data[4];
	IJColor  color[3];
	IJuint step;    // aviliable only for spheres.
};

struct IJTriangle {
	IJVector data[3];
	IJColor  color[3];
	double zbuffer;
};

struct IJPatch {
	IJTriangle *data;
	IJuint size;
};

struct IJCamera {
	IJVector position;
	IJAuxVector direction;
	IJAuxVector upwards;
	IJuint type;
};

struct IJWorld {
	std::vector<IJShape> shapes;
	IJCamera camera;
};

#endif

