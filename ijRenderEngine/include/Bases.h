/*     
*      Basic definitions of constants and structures used in the engine
*                    
*                       Created by 78ij in 2017.11
*/

#ifndef BASE_H
#define BASE_H


#include <cmath>
#include <vector>
#include <new>
#include <iostream>
#include <cstring>
#include <fstream>
#include <algorithm>
#include <exception>
#include <omp.h>

#include "Eigen/Dense"
#include "Eigen/StdVector"
#include "pystring.h"

#ifdef _WIN32
#include <Windows.h>
#else
#include "svpng.h"
#endif


//global configuration of Width and Size and PI
#define WIDTH                                                 800
#define HEIGHT                                                800
#define PI                                                    3.14159265358979f
#define PI2                                                   6.28318530717958f

//types of shape
#define IJ_CUBE                                               0
#define IJ_SPHERE                                             1
#define IJ_TRIANGLE                                           2
#define IJ_POLYGEN                                            3  //4 vertices
#define IJ_OBJECT											  4

//types of camera
#define IJ_PERSPECTIVE                                        4
#define IJ_ORTHOGRAPHIC                                       5

//types of data
#define IJuint                                                unsigned long
#define IJint                                                 long

//essential math data types
typedef Eigen::Matrix<double,4,1,Eigen::DontAlign>            IJVector;
typedef Eigen::Matrix<double, 3, 1, Eigen::DontAlign>         IJAuxVector;
typedef Eigen::Matrix<double,4,4,Eigen::RowMajor>             IJtransform;
typedef BYTE                                                  IJColor;

using std::swap;
using std::max;
using std::cout;
using std::cin;
using std::endl;
using std::getline;
using std::vector;
using std::string;
using std::ifstream;
using std::fstream;
using std::exception;
using pystring::split; 

// Returns 0 when the result is below 0
inline double ClampedCos(double angle) { 
	return std::cos(angle) < 0 ? 0 :
		std::cos(angle);
}  

class IJTriangle {
public:
	IJVector data[3];
	IJColor  color[3][3];
	IJint zbuffer;
};

class IJObject {
public:
	IJTriangle *triangles;
	IJint size;
	string path;
};

class IJShape {
public:
	IJuint type;// cube,sphere,rectangle
	double radius; // aviliable only for spheres.
    IJVector data[4];
	IJColor  color[3];
	IJuint step[2];    // aviliable only for spheres.
	IJObject object;
};


class IJPatch {
public:
	IJTriangle *data;
	IJuint size;
};


class IJCamera {
public:
	IJVector position;
	IJAuxVector direction;
	IJAuxVector upwards;
	IJuint type;
};

class IJLight {
public:
	IJVector position;
};

class IJWorld {
public:
	std::vector<IJShape> shapes;
	IJCamera camera;
	IJLight light;
};


#endif

