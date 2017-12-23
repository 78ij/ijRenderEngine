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


template<class T>
inline T Clampmin(T value,T min) { 
	return value < min ? min : value;
}  

template<class T>
inline T Clampmax(T value, T max) {
	return value > max ? max : value;
}

template<class T>
inline T Clamp(T value, T max, T min) {
	T ret;
	ret = Clampmax(value, max);
	ret = Clampmin(ret, min);
	return ret;
} 

class IJColor {
public:
	IJColor() = default;
	BYTE r;
	BYTE g;
	BYTE b;
	IJColor(BYTE _r, BYTE _g, BYTE _b) {
		r = _r;
		g = _g;
		b = _b;
	}
	inline IJColor &operator+(const IJColor &rhs);
	inline IJColor &operator-(const IJColor &rhs);
	inline IJColor &operator*(const int &rhs);
	inline IJColor &operator*(const double &rhs);
	inline IJColor &operator/(const int &rhs);
	inline IJColor &operator+(const BYTE &rhs);
};

IJColor & IJColor::operator+(const IJColor &rhs) {
	r = Clamp<int>(rhs.r + r ,255, 0);
	g = Clamp<int>(rhs.r + g, 255, 0);
	b = Clamp<int>(rhs.r + b, 255, 0);
	return *this;
}

IJColor & IJColor::operator-(const IJColor &rhs) {
	r = Clamp<int>(r - rhs.r, 255, 0);
	g = Clamp<int>(g - rhs.g, 255, 0);
	b = Clamp<int>(b - rhs.b, 255, 0);
	return *this;
}

IJColor & IJColor::operator*(const int &rhs) {
	r = Clamp<int>(rhs * r, 255, 0);
	g = Clamp<int>(rhs * g, 255, 0);
	b = Clamp<int>(rhs * b, 255, 0);
	return *this;
}

IJColor & IJColor::operator*(const double &rhs) {
	r = Clamp<double>(rhs * r, 255, 0);
	g = Clamp<double>(rhs * g, 255, 0);
	b = Clamp<double>(rhs * b, 255, 0);
	return *this;
} 

IJColor & IJColor::operator/(const int &rhs) {
	r = Clamp<int>(r / rhs, 255, 0);
	g = Clamp<int>(g / rhs, 255, 0);
	b = Clamp<int>(b / rhs, 255, 0);
	return *this;
}

IJColor & IJColor::operator+(const BYTE &rhs) {
	r = Clamp<int>(rhs + r, 255, 0);
	g = Clamp<int>(rhs + g, 255, 0);
	b = Clamp<int>(rhs + b, 255, 0);
	return *this;
}


class IJTriangle {
public:
	IJTriangle() = default;
	IJVector data[3];
	IJColor  color[3];
	float zbuffer;
};

class IJObject {
public:
	IJObject() = default;
	IJTriangle *triangles;
	IJint size;
	string path;
};

class IJShape {
public:
	IJShape() = default;
	IJuint type;// cube,sphere,rectangle
	double radius; // aviliable only for spheres.
    IJVector data[4];
	IJColor  color;
	IJuint step[2];    // aviliable only for spheres.
	IJObject object;
};


class IJPatch {
public:
	IJPatch() = default;
	IJTriangle *data;
	IJuint size;
};


class IJCamera {
public:
	IJCamera() = default;
	IJVector position;
	IJAuxVector direction;
	IJAuxVector upwards;
	double znear;
	double zfar;
	double fov;
	IJuint type;
};

class IJLight {
public:
	IJLight() = default;
	IJVector position;
};

class IJWorld {
public:
	IJWorld() = default;
	std::vector<IJShape> shapes;
	IJCamera camera;
	IJLight light;
};


#endif