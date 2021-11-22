#ifndef VECTOR_H
#define VECTOR_H

#include "../externel/glm/glm.hpp"
#include "../externel/glm/gtc/matrix_transform.hpp"
#include "../externel/glm/gtc/type_ptr.hpp"
#include <iostream>
using namespace std;

class Vector3 {
public:
	Vector3(float x0, float y0, float z0);
	Vector3();
	~Vector3();

	Vector3 operator+(Vector3& vec);

	Vector3 operator-(Vector3& vec);

	Vector3 operator*(float num);

	Vector3 operator/(float num);

	bool operator==(Vector3& vec);

	float dot(Vector3& vec);

	Vector3 cross(Vector3& vec);
	
	Vector3 normalize();

public:
	float x;
	float y;
	float z;
};

class Ray {
public:
	Ray();
	~Ray();
	Ray(Vector3 os, Vector3 ds, float ks);
	bool IsInsert(Ray& r);

public:
	Vector3 o;
	Vector3 d;
	float k;

};

#endif