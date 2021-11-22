#include "vector3.h"

Vector3::Vector3(float x0, float y0, float z0) :x(x0), y(y0), z(z0) {}
Vector3::Vector3() :x(0.0), y(0.0), z(0.0) {}
Vector3::~Vector3() {}

Vector3 Vector3::operator+(Vector3& vec) {
	float x1, y1, z1;
	x1 = this->x + vec.x;
	y1 = this->y + vec.y;
	z1 = this->z + vec.z;

	Vector3 temp(x1, y1, z1);
	return temp;
}

Vector3 Vector3::operator-(Vector3& vec) {
	float x1, y1, z1;
	x1 = this->x - vec.x;
	y1 = this->y - vec.y;
	z1 = this->z - vec.z;

	Vector3 temp(x1, y1, z1);
	return temp;
}

Vector3 Vector3::operator*(float num) {
	float x1, y1, z1;
	x1 = this->x * num;
	y1 = this->y * num;
	z1 = this->z * num;

	Vector3 temp(x1, y1, z1);
	return temp;
}

Vector3 Vector3::operator/(float num) {
	float x1, y1, z1;
	x1 = this->x / num;
	y1 = this->y / num;
	z1 = this->z / num;

	Vector3 temp(x1, y1, z1);
	return temp;
}

bool Vector3::operator==(Vector3& vec) {
	if ((fabs(this->x - vec.x) <= 1e-6) && (fabs(this->y - vec.y) <= 1e-6) && (fabs(this->z - vec.z) <= 1e-6))return true;
	return false;
}

float Vector3::dot(Vector3& vec) {
	auto result = this->x * vec.x + this->y * vec.y + this->z * vec.z;
	return result;
}

Vector3 Vector3::cross(Vector3& vec) {
	float a = this->y * vec.z - this->z * vec.y;
	float b = this->z * vec.x - vec.z * this->x;
	float c = this->x * vec.y - vec.x * this->y;
	Vector3 result(a, b, c);
	return  result;
}

Vector3 Vector3::normalize() {
	float length = sqrt(this->x * this->x + this->y * this->y + this->z * this->z);
	this->x /= length;
	this->y /= length;
	this->z /= length;

	return *this;
}




Ray::Ray() :o(Vector3()), d(Vector3()), k(0) {}
Ray::~Ray() {}
Ray::Ray(Vector3 os, Vector3 ds, float ks) :o(os), d(ds), k(ks) {}
bool Ray::IsInsert(Ray& r) {
	auto P = r.o - this->o;
	auto Q1 = this->d;
	auto Q2 = r.d;
	auto R = this->d.cross(r.d);
	glm::mat3x3 mat_t = {
		{P.x,P.y,P.z},
		{Q2.x,Q2.y,Q2.z},
		{R.x,R.y,R.z}
	};

	glm::mat3x3 mat_s = {
		{P.x,P.y,P.z},
		{Q1.x,Q1.y,Q1.z},
		{R.x,R.y,R.z}
	};

	auto Dt = glm::determinant(mat_t);
	auto Ds = glm::determinant(mat_s);
	if (Dt >= 0.0 && Ds >= 0.0) {

		return true;
	}
	return false;

}

