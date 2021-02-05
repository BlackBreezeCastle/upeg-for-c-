#pragma once
#include<math.h>
#include<algorithm>

#define PI 3.14159265358979323846
double normalized_angle(double angle);

double normalized_rad(double rad);

//三维向量
class _declspec(dllexport) Vector3
{
private:
	double m_x = 0.0;
	double m_y = 0.0;
	double m_z = 0.0;

public:
	Vector3();

	Vector3(const double &x, const double& y, const double& z);

	double x()const;

	double y()const;

	double z()const;

	double magnitude()const;

	Vector3 operator +(const Vector3 &b)const;

	Vector3 operator -(const Vector3 &b)const;

	friend Vector3 operator *(const Vector3 &a ,const double &b);
	friend Vector3 operator *(const double &b, const Vector3 &a);

	Vector3 operator /(const double &b)const;

	bool operator==(const Vector3 &b)const;

	Vector3 normalized()const;

public:
	double static Dot(const Vector3&a,const Vector3&b);

	Vector3 static Cross(const Vector3&a, const Vector3&b);

	//向量夹角（弧度）
	double static Angle(const Vector3&a, const Vector3&b);
};

class _declspec(dllexport)Quaternion
{
private:
	double m_w;

	double m_x;
	double m_y;
	double m_z;
public:
	Quaternion();

	Quaternion(const double &w,const double &x, const double &y, const double &z);

	Quaternion(const Vector3 &pivot, const double &rad);

	Quaternion operator *(const Quaternion &b)const;

	Quaternion inverse()const;

	//右手系轴角顺时针旋转
	Vector3 rotate(const Vector3 &v)const;
};