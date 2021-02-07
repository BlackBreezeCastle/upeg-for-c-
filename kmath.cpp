#include "kmath.h"

double atanh(double x)
{
	return 0.5*log((1.0+x)/(1.0-x));
}


double normalized_angle(double angle)
{
	angle = angle + 180;
	angle = fmod(angle, 360);
	angle = angle - 180;
	return angle;
}
double normalized_rad(double rad)
{
	rad = rad + PI;
	rad = fmod(rad, (2 * PI));
	rad = rad - PI;
	return rad;
}

Vector3::Vector3(void)
{
	m_x = 0.0;
	m_y = 0.0;
	m_z = 0.0;
}

Vector3::Vector3(const double &x, const double &y, const double &z)
{
	m_x = x;
	m_y = y;
	m_z = z;
}

inline double Vector3::x() const
{
	return m_x;
}

inline double Vector3::y() const
{
	return m_y;
}

inline double Vector3::z() const
{
	return m_z;
}

inline double Vector3::magnitude() const
{
	return pow((m_x*m_x + m_y*m_y + m_z*m_z), 0.5);
}

inline Vector3 Vector3::operator +(const Vector3 &b)const
{
	return Vector3(m_x + b.m_x, m_y + b.m_y, m_z + b.m_z);
}

inline Vector3 Vector3::operator-(const Vector3 &b)const
{
	return Vector3(m_x - b.m_x, m_y - b.m_y, m_z - b.m_z);
}

inline Vector3 operator*(const Vector3&a,const double &b)
{
	return Vector3(a.m_x*b, a.m_y*b, a.m_z*b);
}

inline Vector3 operator*(const double &b, const Vector3&a)
{
	return Vector3(a.m_x*b, a.m_y*b, a.m_z*b);
}

inline Vector3 Vector3::operator/(const double &b)const
{
	return Vector3(m_x/b, m_y/b, m_z/b);
}

bool Vector3 ::operator==(const Vector3 &b)const
{
	if (m_x == b.m_x
		&&m_y == b.m_y
		&&m_z == b.m_z)
	{
		return true;
	}
	return false;
}

Vector3 Vector3::normalized()const
{
	return (*this)*(magnitude() > 0 ? 1 / magnitude() : 0.0);
}

inline double Vector3::Dot(const Vector3 &a, const Vector3& b)
{
	return (a.m_x*b.m_x + a.m_y*b.m_y + a.m_z*b.m_z);
}

inline Vector3 Vector3::Cross(const Vector3&a, const Vector3 &b)
{
	double x = a.m_y*b.m_z - a.m_z*b.m_y;
	double y = a.m_z*b.m_x - a.m_x*b.m_z;
	double z = a.m_x*b.m_y - a.m_y*b.m_x;
	return Vector3(x, y, z);
}

inline double Vector3::Angle(const Vector3&a, const Vector3&b)
{
	return acos(std::max(std::min(Vector3::Dot(a, b) / (a.magnitude()*b.magnitude()), 1.0), -1.0));
}

Quaternion::Quaternion()
{
	m_w = -1.0;
	m_x = 0.0;
	m_y = 0.0;
	m_z = 0.0;
}

Quaternion::Quaternion(const double &w, const double &x, const double &y, const double &z)
{
	m_w = w;
	m_x = x;
	m_y = y;
	m_z = z;
}

Quaternion::Quaternion(const Vector3 &pivot, const double &rad)
{
	double theta = rad / 2;
	Vector3 u;
	if (pivot.magnitude()> 0.0)
	{
		u = pivot.normalized();
	}
	else
	{
		u = Vector3(1, 0, 0);
		theta = 0.0;
	}

	m_w = cos(theta);
	m_x = sin(theta)*u.x();
	m_y = sin(theta)*u.y();
	m_z = sin(theta)*u.z();
}



inline Quaternion Quaternion::operator*(const Quaternion &b)const
{
	double w1 = m_w;
	double w2 = b.m_w;
	Vector3 v1 = Vector3(m_x, m_y, m_z);
	Vector3 v2 = Vector3(b.m_x, b.m_y, b.m_z);
	double w3 = w1*w2 - Vector3::Dot(v1, v2);
	Vector3 v3 = Vector3::Cross(v1, v2) + w1*v2 + w2*v1;
	return Quaternion(w3, v3.x(), v3.y(), v3.z());
}

inline Quaternion Quaternion::inverse()const
{
	return Quaternion(m_w, -m_x, -m_y, -m_z);
}

inline Vector3 Quaternion::rotate(const Vector3 &v)const
{
	Vector3 u = Vector3(m_x, m_y, m_z);
	double s = m_w;
	return 2.0*Vector3::Dot(u, v)*u + (s*s - Vector3::Dot(u, u))*v + 2.0*s*Vector3::Cross(u, v);
}