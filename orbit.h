#pragma once
#include<string>
#include<vector>
#include<map>
#ifndef PY_BUILD
#include"kepler.h"
#else
#include"kepler.cpp"
#endif
struct state
{
	Vector3 r;
	Vector3 v;
	std::string body;
	double t;
};
//��������ϵ
using namespace std;

class orbit
{
	friend class bodies;
private:
	//������������
	string m_Body;
	//��������
	double m_Gm;
	//�볤�ᣨ�ף�
	double m_Sem;
	//ƫ����
	double m_Ecc;
	//��ǣ����ȣ�
	double m_Inc;
	//�����㣨���ȣ�
	double m_Lan;
	//���ص���ǣ����ȣ�
	double m_Aop;
	//�����ʼʱ��
	double m_T0;
	//ʱ��Ϊ0��ʱƽ����
	double m_M0;
	//������һ��������µ�Ӱ�����ڣ���ʱ�̣��룩,-1��ʾ������
	double m_t_next;
	//���ڣ��룩
	double m_period;
	//�Ƕ���
	Vector3 m_H;
	//���ص�����
	Vector3 m_Pe;
	//Զ�ص�����
	Vector3 m_Ap;
	//��һ��Ӱ����Ĺ��
	orbit * m_Next;
public:
	orbit();

	orbit(const orbit & ob);

	orbit(const orbit &&ob);

	orbit(Vector3 r, Vector3 v, double t, double gm);

	//orbit(double sem,double ecc,double inc,double lan,double aop,double m0,double t0,double gm);

	~orbit();

	void reset_orbit(Vector3 r, Vector3 v, double t, double gm);

	void reset_orbit(Vector3 r, Vector3 v, double t, string body,const int&round=3);

	void reset_orbit(double sem,double ecc,double inc,double lan,double aop,double m0,double t0,double gm);

	void reset_orbit(double sem,double ecc,double inc,double lan,double aop,double m0,double t0,std::string body,const int &round =2);

	orbit& operator =(const orbit &ob);

	//orbit& operator =(const orbit &&ob);

public:
	Vector3 position_at_t(double t)const;

	Vector3 position_at_f(double f)const;

	double f_at_position(Vector3 pos)const;

	double f_at_t(double t)const;

	double f_at_r(double r)const;

	double t_to_f(double t0,double f)const;

	state Kstate_at_t(double t)const;

	state state_at_t(double t)const;

	state state_at_f(double F)const;

	const orbit * next_orbit()const;

	//double friend find_closest_t(const orbit &a, const orbit &b,int round=1);

public:
	
	void print()const;

	std::string body()const; 

	double period()const;

	Vector3 apoapsis()const;

	Vector3 periapsis()const; 

	Vector3 h()const;

	Vector3 angular_momentum()const;

	double eccentric()const;

	double semimajor_axis()const;

	double longitude_of_ascend_node()const;

	double incliantion()const;

	double argument_of_perigee()const;

	double gravity_parameter()const;

	double mean_anomaly0()const;

	double t0()const;

	double conic_a()const;
	
	double conic_b()const;

	//bƽ�����
	bool b_parameter(double &bt,double &br)const;
private:
	void set_body_name(std::string name);

	void count_next_orbit(const int &round);
};

struct celestial_body
{
	double gm;
	double radius;
	double atmosphere_depth;
	double soi;
	double init_rotation;
	double rotate_speed;
	string name;
	string parent;
	vector<string> satellites;
	Quaternion rotation;
	orbit orbit;
	celestial_body() 
	{ 
		gm = 0; 
		radius = 0; 
		atmosphere_depth = 0;
		soi = 0;
		init_rotation=0; 
		rotate_speed = 0.0;
		name = ""; 
		parent = "";
	}

	Vector3 msl_position(const double &lon,
					 const double &lat,
					 const double &alt,
					 const double &t)
	{
		double r=lat+radius;
		double coslat=cos(lat);
		double x=cos(lon)*coslat*r;
		double z=sin(lon)*coslat*r;
		double y=sin(lat)*r;
		Vector3 ret= Vector3(x,y,z);
		ret=Quaternion(Vector3(0,1,0),-(rotate_speed*t+init_rotation)).rotate(ret);
		return ret;
	}
	/*
	celestial_body(const celestial_body&&b)
	{
		gm=b.gm;
		radius=b.radius;
		atmosphere_depth=b.atmosphere_depth;
		soi=b.soi;
		name=b.name;
		parent=b.parent;
		for(int i=0;i<satellites.size();i++)
		{
			satellites.push_back(b.satellites[i]);
		}
		rotation=b.rotation;
		orbit=b.orbit;
	}
	*/
};

class bodies
{
private:
	map<string, celestial_body>m_bodies;
private:
	bodies();
public:
	static bodies &instance();
	celestial_body operator[](string name);
};
