#pragma once
#include<string>
#include<vector>
#include<map>
#include"kepler.h"
//#include"kepler.cpp"
//#include"config.m_H"
struct state
{
	Vector3 r;
	Vector3 v;
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

	orbit(Vector3 r, Vector3 v, double t, double gm);

	orbit(double sem,double ecc,double lan,double inc,double aop,double m0,double t0,double gm);

	~orbit();

	void reset_orbit(Vector3 r, Vector3 v, double t, double gm);

	void reset_orbit(Vector3 r, Vector3 v, double t, string body,const int&round=3);

	void reset_orbit(double sem,double ecc,double lan,double inc,double aop,double m0,double t0,double gm);

	void reset_orbit(double sem,double ecc,double lan,double inc,double aop,double m0,double t0,std::string body,const int &round =3);

	orbit& operator =(const orbit &ob);

public:
	Vector3 position_at_t(double t)const;

	Vector3 position_at_f(double f)const;

	double f_at_t(double t)const;

	double f_at_r(double r)const;

	double t_to_f(double t0,double f)const;

	state state_at_t(double t)const;

	const orbit * next_orbit()const;

	//double friend find_closest_t(const orbit &a, const orbit &b,int round=1);

public:
	
	void print()const;

	std::string body()const; 

	double period()const;

	Vector3 apoapsis()const;

	Vector3 periapsis()const; 

	double eccentric()const;

	double semimajor_axis()const;

	double longitude_of_ascend_node()const;

	double incliantion()const;

	double argument_of_perigee()const;

	double gravity_parameter()const;

	double mean_anomaly0()const;
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
	string name;
	string parent;
	vector<string> satellites;
	Quaternion rotation;
	orbit orbit;
	celestial_body(){gm=0;radius=0;atmosphere_depth=0,soi=0;name="";parent="";}
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
