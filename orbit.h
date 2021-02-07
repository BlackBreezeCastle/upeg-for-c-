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
//左手坐标系
using namespace std;
class orbit
{
	friend class bodies;
private:
	//中心天体名称
	string m_Body;
	//重力常数
	double m_Gm;
	//半长轴（米）
	double m_Sem;
	//偏心率
	double m_Ecc;
	//倾角（弧度）
	double m_Inc;
	//升交点（弧度）
	double m_Lan;
	//近地点辐角（弧度）
	double m_Aop;
	//轨道起始时间
	double m_T0;
	//时刻为0秒时平近角
	double m_M0;
	//到达下一条轨道（新的影响球内）的时刻（秒）,-1表示不存在
	double m_t_next;
	//周期（秒）
	double m_period;
	//角动量
	Vector3 m_H;
	//近地点向量
	Vector3 m_Pe;
	//远地点向量
	Vector3 m_Ap;
	//下一个影响球的轨道
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
