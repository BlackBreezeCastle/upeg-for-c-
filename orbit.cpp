#include "orbit.h"
#include "config.h"
#include <queue>



inline double orbit_dis(const orbit &a, const orbit &b, const double& t)
{
	return  (a.position_at_t(t) - b.position_at_t(t)).magnitude();
}

double closest_dichotomy(const orbit &a, const orbit &b, const double& t0,const double &tf,double *closest_dis=NULL)
{
	double lt = t0;
	double lr = tf;
	double tmp1;
	double tmp2;
	double disl;
	double disr;
	while (fabs(lt-lr)>1)
	{
		tmp1 = 0.5*(lt + lr);
		tmp2 = tmp1 + 0.1;
		disl = (a.position_at_t(tmp1) - b.position_at_t(tmp1)).magnitude();
		disr = (a.position_at_t(tmp2) - b.position_at_t(tmp2)).magnitude();
		if (disl > disr)
		{
			lt = tmp1;
		}
		else
		{
			lr = tmp1;
		}
	}
	if(closest_dis!=NULL)
	{
		*closest_dis=0.5*(disl+disr);
	}
	return 0.5*(lt + lr);
}

double orbit::period()const 
{
	return m_period;
}

double closest_distance(const orbit &a, const orbit &b, double start_time,double *closest_t=NULL)
{
	double ret=0.0;
	double min_period;
	double max_period;
	if(a.period() < b.period())
	{
		min_period=a.period();
		max_period=b.period();
	}
	else
	{
		min_period=b.period();
		max_period=a.period();		
	}
	
	double closest_dis=INFINITY_DOUBLE;
	while(true)
	{
		double dis;
		double time=closest_dichotomy(a,b,start_time,start_time+min_period,&dis);
		if(closest_dis<dis)
		{
			//printf("%f\n",closest_dis);
			break;
		}
		else
		{	
			if(closest_t!=NULL)
			{
				*closest_t=time;
			}
			closest_dis=dis;
		}
		start_time+=min_period;
	}
	return closest_dis;
}

double solve_soi_bsp(const orbit &planet,const orbit &vessel,const double &min_ut,const double &max_ut,const double &soi)
{
	if(min_ut>=max_ut)
	{
		return INFINITY_DOUBLE;
	}
	double ret=0.5*(min_ut+max_ut);
	double dt=0.25*(max_ut-min_ut);
	double dis1=orbit_dis(planet,vessel,min_ut);
	double dis2=orbit_dis(planet,vessel,max_ut);
	double sig=1;
	if(dis1>dis2)
	{
		sig=-1;
	}
	double dis=INFINITY_DOUBLE;
	int count =0;
	while (fabs(dis-soi)>0.5)
	{
		count++;
		dis=orbit_dis(planet,vessel,ret);
		if(dis<soi)
		{
			ret+=sig*dt;
		}
		else
		{
			ret-=sig*dt;
		}
		dt=dt*0.5;
	}
	return ret;
}

double t_to_capture(const orbit&obt,const double &ut,std::string &new_body)
{
	double ret=INFINITY_DOUBLE;
	new_body="";
	celestial_body body=bodies::instance()[obt.body()];
	if(body.gm==0)
	{
		return ret;
	}

	for(int i=0;i<body.satellites.size();i++)
	{
		celestial_body satelite=bodies::instance()[body.satellites[i]];
		double t_closest;
		double dis_closest=closest_distance(obt,satelite.orbit,ut,&t_closest);
		if(t_closest<ut)
		{
			continue;
		}
		if(dis_closest<satelite.soi)
		{
			double t_next_soi=solve_soi_bsp(satelite.orbit,obt,ut,t_closest,satelite.soi);
			if(t_next_soi<ret)
			{
				ret=t_next_soi;
				new_body=satelite.name;
			}
		}
	}
	return ret;
}

double t_to_escape(const orbit &obt,const double &ut,std::string &new_body)
{
	auto body=bodies::instance()[obt.body()];
	if((obt.eccentric()<1.0&&obt.apoapsis().magnitude()<body.soi)||obt.body()=="Sun")
	{
		return INFINITY_DOUBLE;
	}
	new_body=body.parent;
	double f_escape=obt.f_at_r(body.soi);
	double ret=obt.t_to_f(ut,f_escape);
	if(ret>0)
	{
		return ret+ut;
	}
	else
	{
		return INFINITY_DOUBLE;
	}
}

void orbit::count_next_orbit(const int &round)
{
	if(round<2)
	{
		m_t_next=INFINITY_DOUBLE;
		return;
	}

	orbit *obt=this;
	double ut=m_T0;
	std::string escape_body;
	std::string capture_body;
	double t_escape=t_to_escape((*obt),ut,escape_body);
	double t_capture=t_to_capture((*obt),ut,capture_body);
	if(t_escape==INFINITY_DOUBLE&&t_capture==INFINITY_DOUBLE)
	{
		m_t_next=INFINITY_DOUBLE;
		return ;
	}

	if(t_escape<t_capture)
	{
		m_t_next=t_escape;
		state planet=bodies::instance()[(*obt).body()].orbit.state_at_t(t_escape);
		state vessel=(*obt).state_at_t(t_escape);
		auto tmp=(*obt).position_at_t(t_escape);
		auto dis=tmp.magnitude();
		Vector3 r=vessel.r+planet.r;
		Vector3 v=vessel.v+planet.v;
		(*obt).m_Next=new orbit();
		(*obt).m_Next->reset_orbit(r,v,t_escape,escape_body,round-1);
		return;
	}
	else
	{
		m_t_next=t_capture;
		state planet=bodies::instance()[capture_body].orbit.state_at_t(t_capture);
		state vessel=(*obt).state_at_t(t_capture);
		Vector3 r=vessel.r-planet.r;
		Vector3 v=vessel.v-planet.v;
		(*obt).m_Next=new orbit();
		(*obt).m_Next->reset_orbit(r,v,t_capture,capture_body,round-1);
		return ;
	}
}

orbit::orbit()
{
	m_Body="";
	m_Gm=0.0;
	m_Sem=0.0;	
	m_Ecc=0.0;
	m_Inc=0.0;
	m_Lan=0.0;
	m_Aop=0.0;
	m_T0=0.0;
	m_M0=0.0;
	m_t_next =INFINITY_DOUBLE;
	Vector3 m_H=Vector3(0.0,0.0,0.0);
	Vector3 m_Pe=Vector3(0.0,0.0,0.0);
	m_Next=NULL;
}
orbit::orbit(Vector3 r, Vector3 v, double t, double gm)
{
	m_Next=NULL;
	reset_orbit(r, v, t, gm);
}



orbit::~orbit()
{
	if(m_Next!=NULL)
	{
		delete m_Next;
		m_Next=NULL;
	}
}

void orbit::reset_orbit(Vector3 r, Vector3 v, double t, std::string _body,const int &round)
{
	m_Body = _body;
	reset_orbit(r,v,t,bodies::instance()[m_Body].gm);
	count_next_orbit(round);
}

orbit&orbit::operator=(const orbit &ob)
{
	if(this!=&ob)
	{
		if(this->m_Next!=NULL)
		{
			delete m_Next;
		}
		this->m_Body.clear();
		memcpy(this,&ob,sizeof(ob));
		this->m_Body=ob.body();
		if(ob.m_Next!=NULL)
		{
			this->m_Next=new orbit();
			*(this->m_Next)=*(ob.m_Next);
		}
	}
	return *this;
}

void orbit::set_body_name(std::string name)
{
	m_Body=name;
}

void orbit::reset_orbit(double sem,double ecc,double lan,double inc,double aop,double m0,double t0,double gm)
{
	m_Sem=sem;
	m_Ecc=ecc;
	m_Lan=lan;
	m_Inc=inc;
	m_Aop=aop;
	m_Gm=gm;
	m_T0=t0;
	m_M0=m0;
	//count the period
	m_period = 2 * PI*pow((pow(m_Sem,3) /m_Gm),0.5);

	//
	Vector3 y = Vector3(0.0, 1.0, 0.0);
	Vector3 x = Vector3(1.0, 0.0, 0.0);
	//计算角动量
	m_H=Quaternion(x,-m_Inc).rotate(y);
	m_H=Quaternion(y,-m_Lan).rotate(m_H);
	double E=0.0;
	if(m_Ecc!=1.0)
	{
		E=-0.5*m_Gm/m_Sem;
	}

	double pe=m_Sem*(1-m_Ecc);
	m_H=-sqrt(2.0*(E+gm/pe))*pe*m_H;
	//计算远地点，近地点
	m_Pe=Quaternion(y,-m_Lan).rotate(x);
	m_Pe=Quaternion(m_H,m_Aop).rotate(m_Pe);
	m_Ap=-1.0*m_Sem*(1+m_Ecc)*m_Pe;
	m_Pe=m_Pe*pe;

	double n = pow(fabs(m_Gm / (m_Sem*m_Sem*m_Sem)), 0.5);
}

void orbit::reset_orbit(double sem,double ecc,double lan,double inc,double aop,double m0,double t0,std::string body,const int &round)
{
	m_Body=body;
	reset_orbit(sem,ecc,lan,inc,aop,m0,t0,bodies::instance()[body].gm);
	count_next_orbit(round);
}

void orbit::reset_orbit(Vector3 r, Vector3 v, double t, double gm)
{
	m_Next = NULL;
	m_t_next =INFINITY_DOUBLE;
	m_Gm = gm;
	m_T0=t;
	double F0=0.0;
	//energy
	double E = 0.5*pow(v.magnitude(), 2) - m_Gm / r.magnitude();

	Vector3 y = Vector3(0.0, 1.0, 0.0);
	Vector3 x = Vector3(1.0, 0.0, 0.0);
	m_H = Vector3::Cross(r, v);
	double h2 = Vector3::Dot(m_H, m_H);
	if (E != 0.0)
	{
		m_Sem = -0.5*m_Gm / E;
		m_Ecc = pow(max(1 -h2/ (m_Gm*m_Sem),0.0), 0.5);
	}

	m_period = 2 * PI*pow((pow(m_Sem,3) /m_Gm),0.5);

	if(m_Ecc!=0.0)
	{
		F0 = R_to_F(m_Ecc,m_Sem,r.magnitude());
	}
	else
	{
		F0=0.0;
	}

	if (Vector3::Dot(r, v) < 0)
	{
		F0 = -F0;
	}

	m_Pe = Quaternion(m_H, -F0).rotate(r).normalized();
	m_Ap=-1.0*m_Pe*(1+m_Ecc)*m_Sem;
	m_Pe=m_Pe*(1-m_Ecc)*m_Sem;

	m_Inc = PI - Vector3::Angle(m_H, y);

	Vector3 vlan = Vector3::Cross(y, m_H);
	if(vlan.magnitude()<1e-16)
	{
		vlan=m_Pe;
	}

	m_Lan = Vector3::Angle(vlan, x);
	if (Vector3::Dot(m_H, y) > 0.0)
	{
		m_Lan = -m_Lan;
	}


	m_Aop = Vector3::Angle(m_Pe, vlan);
	if (Vector3::Dot(Vector3::Cross(m_Pe, vlan), m_H) > 0)
	{
		m_Aop = 2*PI-m_Aop;
	}

	double M = E_to_M(m_Ecc, F_to_E(m_Ecc, F0));
	double n = pow(fabs(m_Gm / (m_Sem*m_Sem*m_Sem)), 0.5);
	m_M0=M-n*t;
	//printf("\nC++ E M F\n%.17lf %.17lf %.17lf\n", E, M, m_F0);
	//printf("%.17lf %lf %.17lf %.17lf %.17lf %.17lf\n",m_Ecc,m_Sem,m_Inc,m_Lan,m_Aop,m_F0);
}

Vector3 orbit::position_at_t(double t)const
{
	if(t>m_t_next)
	{
		return m_Next->position_at_t(t);
	}
	double F=f_at_t(t);
	return position_at_f(F);
}

Vector3 orbit::position_at_f(double f)const
{
	return Quaternion(m_H, f).rotate(m_Pe).normalized() * F_to_R(m_Ecc, m_Sem, f);
}

double orbit::f_at_t(double t)const
{
	double n = pow(fabs(m_Gm / (m_Sem*m_Sem*m_Sem)), 0.5);
	double M = n*t+m_M0;
	M = fmod(M, PI * 2);
	if (M > PI&&m_Ecc<1.0)
	{
		M = M - 2*PI;
	}
	double E = M_to_E(m_Ecc, M);
	double F = E_to_F(m_Ecc,E);
	//printf("\nC++ E M F\n%.17lf %.17lf %.17lf\n", E, M, F);
	return F;
}

double orbit::f_at_r(double r)const
{
	return R_to_F(m_Ecc,m_Sem,r);
}

double orbit::t_to_f(double t0,double f)const
{
	double E=F_to_E(m_Ecc,f);
	double M=E_to_M(m_Ecc,E);
	double n = pow(fabs(m_Gm / (m_Sem*m_Sem*m_Sem)), 0.5);
	double t = (M-m_M0)/n-t0;
	return t;
}

state orbit::state_at_t(double t)const
{
	if(t>m_t_next)
	{
		return m_Next->state_at_t(t);
	}
	state ret;
	ret.t = t;
	double F = f_at_t(t);
	ret.r= position_at_f(F);

	Vector3 wr = Vector3::Cross(m_H, ret.r).normalized()*m_H.magnitude()/ret.r.magnitude();
	double tmp = 1 + m_Ecc*cos(F);
	tmp = tmp*tmp;
	double dr_dt = m_Sem*(m_Ecc*m_Ecc-1)*m_Ecc*sin(-F) / tmp;
	ret.v = dr_dt*wr.magnitude()/ret.r.magnitude() *ret.r.normalized() + wr;
	return ret;
}

void orbit::print()const
{
	printf("\n m_Ecc:%.17lf \n m_Sem:%.17lf \n m_Inc:%.17lf \n m_Lan:%.17lf \n m_Aop:%.17lf \n m_M0 %.17lf\n",m_Ecc,m_Sem,m_Inc*180/PI,m_Lan*180/PI,m_Aop*180/PI,m_M0);

}

std::string orbit::body()const
{
	return m_Body;
}

const orbit* orbit::next_orbit()const
{
	return m_Next;
}

Vector3 orbit::apoapsis()const
{
	return m_Ap;
}

Vector3 orbit::periapsis()const
{
	return m_Pe;
}

double orbit::eccentric()const
{
	return m_Ecc;
}

double orbit::semimajor_axis()const
{
	return m_Sem;
}
double orbit::longitude_of_ascend_node()const
{
	return m_Lan;
}

double orbit::incliantion()const
{
	return m_Inc;
}

double orbit::argument_of_perigee()const
{
	return m_Aop;
}

double orbit::gravity_parameter()const
{
	return m_Gm;
}

double orbit::mean_anomaly0()const
{
	return m_M0;
}

bodies::bodies()
{
	Config cfg("solar_config.txt");
	queue<string>q;
	q.push("Sun");
	while (!q.empty())
	{
		string name= q.front();
		auto tmp = cfg[name];
		celestial_body b;
		b.name = name;
		b.gm = tmp["gm"].asDouble();
		b.soi = tmp["soi"].asDouble();
		b.radius = tmp["radius"].asDouble();
		b.atmosphere_depth = tmp["atmosphere_depth"].asDouble();
		b.parent = tmp["parent"].asString();
		b.rotation = Quaternion(tmp["quaternion"][0].asDouble(), tmp["quaternion"][1].asDouble(), tmp["quaternion"][2].asDouble(), tmp["quaternion"][3].asDouble());

		if (name != "Sun")
		{
			double parent_gm = m_bodies.find(b.parent)->second.gm;
			double t0 = tmp["t0"].asDouble();
			Vector3 position(tmp["position"][0].asDouble(), tmp["position"][1].asDouble(), tmp["position"][2].asDouble());
			Vector3 velocity(tmp["velocity"][0].asDouble(), tmp["velocity"][1].asDouble(), tmp["velocity"][2].asDouble());
			b.orbit.reset_orbit(position, velocity, t0, parent_gm);
			b.orbit.set_body_name(b.parent);
		}

		for (int i = 0; i < tmp["satellites"].size();i++)
		{
			q.push(tmp["satellites"][i].asString());
			b.satellites.push_back(tmp["satellites"][i].asString());
		}

		m_bodies.insert(map<string, celestial_body>::value_type(b.name, b));
		q.pop();
	}
}

bodies &bodies::instance()
{
	static bodies m_pbodies;
	return m_pbodies;
}

celestial_body bodies::operator[](string name)
{
	auto it = m_bodies.find(name);
	if (it != m_bodies.end())
	{
		return it->second;
	}
	celestial_body invalid;
	return invalid;
}





