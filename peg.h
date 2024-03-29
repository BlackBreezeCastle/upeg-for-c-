#pragma once
#include "kmath.h"
#include<algorithm>
#include<vector>
#include<float.h>
#include"orbit.h"

#define FIXED_ACCELERATION_MODEL 1
#define MAX_THRUST_MODEL 0
const double g0 = 9.8067;
const double pi = 3.1415926535897932;



class pegas
{
	struct target
	{
		double angle;
		Vector3 normal;
		double radius;
		double velocity;
		target() { memset(this, 0, sizeof(*this)); }
	};
	struct state
	{
		double time;
		double mass;
		Vector3 radius;
		Vector3 velocity;
		state() { memset(this, 0, sizeof(*this)); }
	};

	struct  previous
	{
		Vector3 rbias;
		Vector3 rd;
		Vector3 rgrav;
		double	time;
		Vector3 v;
		Vector3 vgo;
		double	tgo;
		previous() { memset(this, 0, sizeof(*this)); }
	};

	struct stage
	{
		double massWet;
		double massDry;
		double gLim;
		double isp;
		double thrust;
		int mode;
		stage() { memset(this, 0, sizeof(*this)); }
	};
private:
	double m_u;// m_vessel.orbit.body.gravitational_parameter;
	std::vector<stage>m_stages;
	target m_target;
	state m_state;
	previous m_previous;
	double m_gLim;
	double m_output[2];// = (m_vessel.flight().pitch, m_vessel.flight().heading);
	double m_dvToGo;
	double m_final_weight;
public:
	pegas()
	{
		m_gLim = 0;
		m_dvToGo = 0;
		m_final_weight = 0;
	}

	void _upfg(int n)
	{
		auto stages = m_stages;
		auto gamma = m_target.angle;
		auto iy = m_target.normal;
		auto rdval = m_target.radius;
		auto vdval = m_target.velocity;
		auto t = m_state.time;
		auto m = m_state.mass;
		auto r = m_state.radius;
		auto v = m_state.velocity;

		auto rbias = m_previous.rbias;
		auto rd = m_previous.rd;
		auto rgrav = m_previous.rgrav;
		auto tp = m_previous.time;
		auto vprev = m_previous.v;
		auto vgo = m_previous.vgo;

		//1
		std::vector<int>SM;
		std::vector<double>aL;
		std::vector<double>ve;
		std::vector<double>fT;
		std::vector<double>aT;
		std::vector<double>tu;
		std::vector<double>tb;

		for (int i = 0; i < n; i++)
		{
			SM.push_back(stages[i].mode);
			aL.push_back(stages[i].gLim * g0);
			fT.push_back(stages[i].thrust);
			ve.push_back(stages[i].isp * g0);
			aT.push_back(fT[i] / stages[i].massWet);
			tu.push_back(ve[i] / aT[i]);
			if (stages[i].mode == MAX_THRUST_MODEL)
			{
				tb.push_back((stages[i].massWet - stages[i].massDry) * ve[i] / fT[i]);
			}
			else
			{
				tb.push_back(ve[i] * log(stages[i].massWet / stages[i].massDry) / aL[i]);
			}
		}
		//2
		double dt = t - tp;
		auto dvsensed = v - vprev;
		vgo = vgo - dvsensed;

		//3
		if (SM[0] == MAX_THRUST_MODEL)
		{
			aT[0] = fT[0] / m;
		}
		else
		{
			aT[0] = aL[0];
		}

		tu[0] = ve[0] / aT[0];
		double L = 0.0;
		std::vector<double>Li;
		for (int i = 0; i < n - 1; i++)
		{
			if (SM[i] == MAX_THRUST_MODEL)
			{
				Li.push_back(ve[i] * log(tu[i] / (tu[i] - tb[i])));
			}
			else
			{
				Li.push_back(aL[i] * tb[i]);
			}
			L = L + Li[i];
			if (L > vgo.magnitude())
			{
				return _upfg(n - 1);
			}
		}
		Li.push_back(vgo.magnitude() - L);
		int finalIndex = Li.size() - 1;
		if (finalIndex > -1)
		{
			m_final_weight =stages[finalIndex].massWet *exp(- Li[finalIndex]/ve[finalIndex]);
		}
		std::vector<double>tgoi;
		for (int i = 0; i < n; i++)
		{
			if (SM[i] == MAX_THRUST_MODEL)
			{
				tb[i] = tu[i] * (1 - exp((-Li[i] / ve[i])));
			}
			else
			{
				tb[i] = Li[i] / aL[i];
			}
			if (i == 0)
			{
				tgoi.push_back(tb[i]);
			}
			else
			{
				tgoi.push_back(tgoi[i - 1] + tb[i]);
			}
		}
		double tgo = tgoi[n - 1];


		L = 0.0;
		double J = 0.0;
		double S = 0.0;
		double Q = 0.0;
		double H = 0.0;
		double P = 0.0;

		std::vector<double>Ji;
		std::vector<double>Si;
		std::vector<double>Qi;
		std::vector<double>Pi;
		double tgoi1 = 0.0;

		for (int i = 0; i < n; i++)
		{
			if (i > 0)
			{
				tgoi1 = tgoi[i - 1];
			}

			if (SM[i] == MAX_THRUST_MODEL)
			{
				Ji.push_back(tu[i] * Li[i] - ve[i] * tb[i]);
				Si.push_back(-Ji[i] + tb[i] * Li[i]);
				Qi.push_back(Si[i] * (tu[i] + tgoi1) - 0.5 * ve[i] * pow(tb[i], 2));
				Pi.push_back(Qi[i] * (tu[i] + tgoi1) - 0.5 * ve[i] * pow(tb[i], 2) * (tb[i] / 3 + tgoi1));
			}
			else
			{
				Ji.push_back(0.5 * Li[i] * tb[i]);
				Si.push_back(Ji[i]);
				Qi.push_back(Si[i] * (tb[i] / 3 + tgoi1));
				Pi.push_back((1 / 6) * Si[i] * (pow(tgoi[i], 2) + 2 * tgoi[i] * tgoi1 + 3 * pow(tgoi1, 2)));
			}

			Ji[i] = Ji[i] + Li[i] * tgoi1;
			Si[i] = Si[i] + L * tb[i];
			Qi[i] = Qi[i] + J * tb[i];
			Pi[i] = Pi[i] + H * tb[i];

			L = L + Li[i];
			J = J + Ji[i];
			S = S + Si[i];
			Q = Q + Qi[i];
			P = P + Pi[i];
			H = J * tgoi[i] - Q;
		}
		m_dvToGo = L;

		//5
		auto _lambda = vgo.normalized();
		if (m_previous.tgo > 0)
		{
			rgrav = pow((tgo / m_previous.tgo), 2) * rgrav;
		}

		Vector3 rgo = rd - (r + v * tgo + rgrav);
		auto iz = Vector3::Cross(rd, iy).normalized();
		auto rgoxy = rgo - Vector3::Dot(iz, rgo) * iz;
		auto rgoz = (S - Vector3::Dot(_lambda, rgoxy)) / Vector3::Dot(_lambda, iz);
		rgo = rgoxy + rgoz * iz + rbias;
		auto lambdade = Q - S * J / L;
		auto lambdadot = (rgo - S * _lambda) / lambdade;
		Vector3 iF_ = _lambda - lambdadot * J / L;
		iF_ = iF_.normalized();
		double phi = Vector3::Angle(iF_, _lambda);
		auto phidot = -phi * L / J;
		auto vthrust = (L - 0.5 * L * phi * phi - J * phi * phidot - 0.5 * H * phidot * phidot) * _lambda;
		auto rthrust = (S - 0.5 * S * phi * phi - Q * phi * phidot - 0.5 * P * phidot * phidot) * _lambda;
		rthrust = rthrust - (S * phi + Q * phidot) * lambdadot.normalized();
		auto vbias = vgo - vthrust;
		rbias = rgo - rthrust;
		rbias = rbias;
		vbias = vbias;

		//	6;
		//	TODO angle rates;
		auto _up = r.normalized();
		auto _east = Vector3::Cross(_up, Vector3(0, 1, 0)).normalized();
		double pitch = pi / 2 - Vector3::Angle(iF_, _up);
		Vector3 inplane = iF_ - Vector3::Dot(_up, iF_) * _up;
		double yaw = Vector3::Angle(inplane, _east);
		auto tangent = Vector3::Cross(_up, _east);
		if (Vector3::Dot(inplane, tangent) < 0)
		{
			yaw = -yaw;
		}
		yaw = yaw - pi / 2;
		yaw = normalized_rad(yaw) + pi;
		m_output[0] = pitch;
		m_output[1] = yaw;
		//7

		auto rc1 = r - 0.1 * rthrust - (tgo / 30) * vthrust;
		auto vc1 = v + 1.2 * rthrust / tgo - 0.1 * vthrust;
		orbit ob = orbit(rc1, vc1, 0, m_u);
		auto pack = ob.state_at_t(tgo);
		rgrav = pack.r - rc1 - vc1 * tgo;
		Vector3 vgrav = pack.v - vc1;


		/*
		auto w2 = m_u / (pow(rdval, 3));
		auto vd = v + vgo;
		rgrav = -w2*pow(tgo, 2) * ((3 * rd + 7 * r)*0.05 - (2 * vd - 3 * v) / 30);
		auto vgrav = -w2*tgo*((rd + r)*0.5 - (vd - v)*tgo / 12);
		*/


		//	8
		Vector3 rp = r + v * tgo + rgrav + rthrust;
		rp = rp - Vector3::Dot(rp, iy) * iy;
		rd = rdval * rp.normalized();
		auto ix = rd.normalized();
		iz = Vector3::Cross(ix, iy);

		auto vd = (iz * cos(gamma) + ix * sin(gamma)) * vdval;
		vgo = vd - v - vgrav + vbias;

		//print(vbias);
		//print(rbias);
		//print('\n');

		m_previous.rbias = rbias;
		m_previous.rd = rd;
		m_previous.rgrav = rgrav;
		m_previous.time = m_state.time;
		m_previous.v = m_state.velocity;
		m_previous.vgo = vgo;
		m_previous.tgo = tgo;
	}
	void set_target_orbit(double inc, double lan, double radius, double velocity, double angle = 0.0)
	{
		m_target.angle = angle;
		m_target.normal = target_normal_vector(inc, lan);
		m_target.radius = radius;
		m_target.velocity = velocity;
	}

	void init_peg(Vector3 r, Vector3 v, double mass, double t, double gm)
	{
		m_u = gm;
		auto q = Quaternion(m_target.normal, -2 * PI / 30);
		Vector3 rd = q.rotate(r);
		Vector3 vd = Vector3::Cross(rd, m_target.normal).normalized() * m_target.velocity;
		m_previous.rd = rd;
		m_previous.time = t;
		m_previous.v = v;
		m_previous.vgo = vd - v;

		m_state.time = t;
		m_state.mass = mass;
		m_state.radius = r;
		m_state.velocity = v;
	}

	void _add_stage(double massWet, double massDry, double thrust, double isp, double gLim, double mode)
	{
		auto _stage = stage();
		_stage.massWet = massWet;
		_stage.massDry = massDry;
		_stage.gLim = gLim;
		_stage.isp = isp;
		_stage.thrust = thrust;
		_stage.mode = mode;
		m_stages.push_back(_stage);
	}
	//按一级，二级，三级含上级重量添加
	void add_stage(double massWet, double massDry, double thrust, double isp, double gLim = 4.5)
	{
		auto _stage = stage();

		if (thrust == 0 || isp == 0 || massWet == massDry)
		{
			return;
		}

		double mass_tmp = thrust / (gLim * g0);
		if (mass_tmp <= massDry)
		{
			_add_stage(massWet, massDry, thrust, isp, gLim, MAX_THRUST_MODEL);
		}
		else if (mass_tmp >= massWet)
		{
			_add_stage(massWet, massDry, thrust, isp, gLim, FIXED_ACCELERATION_MODEL);
		}
		else
		{
			_add_stage(massWet, mass_tmp, thrust, isp, gLim, MAX_THRUST_MODEL);
			_add_stage(mass_tmp, massDry, thrust, isp, gLim, FIXED_ACCELERATION_MODEL);
		}
	}
	void  update(double ut, double mass, Vector3 r, Vector3 v)
	{
		m_state.time = ut;
		m_state.mass = mass;
		m_state.radius = r;
		m_state.velocity = v;


		int n = m_stages.size();
		m_stages.front().massWet = m_state.mass;
		if (m_state.mass <= m_stages.front().massDry)
		{
			m_stages.pop_back();
			//vessel.control.throttle = 1.0
			//#m_throttle_bias.clear()
			return;
		}


		if (m_stages.front().mode != 0)
		{
			//double acc = vessel.thrust / max(vessel.mass, 0.1)
			//double dacc = acc - g0*m_stages[0].gLim
			//double dacc = max(-1, min(1, dacc))
			//double vessel.control.throttle = 1.0 - m_throttle_bias.integral(0.05*dacc)
		}

		_upfg(n);
		return;
	}


	double time_to_go()
	{
		return m_previous.tgo;
	}

	double __time_to_stage(stage _stage)
	{
		if (_stage.mode == 0)
		{
			double dm = _stage.thrust / (_stage.isp * g0);
			return (_stage.massWet - _stage.massDry) / dm;
		}
		else
		{
			double dv = _stage.isp * log(_stage.massWet / _stage.massDry);
			return dv / _stage.gLim;
		}
	}

	double stages_num()
	{
		return m_stages.size();
	}

	double angle_to_rd()
	{
		return Vector3::Angle(m_previous.rd, m_state.radius);
	}

	//Vector3 rd_position() 
	//{

	//}
	void set_max_g(double g)
	{
		m_gLim = g;
	}

	double dv_to_go()
	{
		return m_dvToGo;
	}

	double final_weight()
	{
		return m_final_weight;
	}

	void vehicle_info()
	{
		//for() 
		//	printf('wet mass%f dry mass%f thurst%f isp%f' ,i.massWet, i.massDry, i.thrust, i.isp);
		//	}
	}

public:

	static Vector3 target_normal_vector(double inc, double lan)
	{
		Vector3 res;
		Vector3 north_pole = Vector3(0, 1, 0);
		Vector3 vernal_vector = Vector3(1, 0, 0);
		Quaternion q = Quaternion(north_pole, -lan);
		Vector3 ascend_node = q.rotate(vernal_vector);
		Vector3 tmp = Vector3::Cross(north_pole, ascend_node).normalized();
		if (abs(inc - PI / 2) < 1e-16)
		{
			res = tmp;
		}
		else
		{
			res = north_pole + (tmp * tan(inc));
		}
		return res.normalized();
	}
};