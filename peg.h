#pragma once
#include "kmath.h"
#include<vector>
const double g0 = 9.8067;
const double pi = 3.141592654;
struct target
{
	double angle = 0.0;
	Vector3 normal = Vector3(0.0, 0.0, 0.0);
	double radius = 0.0;
	double velocity = 0.0;
};

struct state
{
	double time = 0.0;
	double mass = 0.0;
	Vector3 radius = Vector3(0.0, 0.0, 0.0);
	Vector3 velocity = Vector3(0.0, 0.0, 0.0);
};

struct  previous
{
	Vector3 rbias = Vector3(0.0, 0.0, 0.0);
	Vector3 rd = Vector3(0.0, 0.0, 0.0);
	Vector3 rgrav = Vector3(0.0, 0.0, 0.0);
	double	time = 0.0;
	Vector3 v = Vector3(0.0, 0.0, 0.0);
	Vector3 vgo = Vector3(0.0, 0.0, 0.0);
	double	tgo = 0.0;
};

struct stage 
{
	double massWet = 0.0;
	double massDry = 0.0;
	double gLim = 10.0;
	double isp = 0.0;
	double thrust = 0.0;
	double mode = 0;
};
class pegas
{
private:
	//double m_conn = conn;
	//m_vessel = conn.space_center.active_vessel;
	//m_earth = conn.space_center.bodies['Earth'];
	double m_u;// m_vessel.orbit.body.gravitational_parameter;
	//m_reference_frame = m_vessel.orbit.body.non_rotating_reference_frame;
	std::vector<stage>m_stages;;
	target m_target = target();
	state m_state = state();
	previous m_previous = previous();
	double m_tgo = 0.0;
	double m_last_stage_mass = 0.0;
	double m_gLim = 4.5;
	double m_output[2];// = (m_vessel.flight().pitch, m_vessel.flight().heading);

private:
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
		std::vector<double>SM;
		std::vector<double>aL;
		std::vector<double>ve;
		std::vector<double>fT;
		std::vector<double>aT;
		std::vector<double>tu;
		std::vector<double>tb;

		for (int i = 0; i < n; i++)
		{
			SM.push_back(stages[i].mode);
			aL.push_back(stages[i].gLim*g0);
			fT.push_back(stages[i].thrust);
			ve.push_back(stages[i].isp*g0);
			aT.push_back(fT[i] / stages[i].massWet);
			tu.push_back(ve[i] / aT[i]);
			if (stages[i].mode == 0)
			{
				tb.push_back((stages[i].massWet - stages[i].massDry)*ve[i] / fT[i]);
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
		if (SM[0] == 0)
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
			if (SM[i] == 0)
			{
				Li.push_back(ve[i] * log(tu[i] / (tu[i] - tb[i])));
			}
			else
			{
				Li.push_back(aL[i] * tb[i]);
				L = L + Li[i];

				if (L > vgo.magnitude())
				{
					return _upfg(n - 1);
				}
			}
		}
		Li.push_back(vgo.magnitude() - L);
		std::vector<double>tgoi;
		for (int i = 0; i < n; i++)
		{
			if (SM[i] == 0)
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

		double L = 0.0;
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

			if (SM[i] == 0)
			{
				Ji.push_back(tu[i] * Li[i] - ve[i] * tb[i]);
				Si.push_back(-Ji[i] + tb[i] * Li[i]);
				Qi.push_back(Si[i] * (tu[i] + tgoi1) - 0.5*ve[i] * pow(tb[i], 2));
				Pi.push_back(Qi[i] * (tu[i] + tgoi1) - 0.5*ve[i] * pow(tb[i], 2)* (tb[i] / 3 + tgoi1));
			}
			else
			{
				Ji.push_back(0.5*Li[i] * tb[i]);
				Si.push_back(Ji[i]);
				Qi.push_back(Si[i] * (tb[i] / 3 + tgoi1));
				Pi.push_back((1 / 6)*Si[i] * (pow(tgoi[i], 2) + 2 * tgoi[i] * tgoi1 + 3 * pow(tgoi1, 2)));
			}

			Ji[i] = Ji[i] + Li[i] * tgoi1;
			Si[i] = Si[i] + L*tb[i];
			Qi[i] = Qi[i] + J*tb[i];
			Pi[i] = Pi[i] + H*tb[i];

			L = L + Li[i];
			J = J + Ji[i];
			S = S + Si[i];
			Q = Q + Qi[i];
			P = P + Pi[i];
			H = J*tgoi[i] - Q;
		}

		//5
		auto _lambda = vgo.normalized();
		if (m_previous.tgo > 0)
		{
			rgrav = pow((tgo / m_previous.tgo), 2) * rgrav;
		}

		auto rgo = rd - (r + v*tgo + rgrav);
		auto iz = Vector3::Cross(rd, iy).normalized();
		auto rgoxy = rgo - Vector3::Dot(iz, rgo)*iz;
		auto rgoz = (S - Vector3::Dot(_lambda, rgoxy)) / Vector3::Dot(_lambda, iz);
		rgo = rgoxy + rgoz*iz + rbias;
		auto lambdade = Q - S*J / L;
		auto lambdadot = (rgo - S*_lambda) / lambdade;
		auto iF_ = _lambda - lambdadot*J / L;
		iF_ = iF_.normalized();
		auto phi = Vector3::Angle(iF_, _lambda);
		auto phidot = -phi*L / J;
		auto vthrust = (L - 0.5*L*phi**2 - J*phi*phidot - 0.5*H*phidot**2)*_lambda;
		auto rthrust = (S - 0.5*S*phi**2 - Q*phi*phidot - 0.5*P*phidot**2)*_lambda;
		rthrust = rthrust - (S*phi + Q*phidot)*lambdadot.normalized();
		auto vbias = vgo - vthrust;
		rbias = rgo - rthrust;
		rbias = rbias;
		vbias = vbias;

		//	6;
		//	TODO angle rates;
		auto _up = r.normalized();
		auto _east = Vector3::Cross(_up, Vector3(0, 1, 0)).normalized();
		auto pitch = pi / 2 - Vector3::Angle(iF_, _up);
		auto inplane = iF_ - Vector3::Dot(_up, iF_)*_up;
		auto yaw = Vector3::Angle(inplane, _east);
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
		/*
		auto rc1 = r - 0.1*rthrust - (tgo / 30)*vthrust;
		auto vc1 = v + 1.2*rthrust / tgo - 0.1*vthrust;
		auto pack = kepler.state_at(rc1.tuple3(), vc1.tuple3(), m_u, 0, tgo);
		rgrav = Vector3::Tuple3(pack[0]) - rc1 - vc1*tgo;
		vgrav = Vector3::Tuple3(pack[1]) - vc1;
		#print('p', (pack[0].magnitude(), pack[1].magnitude()))

		*/
		auto w2 = m_u / (pow(rdval, 3));
		auto vd = v + vgo;
		auto rgrav = -w2*pow(tgo, 2) * ((3 * rd + 7 * r)*0.05 - (2 * vd - 3 * v) / 30);
		auto vgrav = -w2*tgo*((rd + r)*0.5 - (vd - v)*tgo / 12);


		//	8
		auto rp = r + v*tgo + rgrav + rthrust;
		rp = rp - Vector3::Dot(rp, iy)*iy;
		rd = rdval*rp.normalized();
		auto ix = rd.normalized();
		iz = Vector3::Cross(ix, iy);
#	
		vd = (iz*cos(gamma) + ix*sin(gamma))*vdval;
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
		m_previous.tgo = tgo;;
	}
	void set_target(double inc, double lan, double radius, double velocity, double angle = 0.0)
	{
		m_target.angle = angle;
		m_target.normal; //= target_normal_vector(m_conn, m_earth, inc, lan, m_reference_frame)
		m_target.radius = radius;
		m_target.velocity = velocity;

		r = Vector3::Tuple3(m_vessel.position(m_reference_frame));
		v = Vector3::Tuple3(m_vessel.velocity(m_reference_frame));
		q = Quaternion.PivotRad(m_target.normal, radians(1));
		rd = q.rotate(r);
		vd = Vector3::Cross(rd, m_target.normal).normalized()*velocity;
		m_previous.rd = rd;
		m_previous.time = m_conn.space_center.ut
		m_previous.v = v
		m_previous.vgo = vd - v

		m_state.time = m_conn.space_center.ut
		m_state.mass = m_vessel.mass
		m_state.radius = Vector3::Tuple3(m_vessel.position(m_reference_frame))
		m_state.velocity = Vector3::Tuple3(m_vessel.velocity(m_reference_frame))

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
	void add_stage(double massWet, double massDry, double thrust, double isp, double gLim = 4.5)
	{
		auto _stage = stage();
		m_stages.reverse();
		last_stage_mass = m_last_stage_mass

			if thrust == 0 or isp == 0 or massWet == massDry
				m_last_stage_mass = last_stage_mass + massWet
				return None

				mass_tmp = thrust / (gLim*g0)
				if mass_tmp <= massDry + last_stage_mass
					m_add_stage(massWet + last_stage_mass, massDry + last_stage_mass, thrust, isp, gLim, 0)
					elif mass_tmp >= massWet + last_stage_mass
					m_add_stage(massWet + last_stage_mass, massDry + last_stage_mass, thrust, isp, gLim, 1)
				else
					m_add_stage(mass_tmp, massDry + last_stage_mass, thrust, isp, gLim, 1)
					m_add_stage(massWet + last_stage_mass, mass_tmp, thrust, isp, gLim, 0)
					m_last_stage_mass = last_stage_mass + massWet
					m_stages.reverse()
	}
	void  update(double ut, double mass, Vector3 r, Vector3 v)
	{
		m_state.time = m_conn.space_center.ut
			m_state.mass = m_vessel.mass
			m_state.radius = Vector3::Tuple3(m_vessel.position(m_reference_frame))
			m_state.velocity = Vector3::Tuple3(m_vessel.velocity(m_reference_frame))

			vessel = m_vessel
			n = len(m_stages)
			m_stages[0].massWet = m_state.mass
			if m_state.mass <= m_stages[0].massDry
				m_stages.pop(0)
				vessel.control.throttle = 1.0
				#m_throttle_bias.clear()
				return None

				'''
				if m_stages[0].mode != 0
					acc = vessel.thrust / max(vessel.mass, 0.1)
					dacc = acc - g0*m_stages[0].gLim
					dacc = max(-1, min(1, dacc))
					vessel.control.throttle = 1.0 - m_throttle_bias.integral(0.05*dacc)
					'''

					last_tgo = m_previous.tgo
					_upfg(n)
					return m_output
	}

	void update_stages(double thrustK = 1.0)
	{
		m_last_stage_mass = 0.0
		m_stages = []
		stages = get_stages(m_vessel.parts.root)
		for i in stages
			self.add_stage(i[0], i[1], i[2] * thrustK, i[3], m_gLim)

			def time_to_go(self)
			return m_previous.tgo

			def __time_to_stage(self, stage)
			if stage.mode == 0
				dm = stage.thrust / (stage.isp*g0)
				return (stage.massWet - stage.massDry) / dm
			else
				dv = stage.isp*log(stage.massWet / stage.massDry)
				return dv / stage.gLim

				def time_to_stage(self)
				stage = m_stages[0]
				ret = m_time_to_stage(stage)
				if len(m_stages) > 1
	}
if m_stages[1].massWet == stage.massDry 
	ret = ret + m_time_to_stage(m_stages[1])
	return ret

	def stages_num(self) 
	return len(m_stages)

	def angle_to_rd(self) 
	return Vector3::Angle(m_previous.rd, m_state.radius)

	def rd_position(self) 
	pos = m_previous.rd.tuple3()
	ref = m_reference_frame
	body = m_vessel.orbit.body
	turn_angle = degrees(body.rotational_speed*m_previous.tgo)
	return (body.longitude_at_position(pos, ref) - turn_angle, body.latitude_at_position(pos, ref))

	def set_max_g(self, g) 
	m_gLim = g

	def vehicle_info(self) 
	for i in m_stages 
		print('wet mass%f dry mass%f thurst%f isp%f' % (i.massWet, i.massDry, i.thrust, i.isp))
		}