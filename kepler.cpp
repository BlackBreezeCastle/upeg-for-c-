#include"kepler.h"
/* ---------------------------------------------------------------------
**
** 函数名称: KeplerStart3
** 函数功能: 开普勒方程三阶初值估计
** 输入参数: ecc -> 偏心率 [N/A]
**           M   -> 平近角 [rad]
** 输出参数: 无
** 返回值  : E   -> 偏近角 [rad]
**
** ---------------------------------------------------------------------*/
#define D2PI (3.14159265358979323846)
#define error 1.0e-14

double KeplerStart3(double ecc, double M) {

	double t34, t35, t33, E;

	if (ecc < 0 || ecc >= 1.0) {
		printf("!!!错误的偏心率，该程序只适合椭圆轨道(0 =< ecc < 1.0)!!!\n");
		return -1.0;
	}

	if (M >= D2PI || M < 0) {
		M = fmod(M, D2PI);
	}

	t34 = ecc * ecc;
	t35 = ecc * t34;
	t33 = cos(M);

	E = M + (-0.5*t35 + ecc + (t34 + 3.0 / 2.0*t33*t35) * t33) * sin(M);

	return E;
}

/* ---------------------------------------------------------------------
**
** 函数名称: eps3
** 函数功能: 开普勒方程三阶误差估计
** 输入参数: ecc -> 偏心率             [N/A]
**           M   -> 平近角             [rad]
**           x   -> 偏近角过程值       [rad]
** 输出参数: 无
** 返回值  : E_error -> 偏近角估计误差 [rad]
**
** ---------------------------------------------------------------------*/
double eps3(double ecc, double M, double x) 
{

	double t1, t2, t3, t4, t5, t6, E_error;

	if (ecc < 0.0 || ecc >= 1.0) {
		printf("!!!错误的偏心率，该程序只适合椭圆轨道(0 =< ecc < 1.0)!!!\n");
		return -1.0;
	}

	if (M >= D2PI || M < 0) {
		M = fmod(M, D2PI);
	}

	t1 = cos(x);
	t2 = -1.0 + ecc * t1;
	t3 = sin(x);
	t4 = ecc * t3;
	t5 = -x + t4 + M;
	t6 = t5 / (0.5 * t5 * t4 / t2 + t2);

	E_error = t5 / ((0.5*t3 - 1.0 / 6.0*t1*t6) * ecc * t6 + t2);

	//printf("误差%.40f",E_error);
	return E_error;
}

double E_to_F(double ecc,double E)
{
	if (ecc < 1.0)
	{
		return 2.0*atan(pow((1.0 + ecc) / (1.0 - ecc), 0.5)*tan(0.5*E));
	}
	else if (ecc > 1.0)
	{
		return 2.0*atan(pow((1.0 + ecc) / (ecc - 1.0), 0.5)*tanh(0.5*E));
	}
}

double F_to_E(double ecc, double F)
{
	if (ecc < 1.0)
	{
		return 2.0*atan(pow((1.0 - ecc) / (1.0 + ecc), 0.5)*tan(0.5*F));
	}
	else if (ecc > 1.0)
	{
		return 2.0* atanh(pow((ecc - 1.0) / (1.0 + ecc), 0.5)*tan(0.5*F));
	}
}


/* ---------------------------------------------------------------------
**
** 函数名称: KeplerSolve
** 函数功能: 求解开椭圆下普勒方程
** 输入参数: ecc -> 偏心率 [N/A]
**           M   -> 平近角 [rad]
** 输出参数: 无
** 返回值  : E   -> 偏近角 [rad]
**
** ---------------------------------------------------------------------*/
double kepler_E(double ecc, double M)
{

	double E0, dE, E, Mnorm;
	int    iter_count = 0;

	if (ecc < 0 || ecc >= 1.0) {
		printf("!!!错误的偏心率，该程序只适合椭圆轨道(0 =< ecc < 1.0)!!!\n");
		return -1.0;
	}

	if (M >= D2PI || M < 0) {
		M = fmod(M, D2PI);
	}

	Mnorm = fmod(M, D2PI);

	E0 = KeplerStart3(ecc, Mnorm);
	dE = error + 1.0;

	while (fabs(dE) > error)
	{
		E = E0 - eps3(ecc, Mnorm, E0);
		dE = fabs(E - E0);
		E0 = E;

		iter_count = iter_count + 1.0;

		if (iter_count > 10) {
			printf("太令人惊讶了，KeplerSolve 竟然不收敛！\n");
			return -1.0;
		}

		//printf("迭代次数：iter_count = %d\n", iter_count);
	}

	return E;
}

double kepler_H(double e, double M)
{
	double H = 2.0 * M / e;
    H = log(sqrt(H * H + 1.0) + H);
	double ratio = 1.0;
	while (fabs(ratio) >error)
	{
		ratio = (e*sinh(H) - H - M) / (e*cosh(H) - 1.0);
		H = H - ratio;
	}
	return H;
}

double M_to_E(double ecc, double M)
{
	if (ecc > 1.0)
	{
		return kepler_H(ecc, M);
	}
	else if (ecc < 1.0)
	{
		return kepler_E(ecc, M);
	}
}


double F_to_R(double ecc, double sem, double f)
{
	return sem*(1-ecc*ecc)/(1.0 + ecc*cos(f));
}

double R_to_F(double ecc, double sem, double r)
{
	return acos(std::min(std::max((sem*(1-ecc*ecc)/r - 1.0) / ecc, -1.0 ), 1.0));
}

double E_to_M(double ecc, double E)
{
	if (ecc > 1.0)
	{
		return sinh(E)*ecc - E;
	}
	else if (ecc < 1.0)
	{
		return E - sin(E)*ecc;
	}
}