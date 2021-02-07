#include"kepler.h"
/* ---------------------------------------------------------------------
**
** ��������: KeplerStart3
** ��������: �����շ������׳�ֵ����
** �������: ecc -> ƫ���� [N/A]
**           M   -> ƽ���� [rad]
** �������: ��
** ����ֵ  : E   -> ƫ���� [rad]
**
** ---------------------------------------------------------------------*/
#define D2PI (3.14159265358979323846)
#define error 1.0e-14

double KeplerStart3(double ecc, double M) {

	double t34, t35, t33, E;

	if (ecc < 0 || ecc >= 1.0) {
		printf("!!!�����ƫ���ʣ��ó���ֻ�ʺ���Բ���(0 =< ecc < 1.0)!!!\n");
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
** ��������: eps3
** ��������: �����շ�������������
** �������: ecc -> ƫ����             [N/A]
**           M   -> ƽ����             [rad]
**           x   -> ƫ���ǹ���ֵ       [rad]
** �������: ��
** ����ֵ  : E_error -> ƫ���ǹ������ [rad]
**
** ---------------------------------------------------------------------*/
double eps3(double ecc, double M, double x) 
{

	double t1, t2, t3, t4, t5, t6, E_error;

	if (ecc < 0.0 || ecc >= 1.0) {
		printf("!!!�����ƫ���ʣ��ó���ֻ�ʺ���Բ���(0 =< ecc < 1.0)!!!\n");
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

	//printf("���%.40f",E_error);
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
** ��������: KeplerSolve
** ��������: ��⿪��Բ�����շ���
** �������: ecc -> ƫ���� [N/A]
**           M   -> ƽ���� [rad]
** �������: ��
** ����ֵ  : E   -> ƫ���� [rad]
**
** ---------------------------------------------------------------------*/
double kepler_E(double ecc, double M)
{

	double E0, dE, E, Mnorm;
	int    iter_count = 0;

	if (ecc < 0 || ecc >= 1.0) {
		printf("!!!�����ƫ���ʣ��ó���ֻ�ʺ���Բ���(0 =< ecc < 1.0)!!!\n");
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
			printf("̫���˾����ˣ�KeplerSolve ��Ȼ��������\n");
			return -1.0;
		}

		//printf("����������iter_count = %d\n", iter_count);
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