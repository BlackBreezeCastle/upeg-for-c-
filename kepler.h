#pragma once
#include "kmath.h"
//#include "kmath.cpp"
#include<stdio.h>
/* ---------------------------------------------------------------------
** ��������: E_to_F
** ��������: ��ƫ���Ǽ��������
** �������: ecc -> ƫ����             [N/A]
**           E   -> ƫ����             [rad]
** ����ֵ  : �����                    [rad]
** ---------------------------------------------------------------------*/
double E_to_F(double ecc, double E);

//double E2Fstd(double ecc, double E);

//double E2Fhyp(double eec, double E);

double F_to_E(double ecc, double F);

//double F2Estd(double ecc, double F);

//double F2Ehyp(double ecc, double F);

double F_to_R(double ecc, double sem, double f);

double R_to_F(double ecc, double sem, double r);

double E_to_M(double ecc, double E);


double M_to_E(double ecc, double M);

