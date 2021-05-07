#pragma once
#include<math.h>
#include<algorithm>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cmath>
using namespace std;
#define PI 3.14159265358979323846
#define INFINITY_DOUBLE 1e20
double normalized_angle(double angle);

double normalized_rad(double rad);

//��ά����
class _declspec(dllexport) Vector3
{
private:
	double m_x;
	double m_y;
	double m_z;

public:
	Vector3();

	Vector3(const double &x, const double& y, const double& z);

	double x()const;

	double y()const;

	double z()const;

	double magnitude()const;

	Vector3 operator +(const Vector3 &b)const;

	Vector3 operator -(const Vector3 &b)const;

	friend Vector3 operator *(const Vector3 &a ,const double &b);
	friend Vector3 operator *(const double &b, const Vector3 &a);

	Vector3 operator /(const double &b)const;

	bool operator==(const Vector3 &b)const;

	Vector3 normalized()const;

public:
	double static Dot(const Vector3&a,const Vector3&b);

	Vector3 static Cross(const Vector3&a, const Vector3&b);

	//�����нǣ����ȣ�
	double static Angle(const Vector3&a, const Vector3&b);
};

class _declspec(dllexport)Quaternion
{
private:
	double m_w;

	double m_x;
	double m_y;
	double m_z;
public:
	Quaternion();

	Quaternion(const double &w,const double &x, const double &y, const double &z);

	Quaternion(const Vector3 &pivot, const double &rad);

	Quaternion operator *(const Quaternion &b)const;

	Quaternion inverse()const;

	//����ϵ���˳ʱ����ת
	Vector3 rotate(const Vector3 &v)const;
};

#ifndef _In_opt_
  #define _In_opt_
#endif
#ifndef _Out_
  #define _Out_
#endif
typedef unsigned Index_T;
class Matrix
{
private:
    Index_T m_row, m_col;
    Index_T m_size;
    Index_T m_curIndex;
    double *m_ptr;//����ָ��
public:
    Matrix(Index_T r, Index_T c) :m_row(r), m_col(c)//�Ƿ�����
    {
        m_size = r*c;
        if (m_size>0)
        {
            m_ptr = new double[m_size];
			memset(m_ptr,0,sizeof(double)*m_size);
        }
        else
            m_ptr = NULL;
    };
    Matrix(Index_T r, Index_T c, double val ) :m_row(r), m_col(c)// ����ֵval
    {
        m_size = r*c;
        if (m_size>0)
        {
            m_ptr = new double[m_size];
			memset(m_ptr,0,sizeof(double)*m_size);			
        }
        else
            m_ptr = NULL;
    };
    Matrix(Index_T n) :m_row(n), m_col(n)//������
    {
        m_size = n*n;
        if (m_size>0)
        {
            m_ptr = new double[m_size];
			memset(m_ptr,0,sizeof(double)*m_size);
        }
        else
            m_ptr = NULL;
    };
    Matrix(const Matrix &rhs)//��������
    {
        m_row = rhs.m_row;
        m_col = rhs.m_col;
        m_size = rhs.m_size;
        m_ptr = new double[m_size];
        for (Index_T i = 0; i<m_size; i++)
            m_ptr[i] = rhs.m_ptr[i];
    }
 
    ~Matrix()
    {
        if (m_ptr != NULL)
        {
            delete[]m_ptr;
            m_ptr = NULL;
        }
    }
 
    Matrix  &operator=(const Matrix&);  //������Ա��ָ�������д��ֵ������������ǳ�Ա
    friend istream &operator>>(istream&, Matrix&);
 
    friend ofstream &operator<<(ofstream &out, Matrix &obj);  // ������ļ�
    friend ostream &operator<<(ostream&, Matrix&);          // �������Ļ
    friend Matrix &operator<<(Matrix &mat, const double val);
    friend Matrix& operator,(Matrix &obj, const double val);
    friend Matrix  operator+(const Matrix&, const Matrix&);
    friend Matrix  operator-(const Matrix&, const Matrix&);
    friend Matrix  operator*(const Matrix&, const Matrix&);  //����˷�
    friend Matrix  operator*(double, const Matrix&);  //����˷�
    friend Matrix  operator*(const Matrix&, double);  //����˷�
 
    friend Matrix  operator/(const Matrix&, double);  //���� ���Ե���
 
    Matrix multi(const Matrix&); // ��ӦԪ�����
    Matrix mtanh(); // ��ӦԪ�����
    Index_T row()const{ return m_row; }
    Index_T col()const{ return m_col; }
    Matrix getrow(Index_T index); // ���ص�index ��,������0 ����
    Matrix getcol(Index_T index); // ���ص�index ��
 
    Matrix cov(_In_opt_ bool flag = true);   //Э������ ������������
    double det();   //����ʽ
    Matrix solveAb(Matrix &obj);  // b������������������
    Matrix diag();  //���ضԽ���Ԫ��
    //Matrix asigndiag();  //�Խ���Ԫ��
    Matrix T()const;   //ת��
    void sort(bool);//trueΪ��С����
    Matrix adjoint();
    Matrix inverse();
    void QR(_Out_ Matrix&, _Out_ Matrix&)const;
    Matrix eig_val(_In_opt_ Index_T _iters = 1000);
    Matrix eig_vect(_In_opt_ Index_T _iters = 1000);
 
    double norm1();//1����
    double norm2();//2����
    double mean();// �����ֵ
    double*operator[](Index_T i){ return m_ptr + i*m_col; }//ע��this�����ţ� (*this)[i][j]
    void zeromean(_In_opt_  bool flag = true);//Ĭ�ϲ���Ϊtrue������
    void normalize(_In_opt_  bool flag = true);//Ĭ�ϲ���Ϊtrue������
    Matrix exponent(double x);//ÿ��Ԫ��x����
    Matrix  eye();//�Խ���
    void  maxlimit(double max,double set=0);//�Խ���
};
 

 