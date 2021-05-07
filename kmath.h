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

//三维向量
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

	//向量夹角（弧度）
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

	//右手系轴角顺时针旋转
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
    double *m_ptr;//数组指针
public:
    Matrix(Index_T r, Index_T c) :m_row(r), m_col(c)//非方阵构造
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
    Matrix(Index_T r, Index_T c, double val ) :m_row(r), m_col(c)// 赋初值val
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
    Matrix(Index_T n) :m_row(n), m_col(n)//方阵构造
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
    Matrix(const Matrix &rhs)//拷贝构造
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
 
    Matrix  &operator=(const Matrix&);  //如果类成员有指针必须重写赋值运算符，必须是成员
    friend istream &operator>>(istream&, Matrix&);
 
    friend ofstream &operator<<(ofstream &out, Matrix &obj);  // 输出到文件
    friend ostream &operator<<(ostream&, Matrix&);          // 输出到屏幕
    friend Matrix &operator<<(Matrix &mat, const double val);
    friend Matrix& operator,(Matrix &obj, const double val);
    friend Matrix  operator+(const Matrix&, const Matrix&);
    friend Matrix  operator-(const Matrix&, const Matrix&);
    friend Matrix  operator*(const Matrix&, const Matrix&);  //矩阵乘法
    friend Matrix  operator*(double, const Matrix&);  //矩阵乘法
    friend Matrix  operator*(const Matrix&, double);  //矩阵乘法
 
    friend Matrix  operator/(const Matrix&, double);  //矩阵 除以单数
 
    Matrix multi(const Matrix&); // 对应元素相乘
    Matrix mtanh(); // 对应元素相乘
    Index_T row()const{ return m_row; }
    Index_T col()const{ return m_col; }
    Matrix getrow(Index_T index); // 返回第index 行,索引从0 算起
    Matrix getcol(Index_T index); // 返回第index 列
 
    Matrix cov(_In_opt_ bool flag = true);   //协方差阵 或者样本方差
    double det();   //行列式
    Matrix solveAb(Matrix &obj);  // b是行向量或者列向量
    Matrix diag();  //返回对角线元素
    //Matrix asigndiag();  //对角线元素
    Matrix T()const;   //转置
    void sort(bool);//true为从小到大
    Matrix adjoint();
    Matrix inverse();
    void QR(_Out_ Matrix&, _Out_ Matrix&)const;
    Matrix eig_val(_In_opt_ Index_T _iters = 1000);
    Matrix eig_vect(_In_opt_ Index_T _iters = 1000);
 
    double norm1();//1范数
    double norm2();//2范数
    double mean();// 矩阵均值
    double*operator[](Index_T i){ return m_ptr + i*m_col; }//注意this加括号， (*this)[i][j]
    void zeromean(_In_opt_  bool flag = true);//默认参数为true计算列
    void normalize(_In_opt_  bool flag = true);//默认参数为true计算列
    Matrix exponent(double x);//每个元素x次幂
    Matrix  eye();//对角阵
    void  maxlimit(double max,double set=0);//对角阵
};
 

 