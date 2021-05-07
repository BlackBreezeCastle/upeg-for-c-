#include "kmath.h"

double normalized_angle(double angle)
{
	angle = angle + 180;
	angle = fmod(angle, 360);
	angle = angle - 180;
	return angle;
}
double normalized_rad(double rad)
{
	rad = rad + PI;
	rad = fmod(rad, (2 * PI));
	rad = rad - PI;
	return rad;
}

Vector3::Vector3(void)
{
	m_x = 0.0;
	m_y = 0.0;
	m_z = 0.0;
}

Vector3::Vector3(const double &x, const double &y, const double &z)
{
	m_x = x;
	m_y = y;
	m_z = z;
}

inline double Vector3::x() const
{
	return m_x;
}

inline double Vector3::y() const
{
	return m_y;
}

inline double Vector3::z() const
{
	return m_z;
}

inline double Vector3::magnitude() const
{
	return pow((m_x*m_x + m_y*m_y + m_z*m_z), 0.5);
}

inline Vector3 Vector3::operator +(const Vector3 &b)const
{
	return Vector3(m_x + b.m_x, m_y + b.m_y, m_z + b.m_z);
}

inline Vector3 Vector3::operator-(const Vector3 &b)const
{
	return Vector3(m_x - b.m_x, m_y - b.m_y, m_z - b.m_z);
}

inline Vector3 operator*(const Vector3&a,const double &b)
{
	return Vector3(a.m_x*b, a.m_y*b, a.m_z*b);
}

inline Vector3 operator*(const double &b, const Vector3&a)
{
	return Vector3(a.m_x*b, a.m_y*b, a.m_z*b);
}

inline Vector3 Vector3::operator/(const double &b)const
{
	return Vector3(m_x/b, m_y/b, m_z/b);
}

bool Vector3 ::operator==(const Vector3 &b)const
{
	if (m_x == b.m_x
		&&m_y == b.m_y
		&&m_z == b.m_z)
	{
		return true;
	}
	return false;
}

Vector3 Vector3::normalized()const
{
	return (*this)*(magnitude() > 0 ? 1 / magnitude() : 0.0);
}

inline double Vector3::Dot(const Vector3 &a, const Vector3& b)
{
	return (a.m_x*b.m_x + a.m_y*b.m_y + a.m_z*b.m_z);
}

inline Vector3 Vector3::Cross(const Vector3&a, const Vector3 &b)
{
	double x = a.m_y*b.m_z - a.m_z*b.m_y;
	double y = a.m_z*b.m_x - a.m_x*b.m_z;
	double z = a.m_x*b.m_y - a.m_y*b.m_x;
	return Vector3(x, y, z);
}

inline double Vector3::Angle(const Vector3&a, const Vector3&b)
{
	return acos(std::max(std::min(Vector3::Dot(a, b) / (a.magnitude()*b.magnitude()), 1.0), -1.0));
}

Quaternion::Quaternion()
{
	m_w = -1.0;
	m_x = 0.0;
	m_y = 0.0;
	m_z = 0.0;
}

Quaternion::Quaternion(const double &w, const double &x, const double &y, const double &z)
{
	m_w = w;
	m_x = x;
	m_y = y;
	m_z = z;
}

Quaternion::Quaternion(const Vector3 &pivot, const double &rad)
{
	double theta = rad / 2;
	Vector3 u;
	if (pivot.magnitude()> 0.0)
	{
		u = pivot.normalized();
	}
	else
	{
		u = Vector3(1, 0, 0);
		theta = 0.0;
	}

	m_w = cos(theta);
	m_x = sin(theta)*u.x();
	m_y = sin(theta)*u.y();
	m_z = sin(theta)*u.z();
}



inline Quaternion Quaternion::operator*(const Quaternion &b)const
{
	double w1 = m_w;
	double w2 = b.m_w;
	Vector3 v1 = Vector3(m_x, m_y, m_z);
	Vector3 v2 = Vector3(b.m_x, b.m_y, b.m_z);
	double w3 = w1*w2 - Vector3::Dot(v1, v2);
	Vector3 v3 = Vector3::Cross(v1, v2) + w1*v2 + w2*v1;
	return Quaternion(w3, v3.x(), v3.y(), v3.z());
}

inline Quaternion Quaternion::inverse()const
{
	return Quaternion(m_w, -m_x, -m_y, -m_z);
}

inline Vector3 Quaternion::rotate(const Vector3 &v)const
{
	Vector3 u = Vector3(m_x, m_y, m_z);
	double s = m_w;
	return 2.0*Vector3::Dot(u, v)*u + (s*s - Vector3::Dot(u, u))*v + 2.0*s*Vector3::Cross(u, v);
}

/*
类方法的实现
*/
 
Matrix Matrix::mtanh() // 对应元素 tanh()
{
    Matrix ret(m_row, m_col);
    for (Index_T i = 0; i<ret.m_size; i++)
    {
        ret.m_ptr[i] = tanh(m_ptr[i]);
    }
    return ret;
}
/*
 * 递归调用
 */
double calcDet(Index_T n, double *&aa)
{
    if (n == 1)
        return aa[0];
    double *bb = new double[(n - 1)*(n - 1)];//创建n-1阶的代数余子式阵bb
    double sum = 0.0;
    for (Index_T Ai = 0; Ai<n; Ai++)
    {
        for (Index_T Bi = 0; Bi < n - 1; Bi++)//把aa阵第一列各元素的代数余子式存到bb
        {
            Index_T offset =  Bi < Ai ? 0 : 1; //bb中小于Ai的行，同行赋值，等于的错过，大于的加一
            for (Index_T j = 0; j<n - 1; j++)
            {
                bb[Bi*(n - 1) + j] = aa[(Bi + offset)*n + j + 1];
            }
        }
        int flag = (Ai % 2 == 0 ? 1 : -1);//因为列数为0，所以行数是偶数时候，代数余子式为1.
        sum += flag* aa[Ai*n] * calcDet(n - 1, bb);//aa第一列各元素与其代数余子式积的和即为行列式
    }
    delete[]bb;
    return sum;
}
 
Matrix Matrix::solveAb(Matrix &obj)
{
    Matrix ret(m_row, 1);
    if (m_size == 0 || obj.m_size == 0)
    {
        cout << "solveAb(Matrix &obj):this or obj is null" << endl;
        return ret;
    }
    if (m_row != obj.m_size)
    {
        cout << "solveAb(Matrix &obj):the row of two matrix is not equal!" << endl;
        return ret;
    }
 
    double *Dx = new double[m_row*m_row];
    for (Index_T i = 0; i<m_row; i++)
    {
        for (Index_T j = 0; j<m_row; j++)
        {
            Dx[i*m_row + j] = m_ptr[i*m_row + j];
        }
    }
    double D = calcDet(m_row, Dx);
    if (D == 0)
    {
        cout << "Cramer法则只能计算系数矩阵为满秩的矩阵" << endl;
        return  ret;
    }
 
    for (Index_T j = 0; j<m_row; j++)
    {
        for (Index_T i = 0; i<m_row; i++)
        {
            for (Index_T j = 0; j<m_row; j++)
            {
                Dx[i*m_row + j] = m_ptr[i*m_row + j];
            }
        }
        for (Index_T i = 0; i<m_row; i++)
        {
            Dx[i*m_row + j] = obj.m_ptr[i]; //obj赋值给第j列
        }
 
        //for( int i=0;i<m_row;i++) //print
        //{
        //    for(int j=0; j<m_row;j++)
        //    {
        //        cout<< Dx[i*m_row+j]<<"\t";
        //    }
        //    cout<<endl;
        //}
        ret[j][0] = calcDet(m_row, Dx) / D;
 
    }
 
    delete[]Dx;
    return ret;
}
 
Matrix Matrix::getrow(Index_T index)//返回行
{
    Matrix ret(1, m_col); //一行的返回值
 
    for (Index_T i = 0; i< m_col; i++)
    {
 
         ret[0][i] = m_ptr[(index) *m_col + i] ;
 
    }
    return ret;
}
 
Matrix Matrix::getcol(Index_T index)//返回列
{
    Matrix ret(m_row, 1); //一列的返回值
 
 
    for (Index_T i = 0; i< m_row; i++)
    {
 
        ret[i][0] = m_ptr[i *m_col + index];
 
    }
    return ret;
}
 
Matrix Matrix::exponent(double x)//每个元素x次幂
{
    Matrix ret(m_row, m_col);
    for (Index_T i = 0; i< m_row; i++)
    {
        for (Index_T j = 0; j < m_col; j++)
        {
            ret[i][j]= pow(m_ptr[i*m_col + j],x);
        }
    }
    return ret;
}
void Matrix::maxlimit(double max, double set)//每个元素x次幂
{
 
    for (Index_T i = 0; i< m_row; i++)
    {
        for (Index_T j = 0; j < m_col; j++)
        {
            m_ptr[i*m_col + j] = m_ptr[i*m_col + j]>max ? 0 : m_ptr[i*m_col + j];
        }
    }
 
}
Matrix Matrix::eye()//对角阵
{
 
    for (Index_T i = 0; i< m_row; i++)
    {
        for (Index_T j = 0; j < m_col; j++)
        {
            if (i == j)
            {
                m_ptr[i*m_col + j] = 1.0;
            }
        }
    }
    return *this;
}
void Matrix::zeromean(_In_opt_  bool flag)
{
    if (flag == true) //计算列均值
    {
        double *mean = new double[m_col];
        for (Index_T j = 0; j < m_col; j++)
        {
            mean[j] = 0.0;
            for (Index_T i = 0; i < m_row; i++)
            {
                mean[j] += m_ptr[i*m_col + j];
            }
            mean[j] /= m_row;
        }
        for (Index_T j = 0; j < m_col; j++)
        {
 
            for (Index_T i = 0; i < m_row; i++)
            {
                m_ptr[i*m_col + j] -= mean[j];
            }
        }
        delete[]mean;
    }
    else //计算行均值
    {
        double *mean = new double[m_row];
        for (Index_T i = 0; i< m_row; i++)
        {
            mean[i] = 0.0;
            for (Index_T j = 0; j < m_col; j++)
            {
                mean[i] += m_ptr[i*m_col + j];
            }
            mean[i] /= m_col;
        }
        for (Index_T i = 0; i < m_row; i++)
        {
            for (Index_T j = 0; j < m_col; j++)
            {
                m_ptr[i*m_col + j] -= mean[i];
            }
        }
        delete[]mean;
    }
}
 
void Matrix::normalize(_In_opt_  bool flag)
{
    if (flag == true) //计算列均值
    {
        double *mean = new double[m_col];
 
        for (Index_T j = 0; j < m_col; j++)
        {
            mean[j] = 0.0;
            for (Index_T i = 0; i < m_row; i++)
            {
                mean[j] += m_ptr[i*m_col + j];
            }
            mean[j] /= m_row;
        }
        for (Index_T j = 0; j < m_col; j++)
        {
 
            for (Index_T i = 0; i < m_row; i++)
            {
                m_ptr[i*m_col + j] -= mean[j];
            }
        }
        ///计算标准差
        for (Index_T j = 0; j < m_col; j++)
        {
            mean[j] = 0;
            for (Index_T i = 0; i < m_row; i++)
            {
                mean[j] += pow(m_ptr[i*m_col + j],2);//列平方和
            }
                mean[j] = sqrt(mean[j] / m_row); // 开方
        }
        for (Index_T j = 0; j < m_col; j++)
        {
            for (Index_T i = 0; i < m_row; i++)
            {
                m_ptr[i*m_col + j] /= mean[j];//列平方和
            }
        }
        delete[]mean;
    }
    else //计算行均值
    {
        double *mean = new double[m_row];
        for (Index_T i = 0; i< m_row; i++)
        {
            mean[i] = 0.0;
            for (Index_T j = 0; j < m_col; j++)
            {
                mean[i] += m_ptr[i*m_col + j];
            }
            mean[i] /= m_col;
        }
        for (Index_T i = 0; i < m_row; i++)
        {
            for (Index_T j = 0; j < m_col; j++)
            {
                m_ptr[i*m_col + j] -= mean[i];
            }
        }
        ///计算标准差
        for (Index_T i = 0; i< m_row; i++)
        {
            mean[i] = 0.0;
            for (Index_T j = 0; j < m_col; j++)
            {
                mean[i] += pow(m_ptr[i*m_col + j], 2);//列平方和
            }
            mean[i] = sqrt(mean[i] / m_col); // 开方
        }
        for (Index_T i = 0; i < m_row; i++)
        {
            for (Index_T j = 0; j < m_col; j++)
            {
                m_ptr[i*m_col + j] /= mean[i];
            }
        }
        delete[]mean;
    }
}
 
double Matrix::det()
{
    if (m_col == m_row)
        return calcDet(m_row, m_ptr);
    else
    {
        cout << ("行列不相等无法计算") << endl;
        return 0;
    }
}

istream& operator>>(istream &is, Matrix &obj)
{
    for (Index_T i = 0; i<obj.m_size; i++)
    {
        is >> obj.m_ptr[i];
    }
    return is;
}
 
ostream& operator<<(ostream &out, Matrix &obj)
{
    for (Index_T i = 0; i < obj.m_row; i++) //打印逆矩阵
    {
        for (Index_T j = 0; j < obj.m_col; j++)
        {
            out << (obj[i][j]) << "\t";
        }
        out << endl;
    }
    return out;
}
ofstream& operator<<(ofstream &out, Matrix &obj)//打印逆矩阵到文件
{
    for (Index_T i = 0; i < obj.m_row; i++)
    {
        for (Index_T j = 0; j < obj.m_col; j++)
        {
            out << (obj[i][j]) << "\t";
        }
        out << endl;
    }
    return out;
}
 
Matrix& operator<<(Matrix &obj, const double val)
{
    *obj.m_ptr = val;
    obj.m_curIndex = 1;
    return obj;
}
Matrix& operator,(Matrix &obj, const double val)
{
    if( obj.m_curIndex == 0 || obj.m_curIndex > obj.m_size - 1 )
    {
        return obj;
    }
    *(obj.m_ptr + obj.m_curIndex) = val;
    ++obj.m_curIndex;
    return obj;
}
 
Matrix operator+(const Matrix& lm, const Matrix& rm)
{
    if (lm.m_col != rm.m_col || lm.m_row != rm.m_row)
    {
        Matrix temp(0, 0);
        temp.m_ptr = NULL;
        cout << "operator+(): 矩阵shape 不合适,m_col:"
            << lm.m_col << "," << rm.m_col << ".  m_row:" << lm.m_row << ", " << rm.m_row << endl;
        return temp; //数据不合法时候，返回空矩阵
    }
    Matrix ret(lm.m_row, lm.m_col);
    for (Index_T i = 0; i<ret.m_size; i++)
    {
        ret.m_ptr[i] = lm.m_ptr[i] + rm.m_ptr[i];
    }
    return ret;
}
Matrix operator-(const Matrix& lm, const Matrix& rm)
{
    if (lm.m_col != rm.m_col || lm.m_row != rm.m_row)
    {
        Matrix temp(0, 0);
        temp.m_ptr = NULL;
        cout << "operator-(): 矩阵shape 不合适,m_col:"
            <<lm.m_col<<","<<rm.m_col<<".  m_row:"<< lm.m_row <<", "<< rm.m_row << endl;
 
        return temp; //数据不合法时候，返回空矩阵
    }
    Matrix ret(lm.m_row, lm.m_col);
    for (Index_T i = 0; i<ret.m_size; i++)
    {
        ret.m_ptr[i] = lm.m_ptr[i] - rm.m_ptr[i];
    }
    return ret;
}
Matrix operator*(const Matrix& lm, const Matrix& rm)  //矩阵乘法
{
    if (lm.m_size == 0 || rm.m_size == 0 || lm.m_col != rm.m_row)
    {
        Matrix temp(0, 0);
        temp.m_ptr = NULL;
        cout << "operator*(): 矩阵shape 不合适,m_col:"
            << lm.m_col << "," << rm.m_col << ".  m_row:" << lm.m_row << ", " << rm.m_row << endl;
        return temp; //数据不合法时候，返回空矩阵
    }
    Matrix ret(lm.m_row, rm.m_col);
    for (Index_T i = 0; i<lm.m_row; i++)
    {
        for (Index_T j = 0; j< rm.m_col; j++)
        {
            for (Index_T k = 0; k< lm.m_col; k++)//lm.m_col == rm.m_row
            {
                ret.m_ptr[i*rm.m_col + j] += lm.m_ptr[i*lm.m_col + k] * rm.m_ptr[k*rm.m_col + j];
            }
        }
    }
    return ret;
}
Matrix operator*(double val, const Matrix& rm)  //矩阵乘 单数
{
    Matrix ret(rm.m_row, rm.m_col);
    for (Index_T i = 0; i<ret.m_size; i++)
    {
        ret.m_ptr[i] = val * rm.m_ptr[i];
    }
    return ret;
}
Matrix operator*(const Matrix&lm, double val)  //矩阵乘 单数
{
    Matrix ret(lm.m_row, lm.m_col);
    for (Index_T i = 0; i<ret.m_size; i++)
    {
        ret.m_ptr[i] = val * lm.m_ptr[i];
    }
    return ret;
}
 
Matrix operator/(const Matrix&lm, double val)  //矩阵除以 单数
{
    Matrix ret(lm.m_row, lm.m_col);
    for (Index_T i = 0; i<ret.m_size; i++)
    {
        ret.m_ptr[i] =  lm.m_ptr[i]/val;
    }
    return ret;
}
Matrix Matrix::multi(const Matrix&rm)// 对应元素相乘
{
    if (m_col != rm.m_col || m_row != rm.m_row)
    {
        Matrix temp(0, 0);
        temp.m_ptr = NULL;
        cout << "multi(const Matrix&rm): 矩阵shape 不合适,m_col:"
            << m_col << "," << rm.m_col << ".  m_row:" << m_row << ", " << rm.m_row << endl;
        return temp; //数据不合法时候，返回空矩阵
    }
    Matrix ret(m_row,m_col);
    for (Index_T i = 0; i<ret.m_size; i++)
    {
        ret.m_ptr[i] = m_ptr[i] * rm.m_ptr[i];
    }
    return ret;
 
}
 
Matrix&  Matrix::operator=(const Matrix& rhs)
{
    if (this != &rhs)
    {
        m_row = rhs.m_row;
        m_col = rhs.m_col;
        m_size = rhs.m_size;
        if (m_ptr != NULL)
            delete[] m_ptr;
        m_ptr = new double[m_size];
        for (Index_T i = 0; i<m_size; i++)
        {
            m_ptr[i] = rhs.m_ptr[i];
        }
    }
    return *this;
}
//||matrix||_2  求A矩阵的2范数
double Matrix::norm2()
{
    double norm = 0;
    for (Index_T i = 0; i < m_size; ++i)
    {
        norm += m_ptr[i] * m_ptr[i];
    }
    return (double)sqrt(norm);
}
double Matrix::norm1()
{
    double sum = 0;
    for (Index_T i = 0; i < m_size; ++i)
    {
        sum += abs(m_ptr[i]);
    }
    return sum;
}
double Matrix::mean()
{
    double sum = 0;
    for (Index_T i = 0; i < m_size; ++i)
    {
        sum += (m_ptr[i]);
    }
    return sum/m_size;
}
 
 
 
 
void Matrix::sort(bool flag)
{
    double tem;
    for (Index_T i = 0; i<m_size; i++)
    {
        for (Index_T j = i + 1; j<m_size; j++)
        {
            if (flag == true)
            {
                if (m_ptr[i]>m_ptr[j])
                {
                    tem = m_ptr[i];
                    m_ptr[i] = m_ptr[j];
                    m_ptr[j] = tem;
                }
            }
            else
            {
                if (m_ptr[i]<m_ptr[j])
                {
                    tem = m_ptr[i];
                    m_ptr[i] = m_ptr[j];
                    m_ptr[j] = tem;
                }
            }
 
        }
    }
}
Matrix Matrix::diag()
{
    if (m_row != m_col)
    {
        Matrix m(0);
        cout << "diag():m_row != m_col" << endl;
        return m;
    }
    Matrix m(m_row);
    for (Index_T i = 0; i<m_row; i++)
    {
        m.m_ptr[i*m_row + i] = m_ptr[i*m_row + i];
    }
    return m;
}
Matrix Matrix::T()const
{
    Matrix tem(m_col, m_row);
    for (Index_T i = 0; i<m_row; i++)
    {
        for (Index_T j = 0; j<m_col; j++)
        {
            tem[j][i] = m_ptr[i*m_col + j];// (*this)[i][j]
        }
    }
    return tem;
}
void  Matrix::QR(Matrix &Q, Matrix &R) const
{
    //如果A不是一个二维方阵，则提示错误，函数计算结束
    if (m_row != m_col)
    {
        printf("ERROE: QR() parameter A is not a square matrix!\n");
        return;
    }
    const Index_T N = m_row;
    double *a = new double[N];
    double *b = new double[N];
 
    for (Index_T j = 0; j < N; ++j)  //(Gram-Schmidt) 正交化方法
    {
        for (Index_T i = 0; i < N; ++i)  //第j列的数据存到a，b
            a[i] = b[i] = m_ptr[i * N + j];
 
        for (Index_T i = 0; i<j; ++i)  //第j列之前的列
        {
            R.m_ptr[i * N + j] = 0;  //
            for (Index_T m = 0; m < N; ++m)
            {
                R.m_ptr[i * N + j] += a[m] * Q.m_ptr[m *N + i]; //R[i,j]值为Q第i列与A的j列的内积
            }
            for (Index_T m = 0; m < N; ++m)
            {
                b[m] -= R.m_ptr[i * N + j] * Q.m_ptr[m * N + i]; //
            }
        }
 
        double norm = 0;
        for (Index_T i = 0; i < N; ++i)
        {
            norm += b[i] * b[i];
        }
        norm = (double)sqrt(norm);
 
        R.m_ptr[j*N + j] = norm; //向量b[]的2范数存到R[j,j]
 
        for (Index_T i = 0; i < N; ++i)
        {
            Q.m_ptr[i * N + j] = b[i] / norm; //Q 阵的第j列为单位化的b[]
        }
    }
    delete[]a;
    delete[]b;
}
Matrix Matrix::eig_val(_In_opt_ Index_T _iters)
{
    if (m_size == 0 || m_row != m_col)
    {
        cout << "矩阵为空或者非方阵！" << endl;
        Matrix rets(0);
        return rets;
    }
    //if (det() == 0)
    //{
    //  cout << "非满秩矩阵没法用QR分解计算特征值！" << endl;
    //  Matrix rets(0);
    //  return rets;
    //}
    const Index_T N = m_row;
    Matrix matcopy(*this);//备份矩阵
    Matrix Q(N), R(N);
    /*当迭代次数足够多时,A 趋于上三角矩阵，上三角矩阵的对角元就是A的全部特征值。*/
    for (Index_T k = 0; k < _iters; ++k)
    {
        //cout<<"this:\n"<<*this<<endl;
        QR(Q, R);
        *this = R*Q;
        /*  cout<<"Q:\n"<<Q<<endl;
        cout<<"R:\n"<<R<<endl;  */
    }
    Matrix val = diag();
    *this = matcopy;//恢复原始矩阵；
    return val;
}
Matrix Matrix::eig_vect(_In_opt_ Index_T _iters)
{
    if (m_size == 0 || m_row != m_col)
    {
        cout << "矩阵为空或者非方阵！" << endl;
        Matrix rets(0);
        return rets;
    }
    if (det() == 0)
    {
      cout << "非满秩矩阵没法用QR分解计算特征向量！" << endl;
      Matrix rets(0);
      return rets;
    }
    Matrix matcopy(*this);//备份矩阵
    Matrix eigenValue = eig_val(_iters);
    Matrix ret(m_row);
    const Index_T NUM = m_col;
    double eValue;
    double sum, midSum, diag;
    Matrix copym(*this);
    for (Index_T count = 0; count < NUM; ++count)
    {
        //计算特征值为eValue，求解特征向量时的系数矩阵
        *this = copym;
        eValue = eigenValue[count][count];
 
        for (Index_T i = 0; i < m_col; ++i)//A-lambda*I
        {
            m_ptr[i * m_col + i] -= eValue;
        }
        //cout<<*this<<endl;
        //将 this为阶梯型的上三角矩阵
        for (Index_T i = 0; i < m_row - 1; ++i)
        {
            diag = m_ptr[i*m_col + i];  //提取对角元素
            for (Index_T j = i; j < m_col; ++j)
            {
                m_ptr[i*m_col + j] /= diag; //【i,i】元素变为1
            }
            for (Index_T j = i + 1; j<m_row; ++j)
            {
                diag = m_ptr[j *  m_col + i];
                for (Index_T q = i; q < m_col; ++q)//消去第i+1行的第i个元素
                {
                    m_ptr[j*m_col + q] -= diag*m_ptr[i*m_col + q];
                }
            }
        }
        //cout<<*this<<endl;
        //特征向量最后一行元素置为1
        midSum = ret.m_ptr[(ret.m_row - 1) * ret.m_col + count] = 1;
        for (int m = m_row - 2; m >= 0; --m)
        {
            sum = 0;
            for (Index_T j = m + 1; j < m_col; ++j)
            {
                sum += m_ptr[m *  m_col + j] * ret.m_ptr[j * ret.m_col + count];
            }
            sum = -sum / m_ptr[m *  m_col + m];
            midSum += sum * sum;
            ret.m_ptr[m * ret.m_col + count] = sum;
        }
        midSum = sqrt(midSum);
        for (Index_T i = 0; i < ret.m_row; ++i)
        {
            ret.m_ptr[i * ret.m_col + count] /= midSum; //每次求出一个列向量
        }
    }
    *this = matcopy;//恢复原始矩阵；
    return ret;
}
Matrix Matrix::cov(bool flag)
{
    //m_row 样本数，column 变量数
    if (m_col == 0)
    {
        Matrix m(0);
        return m;
    }
    double *mean = new double[m_col]; //均值向量
 
    for (Index_T j = 0; j<m_col; j++) //init
    {
        mean[j] = 0.0;
    }
    Matrix ret(m_col);
    for (Index_T j = 0; j<m_col; j++) //mean
    {
        for (Index_T i = 0; i<m_row; i++)
        {
            mean[j] += m_ptr[i*m_col + j];
        }
        mean[j] /= m_row;
    }
    Index_T i, k, j;
    for (i = 0; i<m_col; i++) //第一个变量
    {
        for (j = i; j<m_col; j++) //第二个变量
        {
            for (k = 0; k<m_row; k++) //计算
            {
                ret[i][j] += (m_ptr[k*m_col + i] - mean[i])*(m_ptr[k*m_col + j] - mean[j]);
 
            }
            if (flag == true)
            {
                ret[i][j] /= (m_row-1);
            }
            else
            {
                ret[i][j] /= (m_row);
            }
        }
    }
    for (i = 0; i<m_col; i++) //补全对应面
    {
        for (j = 0; j<i; j++)
        {
            ret[i][j] = ret[j][i];
        }
    }
    return ret;
}
 
/*
 * 返回代数余子式
 */
double CalcAlgebraicCofactor( Matrix& srcMat, Index_T ai, Index_T aj)
{
    Index_T temMatLen = srcMat.row()-1;
    Matrix temMat(temMatLen);
    for (Index_T bi = 0; bi < temMatLen; bi++)
    {
        for (Index_T bj = 0; bj < temMatLen; bj++)
        {
            Index_T rowOffset = bi < ai ? 0 : 1;
            Index_T colOffset = bj < aj ? 0 : 1;
            temMat[bi][bj] = srcMat[bi + rowOffset][bj + colOffset];
        }
    }
    int flag = (ai + aj) % 2 == 0 ? 1 : -1;
    return flag * temMat.det();
}
 
/*
 * 返回伴随阵
 */
Matrix Matrix::adjoint()
{
    if (m_row != m_col)
    {
        return Matrix(0);
    }
 
    Matrix adjointMat(m_row);
    for (Index_T ai = 0; ai < m_row; ai++)
    {
        for (Index_T aj = 0; aj < m_row; aj++)
        {
            adjointMat.m_ptr[aj*m_row + ai] = CalcAlgebraicCofactor(*this, ai, aj);
        }
    }
    return adjointMat;
}
 
Matrix Matrix::inverse()
{
    double detOfMat = det();
    if (detOfMat == 0)
    {
        cout << "行列式为0，不能计算逆矩阵。" << endl;
        return Matrix(0);
    }
    return adjoint()/detOfMat;
}