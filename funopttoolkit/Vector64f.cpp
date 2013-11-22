#include <cstring>

#define __VEC64F_EXPORT__
#include "Vector64f.h"
using namespace funopt;

#include "MTRand.h"
#include "funopt_macros.h"

Vector64f::Vector64f() :
    data(0)
{
}

Vector64f::Vector64f(int dim) :
    ndim(dim),
    data(0)
{
    data = new double[dim];
}

Vector64f::Vector64f(double* data_, int dim) :
    ndim(dim),
    data(0)
{
    data = new double[ndim];
    memcpy(data, data_, sizeof(double) * dim);
}

Vector64f::Vector64f(const Vector64f& v) :
    ndim(v.ndim),
    data(0)
{
    data = new double[ndim];
    memcpy(data, v.data, sizeof(double) * ndim);
}

Vector64f::~Vector64f()
{
    delete[] data;
}

Vector64f& Vector64f::operator=(const Vector64f& v)
{
	delete[] data;

    ndim = v.ndim;
    data = new double[ndim];
    memcpy(data, v.data, sizeof(double) * ndim);

    return *this;
}

inline double& Vector64f::operator()(int i) {
    return data[i];
}

inline double Vector64f::operator()(int i) const {
    return data[i];
}

Vector64f& Vector64f::operator+=(const Vector64f& v)
{
    for(int i=0; i<ndim; i++) {
        data[i] += v.data[i];
    }
    return *this;
}

Vector64f& Vector64f::operator-=(const Vector64f& v)
{
    for(int i=0; i<ndim; i++) {
        data[i] -= v.data[i];
    }
    return *this;
}

Vector64f operator+(const Vector64f& v, const Vector64f& u)
{
    Vector64f w = v;
    w += u;
    return w;
}

Vector64f operator-(const Vector64f& v, const Vector64f& u)
{
    Vector64f w = v;
    w -= u;
    return w;
}


Vector64f operator*(const Vector64f& v, double s)
{
	const int ndim = v.dim();
    Vector64f ret(v.dim());
    for(int i=0; i<ndim; i++) {
        ret(i) = v(i) * s;
    }
    return ret;
}

Vector64f operator*(double s, const Vector64f& v)
{
	return v * s;
}

Vector64f operator/(const Vector64f& v, double s)
{
    massert(s != 0.0, "zero division !");
	const int ndim = v.dim();
    Vector64f ret(ndim);
    for(int i=0; i<ndim; i++) {
        ret(i) = v(i) / s;
    }
    return ret;
}

Vector64f Vector64f::rand(int dim)
{
    Vector64f v(dim);
    MTRand rand;
    for(int i=0; i<dim; i++) {
        v(i) = rand.randReal();
    }
    return v;
}

Vector64f Vector64f::zeros(int dim)
{
    Vector64f v(dim);
    for(int i=0; i<dim; i++) {
        v(i) = 0.0;
    }
    return v;
}

int Vector64f::dim() const {
    return ndim;
}

double Vector64f::norm() const {
    return sqrt(norm2());
}

double Vector64f::norm2() const {
    double ret = 0.0;
    for(int i=0; i<ndim; i++) {
        ret += (*this)(i) * (*this)(i);
    }
    return ret;
}

double Vector64f::dot(const Vector64f& v) const
{
    massert(ndim == v.ndim, "vector dimension is different !");

    double ret = 0.0;
    for(int i=0; i<ndim; i++) {
        ret += data[i] * v.data[i];
    }
    return ret;
}

ostream& operator<<(ostream& os, const Vector64f& v)
{
	os << "[ ";
    for(int i=0; i<v.dim(); i++) {
        os << v(i);
        if(i != v.dim()-1) os << ", ";
    }
    os << " ]";
    return os;
}

