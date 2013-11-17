#include <cstring>

#define __EXPORT__
#include "Vector64f.h"
using namespace funopt;

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

double& Vector64f::operator()(int i) {
    return data[i];
}

double Vector64f::operator()(int i) const {
    return data[i];
}

Vector64f Vector64f::operator+(const Vector64f& v) const
{
    Vector64f ret(ndim);
    for(int i=0; i<ndim; i++) {
        ret.data[i] = data[i] + v.data[i];
    }
    return ret;
}

Vector64f Vector64f::operator-(const Vector64f& v) const 
{
    Vector64f ret(ndim);
    for(int i=0; i<ndim; i++) {
        ret.data[i] = data[i] - v.data[i];
    }
    return ret;
}

Vector64f Vector64f::operator*(double s) const
{
    Vector64f ret(ndim);
    for(int i=0; i<ndim; i++) {
        ret.data[i] = data[i] * s;
    }
    return ret;
}

Vector64f Vector64f::operator/(double s) const
{
    massert(s != 0.0, "zero division !");
    Vector64f ret(ndim);
    for(int i=0; i<ndim; i++) {
        ret.data[i] = data[i] / s;
    }
    return ret;
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
