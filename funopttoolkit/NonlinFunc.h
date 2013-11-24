#ifndef _NONLIN_FUNC_H_
#define _NONLIN_FUNC_H_

#include "Vector64f.h"

namespace funopt {
    namespace nonlin {
        class funcNd {

        private:
            double (*_func)(const Vector64f&);
            double (*_grad)(const Vector64f&);

        public:
            funcNd(double (*func)(const Vector64f&), double (*grad)(const Vector64f&) = 0) :
                _func(func),
                _grad(grad)
            {
            }

            funcNd(const funcNd& f) :
                _func(f._func),
                _grad(f._grad)
            {
            }

            funcNd& operator=(const funcNd& f)
            {
                this->_func = f._func;
                this->_grad = f._grad;
                return *this;
            }

            double operator()(const Vector64f& x) const
            {
                return _func(x);
            }

            double grad(const Vector64f& x) const
            {
                return _grad(x);
            }
        };

        class func1dBase {
        public:
            func1dBase() {}
            virtual double operator()(double x) const = 0;
            double deriv(double x) const;
        };

        class func1dNormal : public func1dBase {
        private:
            double (*_func)(double);
            double (*_deriv)(double);

        public:
            func1dNormal(double (*func)(double), double (*deriv)(double) = 0) :
                func1dBase(),
                _func(func),
                _deriv(deriv)
            {
            }

            virtual double operator()(double x) const
            {
                return _func(x);
            }

            virtual double deriv(double x) const
            {
                return _deriv(x);
            }
        };

        class func1dlin : public func1dBase {
        private:
            Vector64f _x, _e;
            funcNd _funcNd;

        public:
            func1dlin(const funcNd& func, const Vector64f& x, const Vector64f& e) :
                func1dBase(),
                _x(x),
                _e(e),
                _funcNd(func)
            {
            }

            virtual double operator()(double x) const
            {
                return _funcNd(_x + x * _e);
            }
        };

        class func1d {
        private:
            func1dBase* fbase;

        public:
            func1d(double (*func)(double), double (*deriv)(double) = 0) :
                fbase(0)
            {
                fbase = new func1dNormal(func, deriv);
            }

            func1d(const funcNd& func, const Vector64f& x, const Vector64f& e) :
                fbase(0)
            {
                fbase = new func1dlin(func, x, e);
            }

            ~func1d() 
            {
                delete fbase;
            }

            double operator()(double x) const
            {
                return (*fbase)(x);
            }
        };       
    }
}

#endif
