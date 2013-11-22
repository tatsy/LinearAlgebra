#ifndef _NONLIN_FUNC_H_
#define _NONLIN_FUNC_H_

namespace funopt {
    class Vector64f;

    namespace nonlin {
        class func1d {
        private:
            double (*_func)(double);
            double (*_deriv)(double);

        public:
            func1d(double (*func)(double), double (*deriv)(double) = 0) :
                _func(func),
                _deriv(deriv)
            {
            }

            func1d(const func1d& f) :
                _func(f._func),
                _deriv(f._deriv)
            {
            }

            func1d& operator=(const func1d& f)
            {
                this->_func  = f._func;
                this->_deriv = f._deriv;
                return *this;
            }

            double operator()(double x) const
            {
                return _func(x);
            }

            double deriv(double x) const
            {
                return _deriv(x);
            }
        };

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

    }
}

#endif
