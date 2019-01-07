/*
 * Class: hyperdual
 *
 * Implementation of hyper-dual numbers
 *
 * Written by: Jeffrey A. Fike
 * Stanford University, Department of Aeronautics and Astronautics
 *
 * Copyright (c) 2006 Jeffrey A. Fike
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 */

#ifndef _hyperdual_h
#define _hyperdual_h

#include <iostream>
#include <math.h>
using namespace std;

class hyperdual{
public:
    double f0,f1,f2,f12;
    //creation operators and function to manually set values
    hyperdual();
    hyperdual(double x1,double x2,double x3,double x4);
    hyperdual(double x1);
    void setvalues(double x1,double x2,double x3,double x4);

    //examine values
    void view(void);
    double real(void);
    double eps1(void);
    double eps2(void);
    double eps1eps2(void);
    friend ostream& operator<<(ostream& output, const hyperdual& rhs);

    //basic manipulation
    hyperdual operator+ () const;
    hyperdual operator+ (const hyperdual rhs) const;
    friend hyperdual operator+ (const double lhs, const hyperdual rhs);
    hyperdual operator- () const;
    hyperdual operator- (const hyperdual rhs) const;
    friend hyperdual operator- (const double lhs, const hyperdual rhs);
    hyperdual operator* (const hyperdual rhs)const;
    friend hyperdual operator* (const double lhs, const hyperdual rhs);
    friend hyperdual operator/ (const hyperdual lhs, const hyperdual rhs);
    friend hyperdual operator/ (const double lhs, const hyperdual rhs);
    friend hyperdual operator/ (const hyperdual lhs, const double rhs);
    hyperdual& operator+= (hyperdual rhs);
    hyperdual& operator-= (hyperdual rhs);
    hyperdual& operator*= (hyperdual rhs);
    hyperdual& operator*= (double rhs);
    hyperdual& operator/= (double rhs);

    //math.h functions
    friend hyperdual pow (hyperdual x, double a);
    friend hyperdual pow (hyperdual x, hyperdual a);
    friend hyperdual exp(hyperdual x);
    friend hyperdual log(hyperdual x);
    friend hyperdual sin(hyperdual x);
    friend hyperdual cos(hyperdual x);
    friend hyperdual tan(hyperdual x);
    friend hyperdual asin(hyperdual x);
    friend hyperdual acos(hyperdual x);
    friend hyperdual atan(hyperdual x);
    friend hyperdual sqrt(hyperdual x);
    friend hyperdual fabs(hyperdual x);
    friend hyperdual max(hyperdual x1, hyperdual x2);
    friend hyperdual max(hyperdual x1, double x2);
    friend hyperdual max(double x1, hyperdual x2);
    friend hyperdual min(hyperdual x1, hyperdual x2);
    friend hyperdual min(hyperdual x1, double x2);
    friend hyperdual min(double x1, hyperdual x2);

    //comparisons
    friend bool operator> (hyperdual lhs, hyperdual rhs);
    friend bool operator> (double lhs, hyperdual rhs);
    friend bool operator> (hyperdual lhs, double rhs);
    friend bool operator>= (hyperdual lhs, hyperdual rhs);
    friend bool operator>= (double lhs, hyperdual rhs);
    friend bool operator>= (hyperdual lhs, double rhs);
    friend bool operator< (hyperdual lhs, hyperdual rhs);
    friend bool operator< (double lhs, hyperdual rhs);
    friend bool operator< (hyperdual lhs, double rhs);
    friend bool operator<= (hyperdual lhs, hyperdual rhs);
    friend bool operator<= (double lhs, hyperdual rhs);
    friend bool operator<= (hyperdual lhs, double rhs);
    friend bool operator== (hyperdual lhs, hyperdual rhs);
    friend bool operator== (double lhs, hyperdual rhs);
    friend bool operator== (hyperdual lhs, double rhs);
    friend bool operator!= (hyperdual lhs, hyperdual rhs);
    friend bool operator!= (double lhs, hyperdual rhs);
    friend bool operator!= (hyperdual lhs, double rhs);
};


hyperdual::hyperdual()
{
    f0 = 0.0;
    f1 = 0.0;
    f2 = 0.0;
    f12 = 0.0;
}
hyperdual::hyperdual(double x1,double x2,double x3,double x4)
{
    f0 = x1;
    f1 = x2;
    f2 = x3;
    f12 = x4;
}
hyperdual::hyperdual(double x1)
{
    f0 = x1;
    f1 = 0.0;
    f2 = 0.0;
    f12 = 0.0;
}
void hyperdual::setvalues(double x1,double x2,double x3,double x4)
{
    f0 = x1;
    f1 = x2;
    f2 = x3;
    f12 = x4;
}

void hyperdual::view(void)
{
    printf("%g  +  %g epsilon1  +  %g epsilon2  +  %g epsilon1 epsilon2\n",f0,f1,f2,f12);
}
double hyperdual::real(void)
{
    return f0;
}
double hyperdual::eps1(void)
{
    return f1;
}
double hyperdual::eps2(void)
{
    return f2;
}
double hyperdual::eps1eps2(void)
{
    return f12;
}
ostream& operator<<(ostream& output, const hyperdual& rhs)
{
    output << "(" << rhs.f0 << ","<< rhs.f1 << ","<< rhs.f2 << ","<< rhs.f12 << ")";
    return output;
}

hyperdual hyperdual::operator+ () const
{
    return *this;
}
hyperdual hyperdual::operator+ (const hyperdual rhs) const
{
    hyperdual temp;
    temp.f0 = f0 + rhs.f0;
    temp.f1 = f1 + rhs.f1;
    temp.f2 = f2 + rhs.f2;
    temp.f12 = f12 + rhs.f12;
    return temp;
}
hyperdual operator+ (const double lhs, const hyperdual rhs)
{
    hyperdual temp;
    temp.f0 = lhs + rhs.f0;
    temp.f1 = rhs.f1;
    temp.f2 = rhs.f2;
    temp.f12 = rhs.f12;
    return temp;
}
hyperdual hyperdual::operator- () const
{
    hyperdual temp;
    temp.f0 = -f0;
    temp.f1 = -f1;
    temp.f2 = -f2;
    temp.f12 = -f12;
    return temp;
}
hyperdual hyperdual::operator- (const hyperdual rhs) const
{
    hyperdual temp;
    temp.f0 = f0 - rhs.f0;
    temp.f1 = f1 - rhs.f1;
    temp.f2 = f2 - rhs.f2;
    temp.f12 = f12 - rhs.f12;
    return temp;
}
hyperdual operator- (const double lhs, const hyperdual rhs)
{
    hyperdual temp;
    temp.f0 = lhs - rhs.f0;
    temp.f1 = -rhs.f1;
    temp.f2 = -rhs.f2;
    temp.f12 = -rhs.f12;
    return temp;
}
hyperdual hyperdual::operator* (const hyperdual rhs) const
{
    hyperdual temp;
    temp.f0 = f0*rhs.f0;
    temp.f1 = f0*rhs.f1 + f1*rhs.f0;
    temp.f2 = f0*rhs.f2 + f2*rhs.f0;
    temp.f12 = f0*rhs.f12 + f1*rhs.f2 + f2*rhs.f1 + f12*rhs.f0;
    return temp;
}
hyperdual operator* (const double lhs, const hyperdual rhs)
{
    hyperdual temp;
    temp.f0 = lhs*rhs.f0;
    temp.f1 = lhs*rhs.f1;
    temp.f2 = lhs*rhs.f2;
    temp.f12 = lhs*rhs.f12;
    return temp;
}
hyperdual operator/ (const hyperdual lhs, const hyperdual rhs)
{
    hyperdual temp,inv;
    inv = pow(rhs,-1);
    temp = lhs*inv;
    return temp;
}
hyperdual operator/ (const double lhs, const hyperdual rhs)
{
    hyperdual temp,inv;
    inv = pow(rhs,-1);
    temp = lhs*inv;
    return temp;
}
hyperdual operator/ (const hyperdual lhs, const double rhs)
{
    hyperdual temp;
    double inv;
    inv = 1.0/rhs;
    temp.f0 = inv*lhs.f0;
    temp.f1 = inv*lhs.f1;
    temp.f2 = inv*lhs.f2;
    temp.f12 = inv*lhs.f12;
    return temp;
}
hyperdual& hyperdual::operator+= (hyperdual rhs)
{
    f0 += rhs.f0;
    f1 += rhs.f1;
    f2 += rhs.f2;
    f12 += rhs.f12;
    return *this;
}
hyperdual& hyperdual::operator-= (hyperdual rhs)
{
    f0 -= rhs.f0;
    f1 -= rhs.f1;
    f2 -= rhs.f2;
    f12 -= rhs.f12;
    return *this;
}
hyperdual& hyperdual::operator*= (hyperdual rhs)
{
    double tf0,tf1,tf2,tf12;
    tf0=f0;
    tf1=f1;
    tf2=f2;
    tf12=f12;
    f0 = tf0*rhs.f0;
    f1 = tf0*rhs.f1 + tf1*rhs.f0;
    f2 = tf0*rhs.f2 + tf2*rhs.f0;
    f12 = tf0*rhs.f12 + tf1*rhs.f2 + tf2*rhs.f1 + tf12*rhs.f0;
    return *this;
}
hyperdual& hyperdual::operator*= (double rhs)
{
    f0 *= rhs;
    f1 *= rhs;
    f2 *= rhs;
    f12 *= rhs;
    return *this;
}
hyperdual& hyperdual::operator/= (double rhs)
{
    f0 /= rhs;
    f1 /= rhs;
    f2 /= rhs;
    f12 /= rhs;
    return *this;
}
hyperdual pow (hyperdual x, double a)
{
    hyperdual temp;
    double deriv,xval,tol;
    xval = x.f0;
    tol = 1e-15;
    if (fabs(xval) < tol)
    {
        if (xval >= 0)
            xval = tol;
        if (xval < 0)
            xval = -tol;
    }
    deriv = a*pow(xval,(a-1));
    //temp.f0 = pow(xval,a);
    temp.f0 = pow(x.f0,a);  //Use actual x value, only use tol for derivs
    temp.f1 = x.f1*deriv;
    temp.f2 = x.f2*deriv;
    temp.f12 = x.f12*deriv + a*(a-1)*x.f1*x.f2*pow(xval,(a-2));

    return temp;
}
hyperdual pow (hyperdual x, hyperdual a)
{
    return exp(a*log(x));
}
hyperdual exp(hyperdual x)
{
    hyperdual temp;
    double deriv;
    deriv = exp(x.f0);
    temp.f0 = deriv;
    temp.f1 = deriv*x.f1;
    temp.f2 = deriv*x.f2;
    temp.f12 = deriv*(x.f12 + x.f1*x.f2);
    return temp;
}
hyperdual log(hyperdual x)
{
    hyperdual temp;
    double deriv1,deriv2;
    deriv1 = x.f1/x.f0;
    deriv2 = x.f2/x.f0;
    temp.f0 = log(x.f0);
    temp.f1 = deriv1;
    temp.f2 = deriv2;
    temp.f12 = x.f12/x.f0 - (deriv1*deriv2);
    return temp;
}
hyperdual sin(hyperdual x)
{
    hyperdual temp;
    double funval,deriv;
    funval = sin(x.f0);
    deriv = cos(x.f0);
    temp.f0 = funval;
    temp.f1 = deriv*x.f1;
    temp.f2 = deriv*x.f2;
    temp.f12 = deriv*x.f12 - funval*x.f1*x.f2;
    return temp;
}
hyperdual cos(hyperdual x)
{
    hyperdual temp;
    double funval,deriv;
    funval = cos(x.f0);
    deriv = -sin(x.f0);
    temp.f0 = funval;
    temp.f1 = deriv*x.f1;
    temp.f2 = deriv*x.f2;
    temp.f12 = deriv*x.f12 - funval*x.f1*x.f2;
    return temp;
}
hyperdual tan(hyperdual x)
{
    hyperdual temp;
    double funval,deriv;
    funval = tan(x.f0);
    deriv  = funval*funval + 1.0;
    temp.f0 = funval;
    temp.f1 = deriv*x.f1;
    temp.f2 = deriv*x.f2;
    temp.f12 = deriv*x.f12 + x.f1*x.f2*(2*funval*deriv);
    return temp;
}
hyperdual asin(hyperdual x)
{
    hyperdual temp;
    double funval,deriv1,deriv;
    funval = asin(x.f0);
    deriv1 = 1.0-x.f0*x.f0;
    deriv = 1.0/sqrt(deriv1);
    temp.f0 = funval;
    temp.f1 = deriv*x.f1;
    temp.f2 = deriv*x.f2;
    temp.f12 = deriv*x.f12 + x.f1*x.f2*(x.f0*pow(deriv1,-1.5));
    return temp;
}
hyperdual acos(hyperdual x)
{
    hyperdual temp;
    double funval,deriv1,deriv;
    funval = acos(x.f0);
    deriv1 = 1.0-x.f0*x.f0;
    deriv = -1.0/sqrt(deriv1);
    temp.f0 = funval;
    temp.f1 = deriv*x.f1;
    temp.f2 = deriv*x.f2;
    temp.f12 = deriv*x.f12 + x.f1*x.f2*(-x.f0*pow(deriv1,-1.5));
    return temp;
}
hyperdual atan(hyperdual x)
{
    hyperdual temp;
    double funval,deriv1,deriv;
    funval = atan(x.f0);
    deriv1 = 1.0+x.f0*x.f0;
    deriv = 1.0/deriv1;
    temp.f0 = funval;
    temp.f1 = deriv*x.f1;
    temp.f2 = deriv*x.f2;
    temp.f12 = deriv*x.f12 + x.f1*x.f2*(-2*x.f0/(deriv1*deriv1));
    return temp;
}
hyperdual sqrt(hyperdual x)
{
    return pow(x,0.5);
}
hyperdual fabs(hyperdual x)
{
    hyperdual temp;
    if (x < 0.0)
        temp = -x;
    else
        temp = x;
    return temp;
}
hyperdual max(hyperdual x1, hyperdual x2)
{
    hyperdual temp;
    if (x1>x2)
        temp = x1;
    else
        temp = x2;
    return temp;
}
hyperdual max(hyperdual x1, double x2)
{
    hyperdual temp;
    if (x1>x2)
        temp = x1;
    else
        temp = x2;
    return temp;
}
hyperdual max(double x1, hyperdual x2)
{
    hyperdual temp;
    if (x1>x2)
        temp = x1;
    else
        temp = x2;
    return temp;
}
hyperdual min(hyperdual x1, hyperdual x2)
{
    hyperdual temp;
    if (x1<x2)
        temp = x1;
    else
        temp = x2;
    return temp;
}
hyperdual min(hyperdual x1, double x2)
{
    hyperdual temp;
    if (x1<x2)
        temp = x1;
    else
        temp = x2;
    return temp;
}
hyperdual min(double x1, hyperdual x2)
{
    hyperdual temp;
    if (x1<x2)
        temp = x1;
    else
        temp = x2;
    return temp;
}

bool operator> (hyperdual lhs, hyperdual rhs)
{
    return (lhs.f0 > rhs.f0);
}
bool operator> (double lhs, hyperdual rhs)
{
    return (lhs > rhs.f0);
}
bool operator> (hyperdual lhs, double rhs)
{
    return (lhs.f0 > rhs);
}
bool operator>= (hyperdual lhs, hyperdual rhs)
{
    return (lhs.f0 >= rhs.f0);
}
bool operator>= (double lhs, hyperdual rhs)
{
    return (lhs >= rhs.f0);
}
bool operator>= (hyperdual lhs, double rhs)
{
    return (lhs.f0 >= rhs);
}
bool operator< (hyperdual lhs, hyperdual rhs)
{
    return (lhs.f0 < rhs.f0);
}
bool operator< (double lhs, hyperdual rhs)
{
    return (lhs < rhs.f0);
}
bool operator< (hyperdual lhs, double rhs)
{
    return (lhs.f0 < rhs);
}
bool operator<= (hyperdual lhs, hyperdual rhs)
{
    return (lhs.f0 <= rhs.f0);
}
bool operator<= (double lhs, hyperdual rhs)
{
    return (lhs <= rhs.f0);
}
bool operator<= (hyperdual lhs, double rhs)
{
    return (lhs.f0 <= rhs);
}
bool operator== (hyperdual lhs, hyperdual rhs)
{
    return (lhs.f0 == rhs.f0);
}
bool operator== (double lhs, hyperdual rhs)
{
    return (lhs == rhs.f0);
}
bool operator== (hyperdual lhs, double rhs)
{
    return (lhs.f0 == rhs);
}
bool operator!= (hyperdual lhs, hyperdual rhs)
{
    return (lhs.f0 != rhs.f0);
}
bool operator!= (double lhs, hyperdual rhs)
{
    return (lhs != rhs.f0);
}
bool operator!= (hyperdual lhs, double rhs)
{
    return (lhs.f0 != rhs);
}

#endif
