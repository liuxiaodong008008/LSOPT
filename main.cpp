#include "matrixlib/matrix.h"
#include "matrixlib/matrix_operation.h"
#include "hyperdual.h"
#include "optimizer.h"
#include "time.h"
#include <iostream>
#include <functional>
#include <iomanip>
using namespace std;

template <typename T>
struct MyFunctor {
    T operator()(const Matrix<T>& m) {
        const T &a = m(0);
        const T &b = m(1);
        const T &xeps = m(idx+2);
        return a*(x+xeps)+b - y;
    }

    double x,y;
    int idx;

    MyFunctor(int idx, double x, double y) : x(x), y(y), idx(idx) {}
};

template <typename T>
struct MyFunctor2 {
    T operator()(const Matrix<T>& m) {
        return m(idx+2);
    }

    int idx;
    MyFunctor2(int idx) : idx(idx) {}
};

int main() {
    Matrix<double> XY(30,2,0.0);

    XY.row( 0)<<1.1,39343.00;
    XY.row( 1)<<1.3,46205.00;
    XY.row( 2)<<1.5,37731.00;
    XY.row( 3)<<2.0,43525.00;
    XY.row( 4)<<2.2,39891.00;
    XY.row( 5)<<2.9,56642.00;
    XY.row( 6)<<3.0,60150.00;
    XY.row( 7)<<3.2,54445.00;
    XY.row( 8)<<3.2,64445.00;
    XY.row( 9)<<3.7,57189.00;
    XY.row(10)<<3.9,63218.00;
    XY.row(11)<<4.0,55794.00;
    XY.row(12)<<4.0,56957.00;
    XY.row(13)<<4.1,57081.00;
    XY.row(14)<<4.5,61111.00;
    XY.row(15)<<4.9,67938.00;
    XY.row(16)<<5.1,66029.00;
    XY.row(17)<<5.3,83088.00;
    XY.row(18)<<5.9,81363.00;
    XY.row(19)<<6.0,93940.00;
    XY.row(20)<<6.8,91738.00;
    XY.row(21)<<7.1,98273.00;
    XY.row(22)<<7.9,101302.00;
    XY.row(23)<<8.2,113812.00;
    XY.row(24)<<8.7,109431.00;
    XY.row(25)<<9.0,105582.00;
    XY.row(26)<<9.5,116969.00;
    XY.row(27)<<9.6,112635.00;
    XY.row(28)<<10.3,122391.00;
    XY.row(29)<<10.5,121872.00;

    cout<<XY<<endl;

    {
        Cost cost(32);
        for(int r=0;r<XY.rows();r++) {
            auto cr = XY.row(r);
            cost.push_back(new_residual<MyFunctor>(cost.vars(), r, cr(0),cr(1)));
            cost.push_back(new_residual<MyFunctor2>(cost.vars(),r));
        }

        BacktrackingStepSolver ss(1.0, 0.5, -1e-15);
        GaussNewtonMethod opt(cost, ss, 1e-6);

        time_t start_time, end_time;
        start_time = clock();
        opt.minimize();
        end_time = clock();
        cout<<"Gauss-Newton Method"<<endl;
        cout << cost.vars()(Range({0,0},{1,2})).t() << endl;
        cout<<"time used: "<<setprecision(9)<<(end_time-start_time+0.0)/CLOCKS_PER_SEC<<endl;
    }

    cout<<"==============="<<endl;

    {
        Cost cost(32);
        for(int r=0;r<XY.rows();r++) {
            auto cr = XY.row(r);
            cost.push_back(new_residual<MyFunctor>(cost.vars(), r, cr(0),cr(1)));
            cost.push_back(new_residual<MyFunctor2>(cost.vars(),r));
        }

        LevenbergMarquardtMethod opt(cost,1e-6,1e4,1e6,0.15);
        time_t start_time, end_time;
        start_time = clock();
        opt.minimize();
        end_time = clock();
        cout<<"Levenberg-Marquardt Method"<<endl;
        cout << cost.vars()(Range({0,0},{1,2})).t() << endl;
        cout<<"time used: "<<setprecision(9)<<(end_time-start_time+0.0)/CLOCKS_PER_SEC<<endl;
    }

    return 0;
}
