//
// Created by liuxi on 2019/1/7.
//

#ifndef LSOPT_PROBLEM_H
#define LSOPT_PROBLEM_H

#include <vector>
#include <functional>
#include "matrixlib/matrix.h"
#include "matrixlib/matrix_operation.h"
#include "hyperdual.h"
#include <cassert>

struct Residual {
    virtual double eval() { return eval(vars()); }
    virtual Matrix<double> grad() { return grad(vars()); }

    virtual double eval(const IMatrix<double>& v) = 0;
    virtual Matrix<double> grad(const IMatrix<double>& v) = 0;
    virtual Matrix<double>& vars() = 0;

    virtual ~Residual() {}
};

template <template <typename > class Func>
struct ResidualFunctor:Residual {
    Func<double>* eval_value;
    Func<hyperdual>* eval_hyperdual;
    Matrix<double>& variables;
    ResidualFunctor(Matrix<double> &variables)
            : variables(variables), eval_value(nullptr), eval_hyperdual(nullptr) {}

    double eval(const IMatrix<double> &v) override {
        return (*eval_value)(v);
    }

    Matrix<double> grad(const IMatrix<double> &v) override {
        int n = v.count();
        Matrix<double> g(n,1,0);
        Matrix<hyperdual> m(n,1,hyperdual());

        for(int i=0;i<n;i++)
            m(i).f0 = v(i);

        for(int i=0;i<n;) {
            if (i+1<n) {
                m(i).f1 = 1.0f;
                m(i+1).f2 = 1.0f;
                hyperdual h = (*eval_hyperdual)(m);
                g(i) = h.eps1();
                g(i+1) = h.eps2();
                m(i).f1 = 0.0f;
                m(i+1).f2 = 0.0f;
                i+=2;
            } else {
                m(i).f1 = 1.0f;
                hyperdual h = (*eval_hyperdual)(m);
                g(i) = h.eps1();
                m(i).f1 = 0.0f;
                i+=1;
            }
        }
        return g;
    }

    Matrix<double> &vars() override {
        return this->variables;
    }

    virtual ~ResidualFunctor() {
        if(eval_value!= nullptr) {
            delete eval_value;
            eval_value = nullptr;
        }

        if(eval_hyperdual != nullptr) {
            delete eval_hyperdual;
            eval_hyperdual = nullptr;
        }
    }
};


template <template <typename > class Func, typename... Args>
Residual* new_residual(Matrix<double> &variables, Args... args){
    auto n = new ResidualFunctor<Func>(variables);
    n->eval_value = new Func<double>(args...);
    n->eval_hyperdual = new Func<hyperdual>(args...);
    return n;
}


struct Cost :std::vector<Residual*> {
    Matrix<double> variables;


    Matrix<double> evalres(const Matrix<double>& v) {
        Matrix<double> r(this->size(),1,0);
        for(int i=0;i<this->size();i++)
            r(i) = this->at(i)->eval(v);
        return r;
    }

    double eval(const Matrix<double>& v) {
        auto e = evalres(v);
        return 0.5*sum(e*e,REDUCE_DIRECTION::ALL)(0);
    }

    Matrix<double> jacobian(const Matrix<double>& v) {
        int cn = v.count();
        int rn = this->size();

        Matrix<double> j(rn,cn,0);
        for(int r=0;r<rn;r++)
            j.row(r) = (*this)[r]->grad(v).t();
        return j;
    }

    Matrix<double> grad(const Matrix<double>& v) {
        return matmul(jacobian(v).t(),evalres(v));
    }


    Matrix<double> evalres() {
        return evalres(this->vars());
    }

    double eval() {
        return eval(this->vars());
    }

    Matrix<double> jacobian() {
        return jacobian(this->vars());
    }

    Matrix<double> grad() {
        return grad(this->vars());
    }

    Matrix<double> &vars() {
        return this->variables;
    }

    Cost(int n) : variables(n,1,0) {}

    Cost(const Matrix<double>& init) : variables(init.count(),1,0) {
        int n = init.count();
        for(int i=0;i<n;i++) variables(i) = init(i);
    }

    Cost(const std::vector<double>& init) : variables(init.size(),1,0) {
        int n = init.size();
        for(int i=0;i<n;i++) variables(i) = init[i];
    }

    Cost(const std::initializer_list<double>& init) : variables(init.size(),1,0) {
        int n = init.size();
        auto it = init.begin();
        for(int i=0;i<n;i++,it++) variables(i) = *it;
    }

    virtual ~Cost() {
        for(auto p:*this) {
            delete p;
            p = nullptr;
        }
    }
};

struct Optimizer{
    virtual void minimize() = 0;
};

struct StepSolver{
    virtual double get_step(Cost & cost, const Matrix<double>& direct, const Matrix<double>& grad) = 0;
};

struct BacktrackingStepSolver:StepSolver {
    double get_step(Cost &cost, const Matrix<double> &direct, const Matrix<double>& grad) override {
        double alpha = this->alpha;
        double rho = this->rho;
        double c = this->c;

        Matrix<double> vars = cost.vars();
        double descent_factor = matmul(grad.t(),direct)(0);

        auto cur_cost = eval(cost,vars);
        int i=0;
        while(true){
            i++;
            auto new_cost = eval(cost,vars+alpha*direct);
            if(new_cost<=cur_cost+c*alpha*descent_factor) break;
            else alpha = alpha*rho;
        }
        return alpha;
    }

    double eval(Cost &cost, const Matrix<double>& v) {
        auto e = cost.evalres();
        return sum(e*e,REDUCE_DIRECTION::ALL)(0);
    }

    double alpha;
    double rho;
    double c;
    BacktrackingStepSolver(double alpha, double rho, double c)
            : alpha(alpha), rho(rho), c(c){}
};

struct GaussNewtonMethod : Optimizer {
    StepSolver& stepSolver;
    Cost& cost;
    double epsilon;

    GaussNewtonMethod(Cost &cost, StepSolver &stepSolver, double epsilon)
    : stepSolver(stepSolver), cost(cost), epsilon(epsilon) {}

    void minimize() override {
        Matrix<double> p(cost.vars().count(),1,0);
        Matrix<double> grad(cost.vars().count(),1,0);
        Matrix<double> jac(cost.size(),cost.vars().count(),0);

        while (true) {
            grad = cost.grad();
            if(magnitude(grad)<=epsilon) break;
            jac = cost.jacobian();
            p = -matmul(inverse(matmul(jac.t(),jac)),grad);
            double alpha = stepSolver.get_step(cost,p,grad);
            cost.vars()+=alpha*p;
        }
    }
};

struct LevenbergMarquardtMethod : Optimizer {
    Cost& cost;
    double epsilon;
    double init_range;
    double max_range;
    double eta;

    LevenbergMarquardtMethod(Cost &cost, double epsilon,
            double init_range, double max_range, double eta)
            : cost(cost),epsilon(epsilon),init_range(init_range),
              max_range(max_range),eta(eta) {}

    void minimize() override {
        Matrix<double> p(cost.vars().count(),1,0);
        Matrix<double> grad(cost.vars().count(),1,0);
        Matrix<double> jac(cost.size(),cost.vars().count(),0);

        double range = this->init_range;
        double cur_cost, next_cost, pm;

        cur_cost = cost.eval();

        while (true) {
            grad = cost.grad();
            double gradm = magnitude(grad);
            if(gradm<=epsilon) break;

            jac = cost.jacobian();
            auto B = matmul(jac.t(),jac);
            auto pb = -matmul(inverse(B),grad);
            pm = magnitude(pb);
            if(pm<=range) p = pb;
            else {
                auto pu = -matmul(grad.t(),grad)(0)*grad/(matmul(matmul(grad.t(),B),grad))(0);
                auto pm = magnitude(pu);
                if(pm>=range) p = range/pm*pu;
                else {
                    auto pbu = pb-pu;
                    double a = matmul(pbu.t(),pbu)(0);
                    double b = matmul(pu.t(),pbu)(0);
                    double c = matmul(pu.t(),pu)(0) - range*range;
                    double t = (-b+sqrt(b*b-4*a*c))/(2*a);
                    p = pu+t*pbu;
                    pm = magnitude(p);
                }
            }
            next_cost = cost.eval(cost.vars()+p);
            double rho = -(cur_cost - next_cost)/
                    (matmul(grad.t(),p)+0.5*matmul(matmul(p.t(),B),p))(0);

            if(rho<0.25)
                range *= 0.25;
            else if(rho>0.75 && abs(pm/range-1)<1e-2 )
                range = min(2*range,this->max_range);

            if(rho>this->eta) {
                cost.vars()+=p;
                cur_cost = next_cost;
            }
        }
    }
};


#endif //LSOPT_PROBLEM_H
