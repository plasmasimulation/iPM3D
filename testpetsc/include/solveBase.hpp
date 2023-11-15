#ifndef SOLVE_BASE_H
#define SOLVE_BASE_H

#include <iostream>
#include <string>
#include <ctime>

#include <petscksp.h>
#include <petscmat.h>
#include <petscvec.h>
#include <petscdmda.h>

#include "vectorBase.hpp"

using namespace std;

class SolveTime {
private:
    clock_t                 startTime, endTime;     // 计时器
    PetscMPIInt             rank = 0;               // mpi
    PetscMPIInt             size = 0;

public:
    SolveTime();

    SolveTime(const SolveTime &) = delete;
    SolveTime(SolveTime &&) = delete;
    SolveTime &operator=(const SolveTime &) = delete;

    ~SolveTime();
};

class SolveBase {
public:
    Vec                     x = nullptr;
    Vec                     r = nullptr;            // r=Ax-b
    KSP                     ksp;                    // 求解器

    PetscMPIInt             rank = 0;               // mpi
    PetscMPIInt             size = 0;

    // 构造 析构
    SolveBase();

    SolveBase(const SolveBase & src) = delete;
    SolveBase(SolveBase &&) = delete;
    SolveBase& operator=(const SolveBase&) = delete;

    ~SolveBase();
    
    // 求解
    void _setA(Mat A);
    void _run(Mat A, Vec b);
    void _run(Vec b);
    void printInfo(Mat A, Vec b);
};

class ScalarSolver {
public:
    DataManage              dm;
    SolveBase               solve;
    Scalar                  b;
    Matrix                  A;

    bool                    isInitVar = false;

    // 默认构造
    ScalarSolver() = default;

    ScalarSolver(const ScalarSolver &) = delete;
    ScalarSolver(ScalarSolver &&) = delete;
    ScalarSolver &operator=(const ScalarSolver &) = delete;

    // 初始化 dm: 网格形状  A: 系数矩阵  b: 源项
    virtual void init() = 0;

    virtual void setb() = 0;

    virtual void setA() {
        if (this->isInitVar)    this->solve._setA(A.data);
        else                    cout << "Eroor: Petsc solver don't init, please call init()" \
                                     << endl;
    }

    virtual void run() {
        if (this->isInitVar)    this->solve._run(this->A.data, this->b.data);
        else                    cout << "Eroor: Petsc solver don't init, please call init()" \
                                     << endl;
    }
};

#endif // SOLVE_BASE_H