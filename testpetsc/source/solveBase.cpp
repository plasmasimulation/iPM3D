#include "solveBase.hpp"

// SolveTime
SolveTime::SolveTime() {
    MPI_Comm_size(PETSC_COMM_WORLD, &(this->size));
    MPI_Comm_rank(PETSC_COMM_WORLD, &(this->rank));

    // 计时 start
    if (0 == this->rank) {
        this->startTime = clock();
    }
}

SolveTime::~SolveTime() {
    // 计时 start
    if (!(this->rank)) {
        this->endTime = clock();
        double consumTime = (double)(this->endTime - this->startTime) / CLOCKS_PER_SEC * 1000;
        cout << "Time consuming: " << consumTime << " ms" << endl;
    }
}

// SolveBase
SolveBase::SolveBase() {
    // MPI
    MPI_Comm_size(PETSC_COMM_WORLD, &(this->size));
    MPI_Comm_rank(PETSC_COMM_WORLD, &(this->rank));

    // ksp
    KSPCreate(PETSC_COMM_WORLD, &(this->ksp));
    KSPSetType(this->ksp, KSPGMRES);
}

// 设置系数矩阵
void SolveBase::_setA(Mat A) {
    KSPSetOperators(this->ksp, A, A);
}

// 求解一次, A不变
void SolveBase::_run(Vec b) {    
    // 创建x
    if (nullptr == this->x) {
        VecDuplicate(b, &(this->x));
    }

    // 求解
    KSPSolve(this->ksp, b, this->x);
}

// 求解一次, A改变
void SolveBase::_run(Mat A, Vec b) {
    this->_setA(A);
    this->_run(b);
}

// 打印迭代次数与误差
void SolveBase::printInfo(Mat A, Vec b) {
    if (nullptr == this->r) {
        VecDuplicate(b, &(this->r));
    }

    // print
    PetscReal norms = -1;
    PetscInt its = -1;
    MatMult(A, this->x, this->r);
    VecAXPY(this->r, -1, b);
    VecNorm(this->r, NORM_2, &norms);
    KSPGetIterationNumber(ksp, &its);

    PetscPrintf(PETSC_COMM_WORLD, "Iterations: %D\n", its);
    PetscPrintf(PETSC_COMM_WORLD, "Norm      : %f\n", norms);
}

SolveBase::~SolveBase() {
    if (nullptr != this->ksp) {
        KSPDestroy(&(this->ksp));
        this->ksp = nullptr;
    }

    if (nullptr != this->x) {
        VecDestroy(&(this->x));
        this->x = nullptr;
    }

    if (nullptr != this->r) {
        VecDestroy(&(this->r));
        this->r = nullptr;
    }
}