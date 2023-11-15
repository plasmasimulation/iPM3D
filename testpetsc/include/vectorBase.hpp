#ifndef VECTOR_BASE_H
#define VECTOR_BASE_H

#include <iostream>
#include <string>
#include <math.h>
#include <sstream>
#include <vector>
#include <fstream>

#include <petscksp.h>
#include <petscmat.h>
#include <petscvec.h>
#include <petscdmda.h>

#define EPS (1e-15)
#define CHECKDMDA(da) if (nullptr == da) {cout << "Error: da is nullptr!" << endl; exit(-1);}

using namespace std;

// init petsc
extern void initPetsc();
extern void finalizePetsc();

// 形状数据管理
class DataManage {
public:
    DM                      da = nullptr;

    PetscInt                Mx = 33;                        // x/r 全局尺寸
    PetscInt                My = 33;                        // y/theta
    PetscInt                Mz = 1;                         // z

    PetscInt                npx = 1;                        // 本地进程数
    PetscInt                npy = 1;
    PetscInt                npz = 1;

    PetscScalar             x0 = 0.0;                       // 坐标值
    PetscScalar             y0 = 0.0;
    PetscScalar             z0 = 0.0;

    PetscScalar             dx = 1.0;                       // 空间步长
    PetscScalar             dy = 1.0;
    PetscScalar             dz = 1.0;

    PetscInt                xstart = 0;                     // 各维度索引范围
    PetscInt                ystart = 0;
    PetscInt                zstart = 0;

    PetscInt                xend = 0;
    PetscInt                yend = 0;
    PetscInt                zend = 0;

    PetscMPIInt             rank = 0;                       // mpi
    PetscMPIInt             size = 0;

    PetscInt                sharedWidth = 1;                // 进程共享数据宽度,在4次及4次以上需要增加

    PetscInt                cod = 0;                        // 0=直角坐标系  1=柱坐标系

    DataManage() = default;

    DataManage(const DataManage&) = delete;

    DataManage& operator=(const DataManage&) = delete;

    bool init(const PetscInt lz[] = NULL, const PetscInt ly[] = NULL, const PetscInt lx[] = NULL);

    bool isInit();

    void showInFile();

    ~DataManage();

private:

    bool                    isInitVar = false;

};

// 标量
class Scalar {
public:
    Vec             data = nullptr;
    DM              da = nullptr;
    DataManage      *dmgPtr = nullptr;
    PetscInt        offset = 0;

    // 构造与析构
    Scalar() = default;

    Scalar(DataManage & dmg);

    Scalar(const Scalar & src);

    Scalar(Scalar && src);

    Scalar& operator=(const Scalar & src);

    Scalar &operator=(Scalar && src);

    void init(DataManage & dmg);

    bool isInit();

    ~Scalar();

    // +-*/操作
    void add(const Scalar & scalar);

    void add(const PetscScalar & value);

    void sub(const Scalar & scalar);

    void sub(const PetscScalar & value);

    void mult(const Scalar & scalar);

    void mult(const PetscScalar & value);

    void div(const Scalar & scalar);

    void div(const PetscScalar & value);

    void set(const PetscScalar & value);

    void setInternal(const PetscScalar & value);

    void setBoundary(const PetscScalar & value);

    void copy(Vec data);

    void setVec(Vec data);

    void showInFile();

private:
    bool isInitVar = false;
};

// 矢量
class Vector {
public:
    vector<Scalar>  data;
    DM              da = nullptr;
    DataManage      *dmgPtr = nullptr;

    // 构造与析构
    Vector() = default;

    Vector(PetscInt dim, DataManage &dmg);

    Vector(const Vector & src);

    Vector(Vector && src);

    Vector& operator=(const Vector & src);

    Vector& operator=(Vector && src);

    void init(PetscInt dim, DataManage &dmg);

    bool isInit();

    // +-*/操作
    void add(const Scalar &scalar);
    
    void add(const Vector &vector);

    void add(const PetscScalar & value);

    void sub(const Scalar & scalar);

    void sub(const Vector & vector);

    void sub(const PetscScalar & value);

    void mult(const Scalar & scalar);

    void mult(const Vector & vector);

    void mult(const PetscScalar & value);

    void div(const Scalar & scalar);

    void div(const Vector & vector);

    void div(const PetscScalar & value);

    void set(const PetscScalar & value);

    void setInternal(const PetscScalar & value);

    void showInFile();

private:
    bool isInitVar = false;
};

// 差分
class BaseDif3D {
public:
    PetscInt lenx = 0;
    vector<vector<PetscInt>> indexLeftPx;
    vector<vector<PetscInt>> indexRightPx;

    PetscInt leny = 0;
    vector<vector<PetscInt>> indexLeftPy;
    vector<vector<PetscInt>> indexRightPy;

    PetscInt lenz = 0;
    vector<vector<PetscInt>> indexLeftPz;
    vector<vector<PetscInt>> indexRightPz;

    PetscInt len = 0;
    vector<vector<PetscInt>> indexAll;

    BaseDif3D() = default;

    BaseDif3D(const BaseDif3D &) = delete;
    BaseDif3D &operator=(const BaseDif3D &) = delete;

    void init(PetscInt Mx, PetscInt My, PetscInt Mz);
    inline PetscScalar px(const PetscScalar ***dataPtr, PetscInt i, PetscInt j, PetscInt k, PetscScalar dx);
    inline PetscScalar py(const PetscScalar ***dataPtr, PetscInt i, PetscInt j, PetscInt k, PetscScalar dy);
    inline PetscScalar pz(const PetscScalar ***dataPtr, PetscInt i, PetscInt j, PetscInt k, PetscScalar dz);
    inline PetscScalar sum(const PetscScalar ***dataPtr, PetscInt i, PetscInt j, PetscInt k);
};

// 矩阵
enum MatrixType {
    NormMatrixType = 0,
    MatrixFreeType,
};

class Matrix {
public:
    Mat             data = nullptr;
    DM              da = nullptr;
    DataManage      *dmgPtr = nullptr;
    MatrixType      type = NormMatrixType;
    PetscInt        matStoreNum = 5;                // 矩阵预分配内存, 用于计算提速

    void            (*myMult)(void) = nullptr;      // 自定义矩阵乘法函数指针
    void            *usrPtr = nullptr;

    // 构造析构
    Matrix() = default;

    Matrix(const Matrix &) = delete;
    Matrix(Matrix &&) = delete;
    Matrix& operator=(const Matrix&) = delete;
    
    Matrix(DataManage & dmg, MatrixType type, PetscInt matStoreNum, void *mymult=nullptr, void *usrPtr=nullptr);
    void init(DataManage & dmg, MatrixType type, PetscInt matStoreNum, void *mymult=nullptr, void *usrPtr=nullptr);
    bool isInit();

    ~Matrix();

private:
    bool isInitVar = false;

};

// 3阶张量
class Tensor3d {
public:
    Scalar xx;
    Scalar xy;
    Scalar xz;

    Scalar yx;
    Scalar yy;
    Scalar yz;

    Scalar zx;
    Scalar zy;
    Scalar zz;

    Tensor3d() = default;

    Tensor3d(const Tensor3d &) = delete;
    Tensor3d(const Tensor3d &&) = delete;
    Tensor3d& operator=(const Tensor3d&) = delete;

    Tensor3d(DataManage & dmg)
    {
       this->xx.init(dmg);
       this->xy.init(dmg);
       this->xz.init(dmg);
       this->yx.init(dmg);
       this->yy.init(dmg);
       this->yz.init(dmg);
       this->zx.init(dmg);
       this->zy.init(dmg);
       this->zz.init(dmg);
    }
};

// DataManage 重载输出运算符
extern ostream & operator<<(ostream & out, const DataManage & dataManage);

// Scalar
extern Scalar add(const Scalar &left, const Scalar &right);

extern Scalar add(const Scalar &left, const PetscScalar &right);

extern Scalar add(const PetscScalar &left, const Scalar &right);

extern Scalar sub(const Scalar &left, const Scalar &right);

extern Scalar sub(const Scalar &left, const PetscScalar &right);

extern Scalar sub(const PetscScalar &left, const Scalar &right);

extern Scalar mult(const Scalar &left, const Scalar &right);

extern Scalar mult(const Scalar &left, const PetscScalar &right);

extern Scalar mult(const PetscScalar &left, const Scalar &right);

extern Scalar div(const Scalar &left, const Scalar &right);

extern Scalar div(const Scalar &left, const PetscScalar &right);

extern Scalar div(const PetscScalar &left, const Scalar &right);

extern ostream & operator<<(ostream & out, const Scalar & scalar);

// Vector
extern Vector add(const Vector &left, const Vector &right);

extern Vector add(const Vector &left, const PetscScalar &right);

extern Vector add(const PetscScalar &left, const Vector &right);

extern Vector sub(const Vector &left, const Vector &right);

extern Vector sub(const Vector &left, const PetscScalar &right);

extern Vector sub(const PetscScalar &left, const Vector &right);

extern Vector mult(const Vector &left, const Vector &right);

extern Vector mult(const Vector &left, const PetscScalar &right);

extern Vector mult(const PetscScalar &left, const Vector &right);

extern Vector div(const Vector &left, const Vector &right);

extern Vector div(const Vector &left, const PetscScalar &right);

extern Vector div(const PetscScalar &left, const Vector &right);

extern ostream & operator<<(ostream & out, const Vector & vector);

// 初始化算子
extern void initOperator(PetscInt Mx, PetscInt My, PetscInt Mz);

extern ostream & operator<<(ostream & out, const BaseDif3D & baseDif3D);

extern void printBaseOperator();

// 给定标量场，计算梯度
extern Vector getGrad(const Scalar &scalarInput);

// 给定矢量场，计算散度
extern Scalar getDiv(const Vector &vectorInput);

// 给定矢量场，计算旋度
extern Vector getCurl(const Vector &vectorInput);

// 给定场的单个分量，计算拉普拉斯算子
extern Scalar getLap(const Scalar &scalarInput);

// 计算矢量的拉普拉斯算子
extern Vector getLap(const Vector &vectorInput);

// 3阶张量与场的3个分量的点乘   Node
extern Vector getTensorDot(const Tensor3d & tensorInput,  const Vector & vectorInput);

// python 画二维向量的3d图
void pyplot(string scriptName, string outFileName, Vec src, PetscInt xsize, PetscInt ysize);

#endif // VECTOR_BASE_H