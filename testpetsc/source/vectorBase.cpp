#include "vectorBase.hpp"

static char help[] = "Solves 3D field.\n\n";

static bool petscIsInit = false;
static bool petscIsFinalize = false;

#define DEBUG_SCALAR_U
#define DEBUG_VECTOR_U

BaseDif3D baseDif;

// init petsc
void initPetsc() {
    if (false == petscIsInit) {
        PetscInitialize(NULL, NULL, (char *)0, help);
        petscIsInit = true;
    }
}

void finalizePetsc() {
    if (false == petscIsFinalize) {
        PetscFinalize();
        petscIsFinalize = true;
    }
}

// DataManage
bool DataManage::init(const PetscInt lz[], const PetscInt ly[], const PetscInt lx[]) {
    if (nullptr != this->da) {
        DMDestroy(&(this->da));
        this->da = nullptr;
    }

    // MPI
    MPI_Comm_size(PETSC_COMM_WORLD, &(this->size));
    MPI_Comm_rank(PETSC_COMM_WORLD, &(this->rank));

    auto stencilType = DMDA_STENCIL_STAR;
    PetscInt count_dim = 0;
    if (this->Mx > 1) count_dim++;
    if (this->My > 1) count_dim++;
    if (this->Mz > 1) count_dim++;
    if (count_dim > 1) stencilType = DMDA_STENCIL_BOX;

    // 创建dmda数据管理器
    if (this->npx * this->npy * this->npz != this->size) {
        if (!this->rank) {
            std::cout << "Error for creating dmda in PETSc." << std::endl;
            std::cout << "The number of image is : " << this->size << std::endl;
            std::cout << "The number of image in each direction is: " 
                      << this->npx << " " << this->npy << " " << this->npz << std::endl;
        }

        DMDACreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                     stencilType, this->Mz, this->My, this->Mx, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE,
                     1, this->sharedWidth, NULL, NULL, NULL, &(this->da));
    }
    else {
        DMDACreate3d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                     stencilType, this->Mz, this->My, this->Mx, this->npz, this->npy, this->npx,
                     1, this->sharedWidth, lz, ly, lx, &(this->da));
    }

    DMSetFromOptions(this->da);
    DMSetUp(this->da);

    // 设置索引范围
    DMDAGetCorners(this->da, &(this->zstart), &(this->ystart), &(this->xstart),
                   &(this->zend), &(this->yend), &(this->xend));

    // 状态
    this->isInitVar = true;

    initOperator(this->Mx, this->My, this->Mz);
}

DataManage::~DataManage() {
    if (nullptr != this->da) {
        DMDestroy(&(this->da));
        this->da = nullptr;
    }
}

bool DataManage::isInit() {
    return this->isInitVar;
}

ostream & operator<<(ostream & out, const DataManage & dataManage) {
    out << "DataManage object shape : " << dataManage.Mx
        << " " << dataManage.My << " " << dataManage.Mz << endl;

    out << "Local process number : " << dataManage.npx
        << " " << dataManage.npy << " " << dataManage.npz << endl;

    out << "Space start value : " << dataManage.x0
        << " " << dataManage.y0 << " " << dataManage.z0 << endl;

    out << "Space step value : " << dataManage.dx
        << " " << dataManage.dy << " " << dataManage.dz << endl;

    out << "Local index start : " << dataManage.xstart
        << " " << dataManage.ystart << " " << dataManage.zstart << endl;

    out << "Local index width : " << dataManage.xend
        << " " << dataManage.yend << " " << dataManage.zend << endl;

    out << "Shared width : " << dataManage.sharedWidth << endl;

    if (dataManage.cod == 0) {
        out << "The cartesian coordinate system is used." << endl;
    }
    else if (dataManage.cod == 1) {
        out << "The column coordinate system is used." << endl;
    }

    out << endl;

    return out;
}

void DataManage::showInFile() {
    string fileName = "dm_shape_";
    fileName += to_string(this->rank);
    fileName += ".dat";

    ofstream out(fileName);

    out << *this;

    out.close();
}

// Scalar
Scalar::Scalar(DataManage & dmg) {
    this->init(dmg);
}

void Scalar::init(DataManage & dmg) {
    if (nullptr != this->data) {
        VecDestroy(&(this->data));
        this->data = nullptr;
    }

    this->da = dmg.da;
    this->dmgPtr = &dmg;

    CHECKDMDA(this->da);
    DMCreateGlobalVector(this->da, &(this->data));

    this->offset = 0;
    this->isInitVar = true;
    
    #ifdef DEBUG_SCALAR
    if (!this->dmgPtr->rank) cout << "Scalar: Construction from DataManage." << endl;
    #endif
}

Scalar::Scalar(const Scalar & src) {
    if (src.isInitVar) {
        if (nullptr != this->data)
        {
            VecDestroy(&(this->data));
            this->data = nullptr;
        }

        this->da = src.da;
        this->dmgPtr = src.dmgPtr;
        this->offset = src.offset;

        CHECKDMDA(this->da);
        DMCreateGlobalVector(this->da, &(this->data));
        VecCopy(src.data, this->data);

        this->isInitVar = true;

        #ifdef DEBUG_SCALAR
        if (!this->dmgPtr->rank) cout << "Scalar: Construction from other scalar." << endl;
        #endif
    }
}

Scalar::Scalar(Scalar && src) {
    if (src.isInitVar) {
        if (nullptr != this->data) {
            VecDestroy(&(this->data));
            this->data = nullptr;
        }

        this->da = src.da;
        this->dmgPtr = src.dmgPtr;
        this->offset = src.offset;

        CHECKDMDA(this->da);
        this->data = src.data;
        src.data = nullptr;

        this->isInitVar = true;
        
        #ifdef DEBUG_SCALAR
        if (!this->dmgPtr->rank) cout << "Scalar: Construction from other scalar by moving." << endl;
        #endif
    }
}

Scalar& Scalar::operator=(const Scalar & src) {
    if (src.isInitVar) {
        if (nullptr != this->data) {
            VecDestroy(&(this->data));
            this->data = nullptr;
        }

        this->da = src.da;
        this->dmgPtr = src.dmgPtr;
        this->offset = src.offset;

        CHECKDMDA(this->da);
        DMCreateGlobalVector(this->da, &(this->data));
        VecCopy(src.data, this->data);

        this->isInitVar = true;

        #ifdef DEBUG_SCALAR
        if (!this->dmgPtr->rank) cout << "Scalar: Construction from other scalar by assigning." << endl;
        #endif
    }
}

Scalar& Scalar::operator=(Scalar && src) {
    if (src.isInitVar) {
        if (nullptr != this->data) {
            VecDestroy(&(this->data));
            this->data = nullptr;
        }

        this->da = src.da;
        this->dmgPtr = src.dmgPtr;
        this->offset = src.offset;

        CHECKDMDA(this->da);
        this->data = src.data;
        src.data = nullptr;

        this->isInitVar = true;
        
        #ifdef DEBUG_SCALAR
        if (!this->dmgPtr->rank) cout << "Scalar: Construction from other scalar by moving." << endl;
        #endif
    }
}

bool Scalar::isInit() {
    return this->isInitVar;
}

Scalar::~Scalar() {
    if (nullptr != this->data)
    {
        VecDestroy(&(this->data));
        this->data = nullptr;
    }
}

void Scalar::add(const Scalar & scalar) {
    VecAXPY(this->data, 1, scalar.data);
}

void Scalar::add(const PetscScalar & value) {
    VecShift(this->data, value);
}

void Scalar::sub(const Scalar & scalar) {
    VecAXPY(this->data, -1, scalar.data);
}

void Scalar::sub(const PetscScalar & value) {
    VecShift(this->data, -1*value);
}

void Scalar::mult(const Scalar & scalar) {
    PetscInt                i, j, k;
    const PetscScalar       ***rightPtr;
    PetscScalar             ***outPtr;
    PetscInt                &xstart = this->dmgPtr->xstart;
    PetscInt                &ystart = this->dmgPtr->ystart;
    PetscInt                &zstart = this->dmgPtr->zstart;
    PetscInt                &xend = this->dmgPtr->xend;
    PetscInt                &yend = this->dmgPtr->yend;
    PetscInt                &zend = this->dmgPtr->zend;

    // 获取指针
    DMDAVecGetArrayRead(scalar.da, scalar.data, &rightPtr);
    DMDAVecGetArray(this->da, this->data, &outPtr);

    for (i = xstart; i < xstart + xend; i++) {
        for (j = ystart; j < ystart + yend; j++) {
            for (k = zstart; k < zstart + zend; k++) {
                outPtr[i][j][k] = outPtr[i][j][k] * rightPtr[i][j][k];
            }
        }
    }

    DMDAVecRestoreArrayRead(scalar.da, scalar.data, &rightPtr);
    DMDAVecRestoreArray(this->da, this->data, &outPtr);
}

void Scalar::mult(const PetscScalar & value) {
    VecScale(this->data, value);
}

void Scalar::div(const Scalar & scalar) {
    PetscInt                i, j, k;
    const PetscScalar       ***rightPtr;
    PetscScalar             ***outPtr;
    PetscInt                &xstart = this->dmgPtr->xstart;
    PetscInt                &ystart = this->dmgPtr->ystart;
    PetscInt                &zstart = this->dmgPtr->zstart;
    PetscInt                &xend = this->dmgPtr->xend;
    PetscInt                &yend = this->dmgPtr->yend;
    PetscInt                &zend = this->dmgPtr->zend;

    // 获取指针
    DMDAVecGetArrayRead(scalar.da, scalar.data, &rightPtr);
    DMDAVecGetArray(this->da, this->data, &outPtr);

    for (i = xstart; i < xstart + xend; i++) {
        for (j = ystart; j < ystart + yend; j++) {
            for (k = zstart; k < zstart + zend; k++) {
                if (abs(rightPtr[i][j][k]) > EPS)
                    outPtr[i][j][k] = outPtr[i][j][k] / rightPtr[i][j][k];
            }
        }
    }

    DMDAVecRestoreArrayRead(scalar.da, scalar.data, &rightPtr);
    DMDAVecRestoreArray(this->da, this->data, &outPtr);
}

void Scalar::div(const PetscScalar & value) {
    if (abs(value) > EPS)
        VecScale(this->data, 1/value);
}

void Scalar::set(const PetscScalar & value) {
    VecSet(this->data, value);
}

void Scalar::setInternal(const PetscScalar & value) {
    PetscInt                i, j, k;
    PetscScalar             ***xPtr;
    PetscInt                &Mx = this->dmgPtr->Mx;
    PetscInt                &My = this->dmgPtr->My;
    PetscInt                &Mz = this->dmgPtr->Mz;
    PetscInt                &xstart = this->dmgPtr->xstart;
    PetscInt                &ystart = this->dmgPtr->ystart;
    PetscInt                &zstart = this->dmgPtr->zstart;
    PetscInt                &xend = this->dmgPtr->xend;
    PetscInt                &yend = this->dmgPtr->yend;
    PetscInt                &zend = this->dmgPtr->zend;

    DMDAVecGetArray(this->da, this->data, &xPtr);

    for (i = xstart; i < xstart + xend; i++) {
        for (j = ystart; j < ystart + yend; j++) {
            for (k = zstart; k < zstart + zend; k++) {
                if ((Mx > 1 && (0 == i || Mx-1 == i))
                    || (My > 1 && (0 == j || My-1 == j))
                    || (Mz > 1 && (0 == k || Mz-1 == k)))
                    continue;
                else
                    xPtr[i][j][k] = value;
            }
        }
    }

    DMDAVecRestoreArray(this->da, this->data, &xPtr);
}

void Scalar::copy(Vec data) {
    VecCopy(data, this->data);
}

void Scalar::setVec(Vec data) {
    VecCopy(this->data, data);
}

ostream & operator<<(ostream & out, const Scalar & scalar) {
    PetscInt            i, j, k;
    const PetscScalar   ***xPtr;
    
    DMDAVecGetArrayRead(scalar.da, scalar.data, &xPtr);

    out << "Scalar: [" << scalar.dmgPtr->rank << "]" << endl;
    for (k = scalar.dmgPtr->zstart; k < scalar.dmgPtr->zstart + scalar.dmgPtr->zend; k++) {
        for (j = scalar.dmgPtr->ystart; j < scalar.dmgPtr->ystart + scalar.dmgPtr->yend; j++) {
            for (i = scalar.dmgPtr->xstart; i < scalar.dmgPtr->xstart + scalar.dmgPtr->xend; i++) {
                out << xPtr[i][j][k] << "\t";
            }
            out << endl;
        }
        out << endl;
    }
    out << endl;

    DMDAVecRestoreArrayRead(scalar.da, scalar.data, &xPtr);

    return out;
}

void Scalar::showInFile() {
    static PetscInt count = 0;
    string fileName = "scalar_data_";
    fileName += to_string(count++);
    fileName += "_";
    fileName += to_string(this->dmgPtr->rank);
    fileName += ".dat";

    ofstream out(fileName);

    out << *this;

    out.close();
}

// Vector
Vector::Vector(PetscInt dim, DataManage &dmg) {
    this->init(dim, dmg);
}

void Vector::init(PetscInt dim, DataManage & dmg) {
    if (dim > 0) {
        this->da = dmg.da;
        this->dmgPtr = &dmg;
        CHECKDMDA(this->da);

        this->data.clear();
        this->data = vector<Scalar>(dim);

        for (auto i = 0; i < dim; i++) {
            this->data[i].init(dmg);
        }

        this->isInitVar = true;

        #ifdef DEBUG_VECTOR
        if (!this->dmgPtr->rank) cout << "Vector: Construction from DataManage." << endl;
        #endif
    }
}

Vector::Vector(const Vector & src) {
    if (src.isInitVar) {
        this->da = src.da;
        this->dmgPtr = src.dmgPtr;
        CHECKDMDA(this->da);

        this->data.clear();
        this->data = vector<Scalar>(src.data.size());
        for (auto i = 0; i < src.data.size(); i++) {
            this->data[i] = src.data[i];
        }

        this->isInitVar = true;

        #ifdef DEBUG_VECTOR
        if (!this->dmgPtr->rank) cout << "Vector: Construction from other vector." << endl;
        #endif
    }
}

Vector::Vector(Vector && src) {
    if (src.isInitVar) {
        this->da = src.da;
        this->dmgPtr = src.dmgPtr;
        CHECKDMDA(this->da);

        this->data.clear();
        this->data = vector<Scalar>(src.data.size());
        for (auto i = 0; i < src.data.size(); i++) {
            this->data[i] = std::move(src.data[i]);
        }

        this->isInitVar = true;

        #ifdef DEBUG_VECTOR
        if (!this->dmgPtr->rank) cout << "Vector: Construction from other vector by moving." << endl;
        #endif
    }
}

Vector& Vector::operator=(const Vector & src) {
    if (src.isInitVar) {
        this->da = src.da;
        this->dmgPtr = src.dmgPtr;
        CHECKDMDA(this->da);

        this->data.clear();
        this->data = vector<Scalar>(src.data.size());
        for (auto i = 0; i < src.data.size(); i++) {
            this->data[i] = src.data[i];
        }

        this->isInitVar = true;

        #ifdef DEBUG_VECTOR
        if (!this->dmgPtr->rank) cout << "Vector: Construction from other vector by assigning." << endl;
        #endif
    }
}

Vector& Vector::operator=(Vector && src) {
    if (src.isInitVar) {
        this->da = src.da;
        this->dmgPtr = src.dmgPtr;
        CHECKDMDA(this->da);

        this->data.clear();
        this->data = vector<Scalar>(src.data.size());
        for (auto i = 0; i < src.data.size(); i++) {
            this->data[i] = std::move(src.data[i]);
        }

        this->isInitVar = true;

        #ifdef DEBUG_VECTOR
        if (!this->dmgPtr->rank) cout << "Vector: Construction from other vector by moving." << endl;
        #endif
    }
}

bool Vector::isInit() {
    return this->isInitVar;
}

void Vector::add(const Scalar & scalar) {
    for (auto i = 0; i < this->data.size(); i++) {
        this->data[i].add(scalar);
    }
}

void Vector::add(const Vector & vector) {
    for (auto i = 0; i < this->data.size(); i++) {
        this->data[i].add(vector.data[i]);
    }
}

void Vector::add(const PetscScalar & value) {
    for (auto i = 0; i < this->data.size(); i++) {
        this->data[i].add(value);
    }
}

void Vector::sub(const Scalar & scalar) {
    for (auto i = 0; i < this->data.size(); i++) {
        this->data[i].sub(scalar);
    }
}

void Vector::sub(const Vector & vector) {
    for (auto i = 0; i < this->data.size(); i++) {
        this->data[i].sub(vector.data[i]);
    }
}

void Vector::sub(const PetscScalar & value)
{
    for (auto i = 0; i < this->data.size(); i++) {
        this->data[i].sub(value);
    }
}

void Vector::mult(const Scalar & scalar) {
    for (auto i = 0; i < this->data.size(); i++) {
        this->data[i].mult(scalar);
    }
}

void Vector::mult(const Vector & vector) {
    for (auto i = 0; i < this->data.size(); i++) {
        this->data[i].mult(vector.data[i]);
    }
}

void Vector::mult(const PetscScalar & value) {
    for (auto i = 0; i < this->data.size(); i++) {
        this->data[i].mult(value);
    }
}

void Vector::div(const Scalar & scalar) {
    for (auto i = 0; i < this->data.size(); i++) {
        this->data[i].div(scalar);
    }
}

void Vector::div(const Vector & vector) {
    for (auto i = 0; i < this->data.size(); i++) {
        this->data[i].div(vector.data[i]);
    }
}

void Vector::div(const PetscScalar & value) {
    for (auto i = 0; i < this->data.size(); i++) {
        this->data[i].div(value);
    }
}

void Vector::set(const PetscScalar & value) {
    for (auto i = 0; i < this->data.size(); i++) {
        this->data[i].set(value);
    }
}

void Vector::setInternal(const PetscScalar & value){
    for (auto i = 0; i < this->data.size(); i++) {
        this->data[i].setInternal(value);
    }
}

ostream & operator<<(ostream & out, const Vector & vector) {
    out << "Vector: [" << vector.dmgPtr->rank << "]" << endl;

    for (auto i = 0; i < vector.data.size(); i++) {
        out << vector.data[i];
    }

    return out;
}

void Vector::showInFile() {
    static PetscInt count = 0;
    string fileName = "vector_data_";
    fileName += to_string(count++);
    fileName += "_";
    fileName += to_string(this->dmgPtr->rank);
    fileName += ".dat";

    ofstream out(fileName);

    out << *this;

    out.close();
}

// Scalar
Scalar add(const Scalar & left, const Scalar & right) {
    Scalar out = left;
    out.add(right);

    return out;
}

Scalar add(const Scalar & left, const PetscScalar & right) {
    Scalar out = left;
    out.add(right);

    return out;
}

Scalar add(const PetscScalar & left, const Scalar & right) {
    Scalar out = right;
    out.add(left);

    return out;
}

Scalar sub(const Scalar & left, const Scalar & right) {
    Scalar out = left;
    out.sub(right);

    return out;
}

Scalar sub(const Scalar & left, const PetscScalar & right) {
    Scalar out = left;
    out.sub(right);

    return out;
}

Scalar sub(const PetscScalar & left, const Scalar & right) {
    Scalar out = right;
    out.mult(-1);
    out.add(left);

    return out;
}

Scalar mult(const Scalar & left, const Scalar & right) {
    Scalar out = left;
    out.mult(right);

    return out;
}

Scalar mult(const Scalar & left, const PetscScalar & right) {
    Scalar out = left;
    out.mult(right);

    return out;
}

Scalar mult(const PetscScalar & left, const Scalar & right) {
    Scalar out = right;
    out.mult(left);

    return out;
}

Scalar div(const Scalar & left, const Scalar & right) {
    Scalar out = left;
    out.div(right);

    return out;
}

Scalar div(const Scalar & left, const PetscScalar & right) {
    Scalar out = left;
    out.div(right);

    return out;
}

Scalar div(const PetscScalar & left, const Scalar & right) {
    Scalar out = right;
    out.set(left);
    out.div(right);

    return out;
}

// Vector
Vector add(const Vector &left, const Vector &right) {
    Vector out = left;
    out.add(right);

    return out;
}

Vector add(const Vector & left, const PetscScalar & right) {
    Vector out = left;
    out.add(right);

    return out;
}

Vector add(const PetscScalar &left, const Vector &right) {
    Vector out = right;
    out.add(left);

    return out;
}

Vector sub(const Vector & left, const Vector & right) {
    Vector out = left;
    out.sub(right);

    return out;
}

Vector sub(const Vector & left, const PetscScalar & right) {
    Vector out = left;
    out.sub(right);

    return out;
}

Vector sub(const PetscScalar & left, const Vector & right) {
    Vector out = right;
    out.mult(-1);
    out.add(left);

    return out;
}

Vector mult(const Vector & left, const Vector & right) {
    Vector out = left;
    out.mult(right);

    return out;
}

Vector mult(const Vector & left, const PetscScalar & right) {
    Vector out = left;
    out.mult(right);

    return out;
}

Vector mult(const PetscScalar & left, const Vector & right) {
    Vector out = right;
    out.mult(left);

    return out;
}

Vector div(const Vector & left, const Vector & right) {
    Vector out = left;
    out.div(right);

    return out;
}

Vector div(const Vector & left, const PetscScalar & right) {
    Vector out = left;
    out.div(right);

    return out;
}

Vector div(const PetscScalar & left, const Vector & right) {
    Vector out = right;
    out.set(left);
    out.div(right);

    return out;
}

// Matrix
Matrix::Matrix(DataManage & dmg, MatrixType type, PetscInt matStoreNum, void *mymult, void *usrPtr) {
    this->init(dmg, type, matStoreNum, mymult, usrPtr);
}

void Matrix::init(DataManage & dmg, MatrixType type, PetscInt matStoreNum, void *mymult, void *usrPtr) {
    if (nullptr != this->data) {
        MatDestroy(&(this->data));
        this->data = nullptr;
        this->myMult = nullptr;
        this->usrPtr = nullptr;
    }

    // da
    CHECKDMDA(dmg.da);
    this->da = dmg.da;
    this->dmgPtr = &dmg;

    // 预分配内存
    this->matStoreNum = matStoreNum;
    
    // 矩阵类型
    this->type = type;
    if (NormMatrixType == this->type) {
        // 系数法
        DMCreateMatrix(this->da, &(this->data));
        MatSetFromOptions(this->data);

        // 预分配内存, 提升性能
        MatMPIAIJSetPreallocation(this->data, this->matStoreNum, NULL, this->matStoreNum, NULL);
        MatSeqAIJSetPreallocation(this->data, this->matStoreNum, NULL);
        MatSeqSBAIJSetPreallocation(this->data, 1, this->matStoreNum, NULL);
        MatMPISBAIJSetPreallocation(this->data, 1, this->matStoreNum, NULL, this->matStoreNum, NULL);
        MatMPISELLSetPreallocation(this->data, this->matStoreNum, NULL, this->matStoreNum, NULL);
        MatSeqSELLSetPreallocation(this->data, this->matStoreNum, NULL);
    }
    else if (MatrixFreeType == this->type) {
        this->myMult = (void (*)(void))mymult;
        this->usrPtr = usrPtr;

        // 创建Mat A shellMat是自定义矩阵方法，可以自定义矩阵乘法
        Scalar temp(dmg);
        PetscInt localSize;
        VecGetLocalSize(temp.data, &localSize);

        PetscInt Mx = this->dmgPtr->Mx;
        PetscInt My = this->dmgPtr->My;
        PetscInt Mz = this->dmgPtr->Mz;

        MatCreateShell(PETSC_COMM_WORLD, localSize, localSize, Mx * My * Mz,
                       Mx * My * Mz, this->usrPtr, &(this->data));
        MatShellSetOperation(this->data, MATOP_MULT, this->myMult);
        MatSetFromOptions(this->data);
    }
    else {
        cout << "matrix type error!" << endl;
        exit(-1);
    }

    this->isInitVar = true;
}

bool Matrix::isInit() {
    return this->isInitVar;
}

Matrix::~Matrix() {
    if (nullptr != this->data) {
        MatDestroy(&(this->data));
        this->data = nullptr;
        this->myMult = nullptr;
        this->usrPtr = nullptr;
    }
}

// 算子
void initOperator(PetscInt Mx, PetscInt My, PetscInt Mz) {
    baseDif.init(Mx, My, Mz);
}

ostream & operator<<(ostream & out, const BaseDif3D & baseDif3D) {
    int size, rank;
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    out << "BaseDif3D i, j, k : [" << rank << "]" << endl;
    out << "len: " << baseDif3D.len << endl;
    for (auto one : baseDif3D.indexAll) {
        for (auto io : one) {
            out << io << " ";
        }
        out << endl;
    }
    cout << endl;

    out << "lenx: " << baseDif3D.lenx << endl;
    for (auto one : baseDif3D.indexLeftPx) {
        for (auto io : one) {
            out << io << " ";
        }
        out << endl;
    }
    cout << endl;
    for (auto one : baseDif3D.indexRightPx) {
        for (auto io : one) {
            out << io << " ";
        }
        out << endl;
    }
    cout << endl;

    out << "leny: " << baseDif3D.leny << endl;
    for (auto one : baseDif3D.indexLeftPy) {
        for (auto io : one) {
            out << io << " ";
        }
        out << endl;
    }
    cout << endl;

    for (auto one : baseDif3D.indexRightPy) {
        for (auto io : one) {
            out << io << " ";
        }
        out << endl;
    }
    cout << endl;

    out << "lenz: " << baseDif3D.lenz << endl;
    for (auto one : baseDif3D.indexLeftPz) {
        for (auto io : one) {
            out << io << " ";
        }
        out << endl;
    }
    cout << endl;

    for (auto one : baseDif3D.indexRightPz) {
        for (auto io : one) {
            out << io << " ";
        }
        out << endl;
    }
    cout << endl;

    return out;
}

void printBaseOperator() {
    int rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    if (!rank) cout << baseDif;

}

void BaseDif3D::init(PetscInt Mx, PetscInt My, PetscInt Mz) {
    if (Mx > 1 || My > 1 || Mz > 1) {
        this->indexAll.clear();
        this->indexLeftPx.clear();
        this->indexLeftPy.clear();
        this->indexLeftPz.clear();
        this->indexRightPx.clear();
        this->indexRightPy.clear();
        this->indexRightPz.clear();

        vector<PetscInt> tmpOne_0 = {0, 0, 0};      // i, j, k for x, y, z
        vector<PetscInt> tmpOne_1 = {0, 0, 1};
        vector<PetscInt> tmpOne_2 = {0, 1, 0};
        vector<PetscInt> tmpOne_3 = {0, 1, 1};
        vector<PetscInt> tmpOne_4 = {1, 0, 0};
        vector<PetscInt> tmpOne_5 = {1, 0, 1};
        vector<PetscInt> tmpOne_6 = {1, 1, 0};
        vector<PetscInt> tmpOne_7 = {1, 1, 1};

        vector<vector<PetscInt>> tempOne = {tmpOne_0, tmpOne_1, tmpOne_2, tmpOne_3,
                                            tmpOne_4, tmpOne_5, tmpOne_6, tmpOne_7};

        if (Mx > 1) {
            this->indexLeftPx.push_back(tmpOne_4);
            this->indexRightPx.push_back(tmpOne_0);

            if (My > 1) {
                this->indexLeftPx.push_back(tmpOne_6);
                this->indexRightPx.push_back(tmpOne_2);
            }

            if (Mz > 1) {
                this->indexLeftPx.push_back(tmpOne_5);
                this->indexRightPx.push_back(tmpOne_1);
            }

            if (My > 1 && Mz > 1) {
                this->indexLeftPx.push_back(tmpOne_7);
                this->indexRightPx.push_back(tmpOne_3);
            }

            this->indexAll.clear();
            for (auto j = 0; j < this->indexLeftPx.size(); j++) {
                this->indexAll.push_back(this->indexLeftPx[j]);
                this->indexAll.push_back(this->indexRightPx[j]);
            }
        }

        if (My > 1)
        {
            this->indexLeftPy.push_back(tmpOne_2);
            this->indexRightPy.push_back(tmpOne_0);

            if (Mx > 1) {
                this->indexLeftPy.push_back(tmpOne_6);
                this->indexRightPy.push_back(tmpOne_4);
            }

            if (Mz > 1) {
                this->indexLeftPy.push_back(tmpOne_3);
                this->indexRightPy.push_back(tmpOne_1);
            }

            if (Mx > 1 && Mz > 1) {
                this->indexLeftPy.push_back(tmpOne_7);
                this->indexRightPy.push_back(tmpOne_5);
            }

            this->indexAll.clear();
            for (auto j = 0; j < this->indexLeftPy.size(); j++) {
                this->indexAll.push_back(this->indexLeftPy[j]);
                this->indexAll.push_back(this->indexRightPy[j]);
            }
        }

        if (Mz > 1) {
            this->indexLeftPz.push_back(tmpOne_1);
            this->indexRightPz.push_back(tmpOne_0);

            if (Mx > 1) {
                this->indexLeftPz.push_back(tmpOne_5);
                this->indexRightPz.push_back(tmpOne_4);
            }

            if (My > 1) {
                this->indexLeftPz.push_back(tmpOne_3);
                this->indexRightPz.push_back(tmpOne_2);
            }

            if (Mx > 1 && My > 1) {
                this->indexLeftPz.push_back(tmpOne_7);
                this->indexRightPz.push_back(tmpOne_6);
            }

            this->indexAll.clear();
            for (auto j = 0; j < this->indexLeftPz.size(); j++) {
                this->indexAll.push_back(this->indexLeftPz[j]);
                this->indexAll.push_back(this->indexRightPz[j]);
            }
        }
    
        this->len = this->indexAll.size();
        this->lenx = this->indexLeftPx.size();
        this->leny = this->indexLeftPy.size();
        this->lenz = this->indexLeftPz.size();
    }
}

inline PetscScalar BaseDif3D::px(const PetscScalar ***dataPtr, PetscInt i, PetscInt j, PetscInt k, PetscScalar dx) {
    PetscScalar out = 0.0;
    PetscInt    ileft, jleft, kleft, iright, jright, kright;

    if (this->lenx > 0) {
        for (auto index = 0; index < this->lenx; index++) {
            ileft = i + this->indexLeftPx[index][0];
            jleft = j + this->indexLeftPx[index][1];
            kleft = k + this->indexLeftPx[index][2];
            iright = i + this->indexRightPx[index][0];
            jright = j + this->indexRightPx[index][1];
            kright = k + this->indexRightPx[index][2];
            out += dataPtr[ileft][jleft][kleft] - dataPtr[iright][jright][kright];
        }

        out /= this->lenx * dx;
    }

    return out ;
}

inline PetscScalar BaseDif3D::py(const PetscScalar ***dataPtr, PetscInt i, PetscInt j, PetscInt k, PetscScalar dy) {
    PetscScalar out = 0.0;
    PetscInt    ileft, jleft, kleft, iright, jright, kright;

    if (this->leny > 0) {
        for (auto index = 0; index < this->leny; index++) {
            ileft = i + this->indexLeftPy[index][0];
            jleft = j + this->indexLeftPy[index][1];
            kleft = k + this->indexLeftPy[index][2];
            iright = i + this->indexRightPy[index][0];
            jright = j + this->indexRightPy[index][1];
            kright = k + this->indexRightPy[index][2];
            out += dataPtr[ileft][jleft][kleft] - dataPtr[iright][jright][kright];
        }

        out /= this->leny * dy;
    }

    return out ;
}

inline PetscScalar BaseDif3D::pz(const PetscScalar ***dataPtr, PetscInt i, PetscInt j, PetscInt k, PetscScalar dz) {
    PetscScalar out = 0.0;
    PetscInt    ileft, jleft, kleft, iright, jright, kright;

    if (this->lenz > 0) {
        for (auto index = 0; index < this->lenz; index++) {
            ileft = i + this->indexLeftPz[index][0];
            jleft = j + this->indexLeftPz[index][1];
            kleft = k + this->indexLeftPz[index][2];
            iright = i + this->indexRightPz[index][0];
            jright = j + this->indexRightPz[index][1];
            kright = k + this->indexRightPz[index][2];
            out += dataPtr[ileft][jleft][kleft] - dataPtr[iright][jright][kright];
        }

        out /= this->lenz * dz;
    }

    return out ;
}

inline PetscScalar BaseDif3D::sum(const PetscScalar ***dataPtr, PetscInt i, PetscInt j, PetscInt k) {
    PetscScalar out = 0.0;
    PetscInt    ii, jj, kk;

    if (this->len > 0) {
        for (auto index : this->indexAll) {
            ii = i + index[0];
            jj = j + index[1];
            kk = k + index[2];
            out += dataPtr[ii][jj][kk];
        }

        out /= this->len;
    }

    return out ;
}

inline PetscInt getStartEnd(PetscInt M, PetscInt & start, PetscInt & end, PetscInt offset) {
    PetscInt off = 0;

    if (M > 1 && offset > 0) {
        off = 1 - offset % 2;
        end = (start + end == M) ? end - offset : end;
        start += static_cast<PetscInt> ((offset-1) / 2.0);
    }

    return off;
}

Scalar getInterpStandard(const Scalar &scalarInput) {
    Scalar              scalarOutput(*(scalarInput.dmgPtr));
    Vec                 scalarLocal, interpLocal;
    PetscInt            i, j, k;

    const PetscScalar   ***scalarPtr;
    PetscScalar         ***interpPtr;

    PetscInt            Mx = scalarInput.dmgPtr->Mx;
    PetscInt            My = scalarInput.dmgPtr->My;
    PetscInt            Mz = scalarInput.dmgPtr->Mz;
    PetscInt            xstart = scalarInput.dmgPtr->xstart;
    PetscInt            ystart = scalarInput.dmgPtr->ystart;
    PetscInt            zstart = scalarInput.dmgPtr->zstart;
    PetscInt            xend = scalarInput.dmgPtr->xend;
    PetscInt            yend = scalarInput.dmgPtr->yend;
    PetscInt            zend = scalarInput.dmgPtr->zend;

    // 输出向量置0, 防止错误累加
    scalarOutput.set(0.0);
    scalarOutput.offset = scalarInput.offset + 1;

    // check
    CHECKDMDA(scalarInput.da)
    CHECKDMDA(scalarOutput.da)

    // 创建本地矩阵，并将全局数据映射到本地
    DMCreateLocalVector(scalarInput.da, &scalarLocal);
    DMGlobalToLocal(scalarInput.da, scalarInput.data, INSERT_VALUES, scalarLocal);

    DMCreateLocalVector(scalarOutput.da, &interpLocal);
    DMGlobalToLocal(scalarOutput.da, scalarOutput.data, INSERT_VALUES, interpLocal);

    // 获取Vec数据指针，需要注意的是输入数据均采用Read模式
    DMDAVecGetArrayRead(scalarInput.da, scalarLocal, &scalarPtr);
    DMDAVecGetArray(scalarOutput.da, interpLocal, &interpPtr);

    int offx = getStartEnd(Mx, xstart, xend, scalarOutput.offset);
    int offy = getStartEnd(My, ystart, yend, scalarOutput.offset);
    int offz = getStartEnd(Mz, zstart, zend, scalarOutput.offset);

    for (i = xstart; i < xstart + xend; i++) {
        for (j = ystart; j < ystart + yend; j++) {
            for (k = zstart; k < zstart + zend; k++) {
                interpPtr[i+offx][j+offy][k+offz] = baseDif.sum(scalarPtr, i, j, k);
            }
        }
    }

    DMDAVecRestoreArrayRead(scalarInput.da, scalarLocal, &scalarPtr);
    DMDAVecRestoreArray(scalarOutput.da, interpLocal, &interpPtr);

    DMLocalToGlobal(scalarOutput.da, interpLocal, ADD_VALUES, scalarOutput.data);

    VecDestroy(&scalarLocal);
    VecDestroy(&interpLocal);

    return scalarOutput;
}

Vector getInterpStandard(const Vector &vectorInput) {
    Vector              vectorOutput(3, *vectorInput.dmgPtr);

    vectorOutput.data[0] = getInterpStandard(vectorInput.data[0]);
    vectorOutput.data[1] = getInterpStandard(vectorInput.data[1]);
    vectorOutput.data[2] = getInterpStandard(vectorInput.data[2]);

    return vectorOutput;
}

Vector getGradStandard(const Scalar &scalarInput) {
    Vector                  vectorOutput(3, *scalarInput.dmgPtr);
    Vec                     scalarLocal;
    Vec                     vectorLocalX, vectorLocalY, vectorLocalZ;
    PetscInt                i, j, k;
    PetscInt                cod = scalarInput.dmgPtr->cod;

    PetscInt                Mx = scalarInput.dmgPtr->Mx;
    PetscInt                My = scalarInput.dmgPtr->My;
    PetscInt                Mz = scalarInput.dmgPtr->Mz;
    PetscInt                xstart = scalarInput.dmgPtr->xstart;
    PetscInt                ystart = scalarInput.dmgPtr->ystart;
    PetscInt                zstart = scalarInput.dmgPtr->zstart;
    PetscInt                xend = scalarInput.dmgPtr->xend;
    PetscInt                yend = scalarInput.dmgPtr->yend;
    PetscInt                zend = scalarInput.dmgPtr->zend;

    PetscScalar             dx = scalarInput.dmgPtr->dx;
    PetscScalar             dy = scalarInput.dmgPtr->dy;
    PetscScalar             dz = scalarInput.dmgPtr->dz;
    PetscScalar             invx;

    const PetscScalar       ***scalarPtr;
    PetscScalar             ***vectorPtrX, ***vectorPtrY, ***vectorPtrZ;

    // 输出向量置0, 防止错误累加
    vectorOutput.set(0.0);
    vectorOutput.data[0].offset = scalarInput.offset + 1;
    vectorOutput.data[1].offset = scalarInput.offset + 1;
    vectorOutput.data[2].offset = scalarInput.offset + 1;

    // check
    CHECKDMDA(scalarInput.da)
    CHECKDMDA(vectorOutput.da)

    // 创建本地矩阵，并将全局数据映射到本地(映射本地后, 可以越本地进程数据范围访问数据, 所越行数即即为DMDACreate3d的s参数)
    DMCreateLocalVector(scalarInput.da, &scalarLocal);
    DMGlobalToLocal(scalarInput.da, scalarInput.data, INSERT_VALUES, scalarLocal);

    DMCreateLocalVector(vectorOutput.da, &vectorLocalX);
    DMGlobalToLocal(vectorOutput.da, vectorOutput.data[0].data, INSERT_VALUES, vectorLocalX);

    DMCreateLocalVector(vectorOutput.da, &vectorLocalY);
    DMGlobalToLocal(vectorOutput.da, vectorOutput.data[1].data, INSERT_VALUES, vectorLocalY);

    DMCreateLocalVector(vectorOutput.da, &vectorLocalZ);
    DMGlobalToLocal(vectorOutput.da, vectorOutput.data[2].data, INSERT_VALUES, vectorLocalZ);

    // 获取Vec数据指针，需要注意的是输入数据均采用Read模式
    DMDAVecGetArrayRead(scalarInput.da, scalarLocal, &scalarPtr);
    DMDAVecGetArray(vectorOutput.da, vectorLocalX, &vectorPtrX);
    DMDAVecGetArray(vectorOutput.da, vectorLocalY, &vectorPtrY);
    DMDAVecGetArray(vectorOutput.da, vectorLocalZ, &vectorPtrZ);

    int offx = getStartEnd(Mx, xstart, xend, vectorOutput.data[0].offset);
    int offy = getStartEnd(My, ystart, yend, vectorOutput.data[0].offset);
    int offz = getStartEnd(Mz, zstart, zend, vectorOutput.data[0].offset);

    if (0 == cod) {                         // 直角坐标系
        for (i = xstart; i < xstart + xend; i++) {
            for (j = ystart; j < ystart + yend; j++) {
                for (k = zstart; k < zstart + zend; k++) {
                    vectorPtrX[i+offx][j+offy][k+offz] = baseDif.px(scalarPtr, i, j, k, dx);
                    vectorPtrY[i+offx][j+offy][k+offz] = baseDif.py(scalarPtr, i, j, k, dy);
                    vectorPtrZ[i+offx][j+offy][k+offz] = baseDif.pz(scalarPtr, i, j, k, dz);
                }
            }
        }
    }
    else if (1 == cod) {                    // 柱坐标系
        for (i = xstart; i < xstart + xend; i++) {
            for (j = ystart; j < ystart + yend; j++) {
                for (k = zstart; k < zstart + zend; k++) {
                    invx = 1.0 / (scalarInput.dmgPtr->x0 + dx * (i + 0.5));

                    vectorPtrX[i+offx][j+offy][k+offz] = baseDif.px(scalarPtr, i, j, k, dx);
                    vectorPtrY[i+offx][j+offy][k+offz] = baseDif.py(scalarPtr, i, j, k, dy) * invx;
                    vectorPtrZ[i+offx][j+offy][k+offz] = baseDif.pz(scalarPtr, i, j, k, dz);
                }
            }
        }
    }

    DMDAVecRestoreArrayRead(scalarInput.da, scalarLocal, &scalarPtr);
    DMDAVecRestoreArray(vectorOutput.da, vectorLocalX, &vectorPtrX);
    DMDAVecRestoreArray(vectorOutput.da, vectorLocalY, &vectorPtrY);
    DMDAVecRestoreArray(vectorOutput.da, vectorLocalZ, &vectorPtrZ);

    // 将本地数据映射回全局
    DMLocalToGlobal(vectorOutput.da, vectorLocalX, ADD_VALUES, vectorOutput.data[0].data);
    DMLocalToGlobal(vectorOutput.da, vectorLocalY, ADD_VALUES, vectorOutput.data[1].data);
    DMLocalToGlobal(vectorOutput.da, vectorLocalZ, ADD_VALUES, vectorOutput.data[2].data);

    VecDestroy(&scalarLocal);
    VecDestroy(&vectorLocalX);
    VecDestroy(&vectorLocalY);
    VecDestroy(&vectorLocalZ);

    return vectorOutput;
}

Scalar getDivStandard(const Vector &vectorInput) {
    Scalar                  scalarOutput(*vectorInput.dmgPtr);
    Vec                     vectorLocalX, vectorLocalY, vectorLocalZ, divLocal;
    PetscInt                i, j, k;
    PetscScalar             invx, vr;
    PetscInt                cod = vectorInput.dmgPtr->cod;

    PetscInt                Mx = vectorInput.dmgPtr->Mx;
    PetscInt                My = vectorInput.dmgPtr->My;
    PetscInt                Mz = vectorInput.dmgPtr->Mz;
    PetscInt                xstart = vectorInput.dmgPtr->xstart;
    PetscInt                ystart = vectorInput.dmgPtr->ystart;
    PetscInt                zstart = vectorInput.dmgPtr->zstart;
    PetscInt                xend = vectorInput.dmgPtr->xend;
    PetscInt                yend = vectorInput.dmgPtr->yend;
    PetscInt                zend = vectorInput.dmgPtr->zend;

    PetscScalar             dx = vectorInput.dmgPtr->dx;
    PetscScalar             dy = vectorInput.dmgPtr->dy;
    PetscScalar             dz = vectorInput.dmgPtr->dz;

    PetscScalar             ***divPtr;
    const PetscScalar       ***vectorPtrX, ***vectorPtrY, ***vectorPtrZ;

    // 输出向量置0, 防止错误累加
    scalarOutput.set(0.0);
    scalarOutput.offset = vectorInput.data[0].offset + 1;

    // check
    CHECKDMDA(vectorInput.da)
    CHECKDMDA(scalarOutput.da)

    // 创建本地矩阵，并将全局数据映射到本地
    DMCreateLocalVector(vectorInput.da, &vectorLocalX);
    DMGlobalToLocal(vectorInput.da, vectorInput.data[0].data, INSERT_VALUES, vectorLocalX);

    DMCreateLocalVector(vectorInput.da, &vectorLocalY);
    DMGlobalToLocal(vectorInput.da, vectorInput.data[1].data, INSERT_VALUES, vectorLocalY);

    DMCreateLocalVector(vectorInput.da, &vectorLocalZ);
    DMGlobalToLocal(vectorInput.da, vectorInput.data[2].data, INSERT_VALUES, vectorLocalZ);

    DMCreateLocalVector(scalarOutput.da, &divLocal);
    DMGlobalToLocal(scalarOutput.da, scalarOutput.data, INSERT_VALUES, divLocal);

    // 获取Vec数据指针，需要注意的是输入数据均采用Read模式
    DMDAVecGetArrayRead(vectorInput.da, vectorLocalX, &vectorPtrX);
    DMDAVecGetArrayRead(vectorInput.da, vectorLocalY, &vectorPtrY);
    DMDAVecGetArrayRead(vectorInput.da, vectorLocalZ, &vectorPtrZ);
    DMDAVecGetArray(scalarOutput.da, divLocal, &divPtr);

    int offx = getStartEnd(Mx, xstart, xend, scalarOutput.offset);
    int offy = getStartEnd(My, ystart, yend, scalarOutput.offset);
    int offz = getStartEnd(Mz, zstart, zend, scalarOutput.offset);

    if (0 == cod) {                     // 直角坐标系
        for (i = xstart; i < xstart + xend; i++) {
            for (j = ystart; j < ystart + yend; j++) {
                for (k = zstart; k < zstart + zend; k++) {
                    divPtr[i+offx][j+offy][k+offz] = baseDif.px(vectorPtrX, i, j, k, dx) +
                                      baseDif.py(vectorPtrY, i, j, k, dy) +
                                      baseDif.pz(vectorPtrZ, i, j, k, dz);
                }
            }
        }
    }
    else if (1 == cod) {                // 柱坐标系
        for (i = xstart; i < xstart + xend; i++) {
            for (j = ystart; j < ystart + yend; j++) {
                for (k = zstart; k < zstart + zend; k++) {
                    invx = 1.0 / (vectorInput.dmgPtr->x0 + dx * (i + 0.5));
                    vr = baseDif.sum(vectorPtrX, i, j, k);

                    divPtr[i+offx][j+offy][k+offz] = baseDif.px(vectorPtrX, i, j, k, dx) + vr*invx +
                                      baseDif.py(vectorPtrY, i, j, k, dy) * invx +
                                      baseDif.pz(vectorPtrZ, i, j, k, dz);
                }
            }
        }
    }

    DMDAVecRestoreArrayRead(vectorInput.da, vectorLocalX, &vectorPtrX);
    DMDAVecRestoreArrayRead(vectorInput.da, vectorLocalY, &vectorPtrY);
    DMDAVecRestoreArrayRead(vectorInput.da, vectorLocalZ, &vectorPtrZ);
    DMDAVecRestoreArray(scalarOutput.da, divLocal, &divPtr);

    DMLocalToGlobal(scalarOutput.da, divLocal, ADD_VALUES, scalarOutput.data);

    VecDestroy(&vectorLocalX);
    VecDestroy(&vectorLocalY);
    VecDestroy(&vectorLocalZ);
    VecDestroy(&divLocal);

    return scalarOutput;
}

Vector getCurlStandard(const Vector & vectorInput) {
    Vector                  vectorOutput(3, *vectorInput.dmgPtr);
    Vec                     vectorLocalX, vectorLocalY, vectorLocalZ;
    Vec                     curlLocalX, curlLocalY, curlLocalZ;
    PetscInt                i, j, k;
    PetscScalar             invx, vth;
    PetscInt                cod = vectorInput.dmgPtr->cod;

    PetscInt                Mx = vectorInput.dmgPtr->Mx;
    PetscInt                My = vectorInput.dmgPtr->My;
    PetscInt                Mz = vectorInput.dmgPtr->Mz;
    PetscInt                xstart = vectorInput.dmgPtr->xstart;
    PetscInt                ystart = vectorInput.dmgPtr->ystart;
    PetscInt                zstart = vectorInput.dmgPtr->zstart;
    PetscInt                xend = vectorInput.dmgPtr->xend;
    PetscInt                yend = vectorInput.dmgPtr->yend;
    PetscInt                zend = vectorInput.dmgPtr->zend;

    PetscScalar             dx = vectorInput.dmgPtr->dx;
    PetscScalar             dy = vectorInput.dmgPtr->dy;
    PetscScalar             dz = vectorInput.dmgPtr->dz;

    PetscScalar             ***curlPtrX, ***curlPtrY, ***curlPtrZ;
    const PetscScalar       ***vectorPtrX, ***vectorPtrY, ***vectorPtrZ;

    // 输出向量置0, 防止错误累加
    vectorOutput.set(0.0);
    vectorOutput.data[0].offset = vectorInput.data[0].offset + 1;
    vectorOutput.data[1].offset = vectorInput.data[1].offset + 1;
    vectorOutput.data[2].offset = vectorInput.data[2].offset + 1;

    // check
    CHECKDMDA(vectorInput.da)
    CHECKDMDA(vectorOutput.da)

    // 创建本地矩阵，并将全局数据映射到本地
    DMCreateLocalVector(vectorInput.da, &vectorLocalX);
    DMGlobalToLocal(vectorInput.da, vectorInput.data[0].data, INSERT_VALUES, vectorLocalX);

    DMCreateLocalVector(vectorInput.da, &vectorLocalY);
    DMGlobalToLocal(vectorInput.da, vectorInput.data[1].data, INSERT_VALUES, vectorLocalY);

    DMCreateLocalVector(vectorInput.da, &vectorLocalZ);
    DMGlobalToLocal(vectorInput.da, vectorInput.data[2].data, INSERT_VALUES, vectorLocalZ);

    DMCreateLocalVector(vectorOutput.da, &curlLocalX);
    DMGlobalToLocal(vectorOutput.da, vectorOutput.data[0].data, INSERT_VALUES, curlLocalX);

    DMCreateLocalVector(vectorOutput.da, &curlLocalY);
    DMGlobalToLocal(vectorOutput.da, vectorOutput.data[1].data, INSERT_VALUES, curlLocalY);

    DMCreateLocalVector(vectorOutput.da, &curlLocalZ);
    DMGlobalToLocal(vectorOutput.da, vectorOutput.data[2].data, INSERT_VALUES, curlLocalZ);

    // 获取Vec数据指针，需要注意的是输入数据均采用Read模式
    DMDAVecGetArrayRead(vectorInput.da, vectorLocalX, &vectorPtrX);
    DMDAVecGetArrayRead(vectorInput.da, vectorLocalY, &vectorPtrY);
    DMDAVecGetArrayRead(vectorInput.da, vectorLocalZ, &vectorPtrZ);
    DMDAVecGetArray(vectorOutput.da, curlLocalX, &curlPtrX);
    DMDAVecGetArray(vectorOutput.da, curlLocalY, &curlPtrY);
    DMDAVecGetArray(vectorOutput.da, curlLocalZ, &curlPtrZ);

    int offx = getStartEnd(Mx, xstart, xend, vectorOutput.data[0].offset);
    int offy = getStartEnd(My, ystart, yend, vectorOutput.data[0].offset);
    int offz = getStartEnd(Mz, zstart, zend, vectorOutput.data[0].offset);

    if (0 == cod) {                     // 直角坐标系
        for (i = xstart; i < xstart + xend; i++) {
            for (j = ystart; j < ystart + yend; j++) {
                for (k = zstart; k < zstart + zend; k++) {
                    curlPtrX[i+offx][j+offy][k+offz] = baseDif.py(vectorPtrZ, i, j, k, dy) - baseDif.pz(vectorPtrY, i, j, k, dz);
                    curlPtrY[i+offx][j+offy][k+offz] = baseDif.pz(vectorPtrX, i, j, k, dz) - baseDif.px(vectorPtrZ, i, j, k, dx);
                    curlPtrZ[i+offx][j+offy][k+offz] = baseDif.px(vectorPtrY, i, j, k, dx) - baseDif.py(vectorPtrX, i, j, k, dy);
                }
            }
        }
    }
    else if (1 == cod) {                // 柱坐标系
        for (i = xstart; i < xstart + xend; i++) {
            for (j = ystart; j < ystart + yend; j++) {
                for (k = zstart; k < zstart + zend; k++) {
                    invx = 1.0 / (vectorInput.dmgPtr->x0 + dx * (i + 0.5));
                    vth = baseDif.sum(vectorPtrY, i, j, k);

                    curlPtrX[i+offx][j+offy][k+offz] = baseDif.py(vectorPtrZ, i, j, k, dy) * invx - baseDif.pz(vectorPtrY, i, j, k, dz);
                    curlPtrY[i+offx][j+offy][k+offz] = baseDif.pz(vectorPtrX, i, j, k, dz) - baseDif.px(vectorPtrZ, i, j, k, dx);
                    curlPtrZ[i+offx][j+offy][k+offz] = vth * invx + baseDif.px(vectorPtrY, i, j, k, dx) - baseDif.py(vectorPtrX, i, j, k, dy) * invx;
                }
            }
        }
    }

    DMDAVecRestoreArrayRead(vectorInput.da, vectorLocalX, &vectorPtrX);
    DMDAVecRestoreArrayRead(vectorInput.da, vectorLocalY, &vectorPtrY);
    DMDAVecRestoreArrayRead(vectorInput.da, vectorLocalZ, &vectorPtrZ);
    DMDAVecRestoreArray(vectorOutput.da, curlLocalX, &curlPtrX);
    DMDAVecRestoreArray(vectorOutput.da, curlLocalY, &curlPtrY);
    DMDAVecRestoreArray(vectorOutput.da, curlLocalZ, &curlPtrZ);

    DMLocalToGlobal(vectorOutput.da, curlLocalX, ADD_VALUES, vectorOutput.data[0].data);
    DMLocalToGlobal(vectorOutput.da, curlLocalY, ADD_VALUES, vectorOutput.data[1].data);
    DMLocalToGlobal(vectorOutput.da, curlLocalZ, ADD_VALUES, vectorOutput.data[2].data);

    VecDestroy(&vectorLocalX);
    VecDestroy(&vectorLocalY);
    VecDestroy(&vectorLocalZ);
    VecDestroy(&curlLocalX);
    VecDestroy(&curlLocalY);
    VecDestroy(&curlLocalZ);

    return vectorOutput;
}

Vector getGrad(const Scalar &scalarInput) {
    return getInterpStandard(getGradStandard(scalarInput));
}

Scalar getDiv(const Vector &vectorInput) {
    return getInterpStandard(getDivStandard(vectorInput));
}

Vector getCurl(const Vector & vectorInput) {
    return getInterpStandard(getCurlStandard(vectorInput));
}

Scalar getLap(const Scalar &scalarInput) {
    return getDivStandard(getGradStandard(scalarInput));
}

Vector getLap(const Vector & vectorInput) {
    return sub(getGradStandard(getDivStandard(vectorInput)), getCurlStandard(getCurlStandard(vectorInput)));
}

Vector getTensorDot(const Tensor3d &tensorInput, const Vector &vectorInput) {
    Vector                  vectorOutput(3, *vectorInput.dmgPtr);

    // check
    CHECKDMDA(tensorInput.xx.da)

    vectorOutput.data[0].set(0.0);
    vectorOutput.data[0].add(mult(tensorInput.xx, vectorInput.data[0]));
    vectorOutput.data[0].add(mult(tensorInput.xy, vectorInput.data[1]));
    vectorOutput.data[0].add(mult(tensorInput.xz, vectorInput.data[2]));

    vectorOutput.data[1].set(0.0);
    vectorOutput.data[1].add(mult(tensorInput.yx, vectorInput.data[0]));
    vectorOutput.data[1].add(mult(tensorInput.yy, vectorInput.data[1]));
    vectorOutput.data[1].add(mult(tensorInput.yz, vectorInput.data[2]));

    vectorOutput.data[2].set(0.0);
    vectorOutput.data[2].add(mult(tensorInput.zx, vectorInput.data[0]));
    vectorOutput.data[2].add(mult(tensorInput.zy, vectorInput.data[1]));
    vectorOutput.data[2].add(mult(tensorInput.zz, vectorInput.data[2]));

    return vectorOutput;
}

void pyplot(string scriptName, string outFileName, Vec src, PetscInt xsize, PetscInt ysize) {
    // 文件写入
    PetscViewer myViewer;
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, outFileName.c_str(), &myViewer);
    VecView(src, myViewer);
    
    PetscMPIInt rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    if (!(rank)) {
        ostringstream os;
        os << scriptName << " " << outFileName << " " << xsize << " " << ysize;
        string cmd = os.str();
        system(cmd.c_str());
    }
}

