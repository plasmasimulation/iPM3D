#include"material.h"
using namespace std;
// class Material
// { public:
//     double epsilon ; //介电常数
//     Material();
// };
Material::Material(){
    epsilon=8.8542E-12;
}

// class Metal : public Material
// { public:
//     double charge;
//     double charge_onestep;
//     double capacitance;
//     double voltage;
//     Metal();
// };
Metal::Metal(){
    epsilon=8.8542E-12;
    charge=0;
    charge_onestep=0;
    capacitance=0;
    voltage=0;
}

// class Dielectric  :public Material
// {public:
//     double permittivity;
//     double permeability;
//     double conductivity;
//     Dielectric();
// };
Dielectric::Dielectric(){
    epsilon=8.8542E-12*2; //假定介质相对介电常数为2
    permittivity=0;
    permeability=0;
    conductivity=0;
}