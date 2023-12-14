
class Material
{ public:
    double epsilon ; //介电常数
    Material();
};
class Metal : public Material
{public:
     double charge;
    double charge_onestep;
    double capacitance;
    double voltage;
    Metal();
};

class  Dielectric :public Material
{public:
    double permittivity;
    double permeability;
    double conductivity;
    Dielectric();
};
 void load_material(int data[5][5][5]);
 int create_material_file();
