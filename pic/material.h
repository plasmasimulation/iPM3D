
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
 void load_material(int data[65][65][65],int coord_x, int coord_y, int coord_z,
                   int width_x, int width_y, int width_z);
 int create_material_file();
