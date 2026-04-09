#include <iostream>
#include <cmath>
#include <iomanip>
#include <sstream>

class CMyVektor {
private:
    size_t dim;
    double* valeurs;

public:
    CMyVektor(int dimension, double* data) {
        dim = dimension;
        valeurs = new double[dim];
        for (size_t i = 0; i < dim; i++) {
            valeurs[i] = data[i];
        }
    }

    CMyVektor(int dimension) {
        dim = dimension;
        valeurs = new double[dim]{};
    }




    CMyVektor(const CMyVektor& other) {
    dim = other.dim;
    valeurs = new double[dim];
    for (size_t i = 0; i < dim; ++i) {
        valeurs[i] = other.valeurs[i];
    }
}

   
    size_t get_DV() const { return dim; }

    void set_KV(size_t index, double Komponente) { valeurs[index] = Komponente; }

    double get_KV(size_t index) const { return valeurs[index]; }

    double operator[](size_t i) const { return get_KV(i); }

    double& operator[](size_t i) { return valeurs[i]; }

    friend CMyVektor operator+(const CMyVektor& a, const CMyVektor& b) {
        size_t d = a.get_DV();
        CMyVektor result(d);
        for (size_t i = 0; i < d; ++i) {
            result[i] = a[i] + b[i];
        }
        return result;

    }

    friend CMyVektor operator*(const double lambda ,CMyVektor  a) {
        size_t d = a.get_DV();
        CMyVektor result(d);
        for (size_t i = 0; i < d; ++i) {
            result[i] = lambda*a[i];
        }
        return result;
    }

}

;


double norm(CMyVektor OM) {
    double sum = 0.0;
    for (size_t i = 0; i < OM.get_DV(); ++i) {
        sum += OM[i] * OM[i];
    }
    return sqrt(sum);
}

class C_DGLSolver {
    private:
    CMyVektor (*f_DGL_System)(CMyVektor y, double x);
    double (*f_DGL_nterOrdnung)(CMyVektor y, double x);
    bool type;
    CMyVektor ableitungen(CMyVektor y, double x){
        if(type==true){
            return f_DGL_System(y,x);
        }else{
            size_t size =y.get_DV();
            CMyVektor erg(size);
            for(size_t i=0;i<size-1;i++){
                erg[i]=y[i+1];
            }
            erg[size-1]=f_DGL_nterOrdnung(y, x);
            return erg;
        }
    }
    public:
     C_DGLSolver (CMyVektor (*f)(CMyVektor y, double x)){
        
            f_DGL_System=f;
            type=true;
        }

     C_DGLSolver (double (*f)(CMyVektor y, double x)){
        f_DGL_nterOrdnung=f;
        type=false;
     } 

    void Euler_methode(double x_start, double x_end, size_t schritt, CMyVektor & y_start, bool anzeigen = true) {
    double h = (x_end - x_start) / schritt;

    for (size_t i = 0; i < schritt; i++) {
        if (anzeigen) {
            std::cout << "Schritt " << i << ":\n";
            std::cout << "\t x = " << std::fixed << std::setprecision(6) << x_start << "\n";

            std::cout << "\t y = ( ";
            for (size_t j = 0; j < y_start.get_DV(); ++j) {
                std::cout << y_start[j];
                if (j < y_start.get_DV() - 1) std::cout << "; ";
            }
            std::cout << " )\n";

            CMyVektor y_prime = ableitungen(y_start, x_start);
            std::cout << "\t y' = ( ";
            for (size_t j = 0; j < y_prime.get_DV(); ++j) {
                std::cout << y_prime[j];
                if (j < y_prime.get_DV() - 1) std::cout << "; ";
            }
            std::cout << " )\n\n";
        }

        y_start = y_start + h * ableitungen(y_start, x_start);
        x_start += h;
    }

    if (anzeigen) {
        std::cout << "Ende bei\n";
        std::cout << "\t x = " << x_start << "\n";
        std::cout << "\t y = ( ";
        for (size_t j = 0; j < y_start.get_DV(); ++j) {
            std::cout << y_start[j];
            if (j < y_start.get_DV() - 1) std::cout << "; ";
        }
        std::cout << " )\n";
    }
}



     void Heun_methode(double x_start, double x_end, size_t schritt, CMyVektor & y_start, bool anzeigen = true) {
    double h = (x_end - x_start) / schritt;

    for (size_t i = 0; i < schritt; i++) {
        if (anzeigen) {
            std::cout << "Schritt " << i << ":\n";
            std::cout << "\t x = " << std::fixed << std::setprecision(6) << x_start << "\n";

            std::cout << "\t y = ( ";
            for (size_t j = 0; j < y_start.get_DV(); ++j) {
                std::cout << y_start[j];
                if (j < y_start.get_DV() - 1) std::cout << "; ";
            }
            std::cout << " )\n";
        }

        CMyVektor y_org = y_start;
        CMyVektor y_test = y_org + h * ableitungen(y_org, x_start);
        double x_neu = x_start + h;
        CMyVektor y_mitte = 0.5 * (ableitungen(y_org, x_start) + ableitungen(y_test, x_neu));
        CMyVektor derivee_moyenne = y_mitte;

        if (anzeigen) {
            std::cout << "\t y' = ( ";
            for (size_t j = 0; j < derivee_moyenne.get_DV(); ++j) {
                std::cout << derivee_moyenne[j];
                if (j < derivee_moyenne.get_DV() - 1) std::cout << "; ";
            }
            std::cout << " )\n\n";
        }

        y_start = y_org + h * y_mitte;
        x_start += h;
    }

    if (anzeigen) {
        std::cout << "Ende bei\n";
        std::cout << "\t x = " << x_start << "\n";
        std::cout << "\t y = ( ";
        for (size_t j = 0; j < y_start.get_DV(); ++j) {
            std::cout << y_start[j];
            if (j < y_start.get_DV() - 1) std::cout << "; ";
        }
        std::cout << " )\n";
    }
}


};

CMyVektor systm(CMyVektor y ,double x){
    CMyVektor erg(2);
    erg[0]= 2*y[1]-x*y[0];
    erg[1]=y[0]*y[1]-2*pow(x,3);
    return erg;

}
double dritte_ord_sys(CMyVektor y, double x){
    return(2*x*y[1]*y[2]+2*pow(y[0],2)*y[1]);
}

    int main() {
    CMyVektor y_start(2);
    y_start[0] = 0.0;
    y_start[1] = 1.0;

    C_DGLSolver monSolveur(systm);//objek _DGLSolver

        double wert = 0.5;
        int array [] = {10, 100, 1000, 10000};

    for (size_t i=0;i<4;i++) {
        CMyVektor y(3);
        y[0] = 1.0;   
        y[1] = -1.0; 
        y[2] = 2.0;   
        C_DGLSolver solver_nte_ordnung(dritte_ord_sys);

        CMyVektor y_euler = y;
        CMyVektor y_heun = y;

        
        solver_nte_ordnung.Euler_methode(1.0, 2.0, 100, y_euler,true);
        double yEuler = y_euler[0];
        double abweichung_Euler = yEuler - wert;

        solver_nte_ordnung.Heun_methode(1.0, 2.0, 100, y_heun,true);
        double yHeun = y_heun[0];
        double abweichung_Heun = yHeun - wert;

        //std::cout << "Abweichung bei Euler bei " << array[i] << " Schritten: " << abweichung_Euler << "\n";
        //std::cout << "Abweichung bei Heun  bei " << array[i]<< " Schritten: " <<abweichung_Heun  << "\n\n";
}

    system("pause");

    return 0;
}



