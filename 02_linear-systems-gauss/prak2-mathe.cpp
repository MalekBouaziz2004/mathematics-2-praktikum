#include <iostream>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <vector>


//praktikum 1
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

    friend CMyVektor operator*(CMyVektor  a, const double lambda) {
        size_t d = a.get_DV();
        CMyVektor result(d);
        for (size_t i = 0; i < d; ++i) {
            result[i] = lambda*a[i];
        }
        return result;
    }

}

;


double norm(CMyVektor v) {
    double sum = 0.0;
    for (size_t i = 0; i < v.get_DV(); i++) {
        sum += v[i] * v[i];
    }
    return sqrt(sum);
}


// praktikum 2
class CMyMatrix{
    private:
    size_t dim_spalte,dim_zeilen;
    std::vector<std::vector<double>> eintrag;
    public:
    CMyMatrix(size_t ligne, size_t colonne){
        dim_spalte=colonne;
        dim_zeilen=ligne;
        eintrag.resize(dim_zeilen);
        for (size_t i = 0; i < dim_zeilen; ++i) {
    eintrag[i].resize(dim_spalte, 0.0);}

        
    }
    //get_zeilen
    size_t get_zeilen() const { return dim_zeilen; }
    //get_spalte
    size_t get_spalte() const { return dim_spalte; }
    //setter-komponent
    void set_K(size_t ligne, size_t colonne,double eingabe){eintrag[ligne][colonne]=eingabe;}//ligne < dim_zeilen && colonne < dim_spalt
    //getter-komponent
    double get_K(size_t ligne, size_t colonne) const{return eintrag[ligne][colonne];}
    //Methode CMyMatrix invers()
    CMyMatrix invers(){
        if(dim_zeilen != 2 || dim_spalte != 2){
            std::cout<< "Fehlermeldung, Muss 2*2 Matrix geben"<<std::endl;
            exit(1);
        }else{
            double a = eintrag[0][0];
            double b=eintrag[0][1];
            double c= eintrag[1][0];
            double d= eintrag[1][1];
            double determinant=a * d - b * c;
            if(determinant==0){
                std::cout<<"erreur"<<std::endl;
                exit(1);
            }else{
                CMyMatrix inversion(2,2);
                inversion.set_K(0,0,d/determinant);
                inversion.set_K(0,1,-b/determinant);
                inversion.set_K(1,0,-c/determinant);
                inversion.set_K(1,1,a/determinant);
                 
                 return inversion;

            }

        }
        

    }

   friend CMyVektor operator*(CMyMatrix A, CMyVektor x) {
    CMyVektor erg(A.dim_zeilen);

    size_t i = 0;
    while (i < A.dim_zeilen) {
        size_t k = 0;
        while (k < A.dim_spalte) {
            erg[i] += A.eintrag[i][k] * x[k];
            k += 1;
        }
        i += 1;
    }
    return erg;
}


};


CMyMatrix jacobi(CMyVektor x, CMyVektor (*funktion)(CMyVektor x)) {
    size_t m = funktion(x).get_DV(); 

    size_t n = x.get_DV();           

    CMyMatrix j(m, n);
    CMyVektor f_x = funktion(x);

    double h = 1e-4;

    for (size_t j_index = 0; j_index < n; j_index++) {
        CMyVektor x_nv = x; // kopy von x
        x_nv[j_index] += h; // wacklei

        CMyVektor f_xv = funktion(x_nv); 

        for (size_t i = 0; i < m; i++) {
            double ableitung = (f_xv[i] - f_x[i]) / h;
            j.set_K(i, j_index, ableitung);
        }
    }

    return j;
}

CMyVektor newton_verfahren(CMyVektor x, CMyVektor (*funktion)(CMyVektor)) {
    size_t schritt = 0;

    while (schritt < 50) {
        CMyVektor f_x = funktion(x);

        std::cout << "Schritt " << schritt << ":\n";

        std::cout << "\t x = ( ";
        for (size_t i = 0; i < x.get_DV(); i++) {
            std::cout << std::fixed << std::setprecision(6) << x[i];
            if (i < x.get_DV() - 1) std::cout << "; ";
        }
        std::cout << " )\n";

        std::cout << "\t f(x) = ( ";
        for (size_t i = 0; i < f_x.get_DV(); i++) {
            std::cout << std::fixed << std::setprecision(6) << f_x[i];
            if (i < f_x.get_DV() - 1) std::cout << "; ";
        }
        std::cout << " )\n";

        CMyMatrix J = jacobi(x, funktion);
        CMyMatrix J_inv = J.invers();

        std::cout << "\t f'(x) = \n";
        for (size_t i = 0; i < J.get_zeilen(); i++) {
            std::cout << "\t\t( ";
            for (size_t j = 0; j < J.get_spalte(); j++) {
                std::cout << std::fixed << std::setprecision(6) << J.get_K(i, j);
                if (j < J.get_spalte() - 1) std::cout << "; ";
            }
            std::cout << " )\n";
        }

        std::cout << "\t (f'(x))^(-1) = \n";
        for (size_t i = 0; i < J_inv.get_zeilen(); i++) {
            std::cout << "\t\t( ";
            for (size_t j = 0; j < J_inv.get_spalte(); j++) {
                std::cout << std::fixed << std::setprecision(6) << J_inv.get_K(i, j);
                if (j < J_inv.get_spalte() - 1) std::cout << "; ";
            }
            std::cout << " )\n";
        }

        CMyVektor delta_x = J_inv * f_x * (-1);

        std::cout << "\t dx = ( ";
        for (size_t i = 0; i < delta_x.get_DV(); i++) {
            std::cout << std::fixed << std::setprecision(6) << delta_x[i];
            if (i < delta_x.get_DV() - 1) std::cout << "; ";
        }
        std::cout << " )\n";

        std::cout << "\t ||f(x)|| = " << std::fixed << std::setprecision(6) << norm(f_x) << "\n\n";

        x = x + delta_x;

       
        CMyVektor f_neu = funktion(x);
        if (norm(f_neu) < 1e-5) {
            std::cout << "Ende wegen ||f(x)||<1e-5 bei\n";
            std::cout << "\t x = ( ";
            for (size_t i = 0; i < x.get_DV(); i++) {
                std::cout << std::fixed << std::setprecision(6) << x[i];
                if (i < x.get_DV() - 1) std::cout << "; ";
            }
            std::cout << " )\n";
            std::cout << "\t f(x) = ( ";
            for (size_t i = 0; i < f_neu.get_DV(); i++) {
                std::cout << std::fixed << std::setprecision(6) << f_neu[i];
                if (i < f_neu.get_DV() - 1) std::cout << "; ";
            }
            std::cout << " )\n";
            std::cout << "\t ||f(x)|| = " << std::fixed << std::setprecision(10)<< norm(f_neu) << "\n";
            return x;
        }

        schritt++;
    }

    std::cout << "Maximale Schritte erreicht. Ergebnis:\n";
    CMyVektor f_x = funktion(x);
    std::cout << "\t x = ( ";
    for (size_t i = 0; i < x.get_DV(); i++) {
        std::cout << std::fixed << std::setprecision(6) << x[i];
        if (i < x.get_DV() - 1) std::cout << "; ";
    }
    std::cout << " )\n";
    std::cout << "\t f(x) = ( ";
    for (size_t i = 0; i < f_x.get_DV(); i++) {
        std::cout << std::fixed << std::setprecision(6) << f_x[i];
        if (i < f_x.get_DV() - 1) std::cout << "; ";
    }
    std::cout << " )\n";
    std::cout << "\t ||f(x)|| = " << std::fixed << std::setprecision(6) << norm(f_x) << "\n";

    return x;
}




CMyVektor f(CMyVektor x) {
    CMyVektor result(3);
    result[0] = x[0] * x[1] * exp(x[2]);
    result[1] = x[1] * x[2] * x[3];
    result[2] = x[3];
    return result;
}
CMyVektor g(CMyVektor OM ) {
    CMyVektor erg(2);
    double x = OM[0];
    double y = OM[1];

    erg[0] = pow(x, 3) * pow(y, 3) - 2 * y;
    erg[1] = x - 2;

    return erg;
}
int main() {
    CMyVektor x(4);
    x[0] = 1;
    x[1] = 2;
    x[2] = 0;
    x[3] = 3;
    CMyMatrix J = jacobi(x, f);
    for (size_t i = 0; i < J.get_zeilen(); i++) {
    for (size_t j = 0; j < J.get_spalte(); j++) {
        std::cout << std::fixed << std::setprecision(6) << J.get_K(i, j) << "; ";
    }
    std::cout << std::endl;}

    CMyVektor OM(2);
    OM[0]=1;
    OM[1]=1;
    CMyVektor erg = newton_verfahren(OM,g);

    system("pause");
    return 0;
    
}