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


double norm(CMyVektor OM) {
    double sum = 0.0;
    for (size_t i = 0; i < OM.get_DV(); ++i) {
        sum += OM[i] * OM[i];
    }
    return sqrt(sum);
}

CMyVektor gradient(CMyVektor OM, double (*funktion)(CMyVektor OM))
    {
        double h = 1e-8;// abweischung ^

        CMyVektor grad(OM.get_DV());

        for (size_t i = 0; i < OM.get_DV(); i++)
        {
            CMyVektor neu_OM= OM;// copy 
            neu_OM[i] +=h;
            grad[i] = (funktion(neu_OM) - funktion(OM))/h;
        }
        return grad;
    }


void gradientenVerfahren(CMyVektor OM, double (*funktion)(CMyVektor), double lambda = 1.0) {
    const double eps = 1e-5;
    const int max_steps = 25;

    for (int schritt = 0; schritt < max_steps; ++schritt) {
        double f_OM = funktion(OM);
        CMyVektor grad = gradient(OM, funktion);
        double grad_norm = norm(grad);

        // Test d'arrêt anticipé
        if (grad_norm < eps) {
            std::cout << "\nEnde wegen ||grad f(x)||<" << eps << " bei \n";
            std::cout << "\t x = ( ";
            for (size_t i = 0; i < OM.get_DV(); ++i) {
                std::cout << std::fixed << std::setprecision(6) << OM[i];
                if (i < OM.get_DV() - 1) std::cout << "; ";
            }
            std::cout << ")\n";
            std::cout << "\t lambda = " << lambda << "\n";
            std::cout << "\t f(x) = " << f_OM << "\n";
            std::cout << "\t grad f(x) = ( ";
            for (size_t i = 0; i < grad.get_DV(); ++i) {
                std::cout << grad[i];
                if (i < grad.get_DV() - 1) std::cout << "; ";
            }
            std::cout << ")\n";
            std::cout << "\t ||grad f(x)|| = " << grad_norm << "\n";
            break;
        }

        std::cout << "\nSchritt " << schritt << ":\n";
        std::cout << "\t x = ( ";
        for (size_t i = 0; i < OM.get_DV(); ++i) {
            std::cout << std::fixed << std::setprecision(6) << OM[i];
            if (i < OM.get_DV() - 1) std::cout << "; ";
        }
        std::cout << ")\n";
        std::cout << "\t lambda = " << lambda << "\n";
        std::cout << "\t f(x) = " << f_OM << "\n";

        std::cout << "\t grad f(x) = ( ";
        for (size_t i = 0; i < grad.get_DV(); ++i) {
            std::cout << grad[i];
            if (i < grad.get_DV() - 1) std::cout << "; ";
        }
        std::cout << ")\n";
        std::cout << "\t ||grad f(x)|| = " << grad_norm << "\n";

        //CMyVektor OM_neu = OM;
CMyVektor OM_neu = OM +grad*lambda;
       // for (size_t i = 0; i < OM.get_DV(); ++i)// kann ich den for schleife spar
          //  OM_neu[i] += lambda * grad[i];

        double f_OM_neu = funktion(OM_neu);

        std::cout << "\n\t x_neu = ( ";
        for (size_t i = 0; i < OM_neu.get_DV(); ++i) {
            std::cout << OM_neu[i];
            if (i < OM_neu.get_DV() - 1) std::cout << "; ";
        }
        std::cout << ")\n";
        std::cout << "\t f(x_neu) = " << f_OM_neu << "\n";

        if (f_OM_neu <= f_OM) {
            do {
                lambda *= 0.5;
                for (size_t i = 0; i < OM.get_DV(); ++i)
                    OM_neu[i] = OM[i] + lambda * grad[i];
                f_OM_neu = funktion(OM_neu);

                std::cout << "\n\t halbiere Schrittweite (lambda = " << lambda << "):\n";
                std::cout << "\t x_neu = ( ";
                for (size_t i = 0; i < OM_neu.get_DV(); ++i) {
                    std::cout << OM_neu[i];
                    if (i < OM_neu.get_DV() - 1) std::cout << "; ";
                }
                std::cout << ")\n";
                std::cout << "\t f(x_neu) = " << f_OM_neu << "\n";
            } while (f_OM_neu <= f_OM);

            OM = OM_neu;
        } else {
            CMyVektor OM_test = OM;
            for (size_t i = 0; i < OM.get_DV(); ++i)
                OM_test[i] += 2.0 * lambda * grad[i];

            double f_OM_test = funktion(OM_test);

            std::cout << "\n\t Test mit doppelter Schrittweite (lambda = " << 2.0 * lambda << "):\n";
            std::cout << "\t x_test = ( ";
            for (size_t i = 0; i < OM_test.get_DV(); ++i) {
                std::cout << OM_test[i];
                if (i < OM_test.get_DV() - 1) std::cout << "; ";
            }
            std::cout << ")\n";
            std::cout << "\t f(x_test) = " << f_OM_test << "\n";

            if (f_OM_test > f_OM_neu) {
                std::cout << "\t verdoppele Schrittweite!\n";
                lambda *= 2.0;
                OM = OM_test;
            } else {
                std::cout << "\t behalte alte Schrittweite!\n";
                OM = OM_neu;
            }
        }

        if (schritt == max_steps - 1) {
            std::cout << "\nEnde wegen Schrittanzahl = " << max_steps << " bei \n";
            std::cout << "\t x = ( ";
            for (size_t i = 0; i < OM.get_DV(); ++i) {
                std::cout << std::fixed << std::setprecision(6) << OM[i];
                if (i < OM.get_DV() - 1) std::cout << "; ";
            }
            std::cout << ")\n";
            std::cout << "\t lambda = " << lambda << "\n";
            std::cout << "\t f(x) = " << funktion(OM) << "\n";
            std::cout << "\t grad f(x) = ( ";
            for (size_t i = 0; i < grad.get_DV(); ++i) {
                std::cout << grad[i];
                if (i < grad.get_DV() - 1) std::cout << "; ";
            }
            std::cout << ")\n";
            std::cout << "\t ||grad f(x)|| = " << grad_norm << "\n";
            break;
        }
    }
}


// Fonction cible
double f(CMyVektor x) {
    return sin(x[0] * x[1]) + sin(x[0]) + cos(x[1]);
}

double g(CMyVektor x) {
    return -(2*x[0]*x[0] - 2*x[0]*x[1] + x[1]*x[1] + x[2]*x[2] - 2*x[0] - 4*x[2]);
}

int main() {
    
  double data1[] = {0.2, -2.1};
    CMyVektor funk1(2, data1);
    gradientenVerfahren(funk1, f);
    

    double data2[] = {0.0, 0.0, 0.0};
    CMyVektor funk2(3, data2);
    gradientenVerfahren(funk2, g, 0.1); 

    system("pause");
    
    return 0;
}
