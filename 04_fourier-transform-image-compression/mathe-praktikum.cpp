#include <iostream>
#include<vector>
#include<cmath>
#include <fstream>   
#include <vector>
#include<string>
using namespace std;
class CKomplex{
    private:
    double rteil;
    double imgteil;
    public:
    CKomplex() : rteil(0), imgteil(0) {}//Standard-Konstruktor
    CKomplex(double a, double b){rteil=a;imgteil=b;};
    CKomplex(double phi){rteil=std::cos(phi);imgteil=std::sin(phi);};
    double get_re()const{return rteil;};
    double get_im()const{return imgteil;};
    double abs()const;
};

    CKomplex operator*(const CKomplex& a, const CKomplex& b){
    double real = a.get_re() * b.get_re() - a.get_im() * b.get_im();
    double imag = a.get_re() * b.get_im() + a.get_im() * b.get_re();
    return CKomplex(real, imag);
    }
    CKomplex operator+(const CKomplex& a, const CKomplex& b){
        double real= a.get_re()+b.get_re();
        double imag= a.get_im()+b.get_im();
        return CKomplex(real,imag);
    }
    CKomplex operator*(const CKomplex & a, double b){
        double real= a.get_re()*b;
        double imag= a.get_im()*b;
        return CKomplex(real,imag);
    }
    double CKomplex:: abs()const{
        return sqrt(rteil*rteil+imgteil*imgteil);
    } 



vector<CKomplex> werte_einlesen(const std::string dateiname)
{
	int i, N, idx;
	double re, im;
	vector<CKomplex> werte;
		// File oeffnen
	ifstream fp;
	fp.open(dateiname);
		// Dimension einlesen
	fp >> N;
		// Werte-Vektor anlegen
	werte.resize(N);
	CKomplex null(0,0);
	for (i = 0; i<N; i++)
		werte[i] = null;
		// Eintraege einlesen und im Werte-Vektor ablegen
	while (!fp.eof())
	{
		fp >> idx >> re >> im;
		CKomplex a(re,im);
		werte[idx] = a;
	}
		// File schliessen
	fp.close();

	return werte;
}

void werte_ausgeben(const std::string dateiname, vector<CKomplex> werte, double epsilon=-1.0)
{
	int i;
	int N = werte.size();
		// File oeffnen
	ofstream fp;
	fp.open(dateiname);
		// Dimension in das File schreiben
	fp << N << endl;
		// Eintraege in das File schreiben
	fp.precision(10);
	for (i = 0; i < N; i++)
		if (werte[i].abs() > epsilon)
			fp << i << "\t" << werte[i].get_re() << "\t" << werte[i].get_im() << endl;
		// File schliessen
	fp.close();
}


vector<CKomplex> transformiere(const std::vector<CKomplex>& eingabe, bool rueckwaerts) {
    int N = eingabe.size();
    std::vector<CKomplex> ergebnis(N);

    const double PI_KONSTANTE = 3.141592653589793;
    double norm_faktor = 1.0 / std::sqrt(static_cast<double>(N));

    // Erzeuge jeden Zielwert im Ergebnisvektor
    for (int ziel = 0; ziel < N; ++ziel) {
        CKomplex summe(0.0, 0.0);  

        //Summiere über den Eingangsvektor
        for (int quelle = 0; quelle < N; ++quelle) {
            double winkel = (2.0 * PI_KONSTANTE * quelle * ziel) / N; 

            if (!rueckwaerts) {
                winkel = -winkel;
            }

            double realanteil = std::cos(winkel);
            double imaganteil = std::sin(winkel);
            CKomplex exponential(realanteil, imaganteil);

            
            summe = summe + (eingabe[quelle] * exponential);
        }

        ergebnis[ziel] = summe * norm_faktor;
    }

    return ergebnis;
}



// Funktion zur Berechnung der maximalen Abweichung
double maximale_abweichung(const vector<CKomplex>& original, const vector<CKomplex>& rekonstruiert) {
    double max_diff = 0.0;
    for (size_t i = 0; i < original.size(); ++i) {
        double differenz = (original[i] + (rekonstruiert[i] * -1.0)).abs(); // |original - rekonstruiert|
        if (differenz > max_diff)
            max_diff = differenz;
    }
    return max_diff;
}

// Diese Funktion testet einen Datensatz mit mehreren epsilons
void teste_datei(const string& dateiname, const string& label) {
    vector<CKomplex> original = werte_einlesen(dateiname);
    vector<CKomplex> spektrum = transformiere(original, false);//Fourier-Transformation anwenden → ergibt Frequenzspektrum

    cout << "\nBei " << label << ":\n";

    // Standard-Epsilon-Fall (default = -1.0)
    string std_filename = "temp_default.txt";
    werte_ausgeben(std_filename, spektrum, -1.0);
    vector<CKomplex> geladen_std = werte_einlesen(std_filename);
    vector<CKomplex> rueck_std = transformiere(geladen_std, true);
    double fehler_std = maximale_abweichung(original, rueck_std);
    cout << "Maximale Abweichung bei Standard-Epsilon: ca. " << scientific << fehler_std << "\n";

    // Weitere epsilon-Werte testen
    vector<double> epsilons = {0.001, 0.01, 0.1, 1.0};
    for (double eps : epsilons) {
        string fname = "temp_eps_" + to_string(eps) + ".txt";
        werte_ausgeben(fname, spektrum, eps);
        vector<CKomplex> geladen = werte_einlesen(fname);
        vector<CKomplex> rueck = transformiere(geladen, true);
        double fehler = maximale_abweichung(original, rueck);
        cout << "Maximale Abweichung bei epsilon=" << eps << ": " << fehler << "\n";
    }
}//rein → transformieren → filtern → rückrechnen → vergleichen.

void verarbeite_roi_mit_epsilons() {
    string basis = "roi-mohammedVI";
    string eingabe_datei = basis + ".txt";

    vector<CKomplex> original = werte_einlesen(eingabe_datei);
    vector<CKomplex> spektrum = transformiere(original, false); // Fourier hin

    
    werte_ausgeben(basis + "_standard_hin.txt", spektrum, -1.0);

    // Definierte epsilons (wie in der Aufgabenstellung)
    vector<double> epsilons = {-1,10, 30, 100, 300, 1000};

    for (double eps : epsilons) {
        // Spektrum filtern und speichern
        string spektrum_datei = basis + "_hin_eps_" + to_string((int)eps) + ".txt";
        werte_ausgeben(spektrum_datei, spektrum, eps);

        // Rücktransformation
        vector<CKomplex> geladen = werte_einlesen(spektrum_datei);
        vector<CKomplex> rueck = transformiere(geladen, true);

        // Rücktransformierte Werte speichern
        string rueck_datei = "rueck_" + basis + "_" + to_string((int)eps) + ".txt";
        werte_ausgeben(rueck_datei, rueck);

        //double fehler = maximale_abweichung(original, rueck);
        //cout << "Epsilon " << eps << " → maximale Abweichung: " << fehler << endl;
    }

    cout << "Alle Rücktransformationen abgeschlossen. Nun mit image2array.exe zurück in Bilder umwandeln." << endl;
}






int main() {
    teste_datei("Daten_original1.txt", "Daten_original1.txt");
    teste_datei("Daten_original2.txt", "Daten_original2.txt");
    verarbeite_roi_mit_epsilons();
    system("pause");
    return 0;
}
