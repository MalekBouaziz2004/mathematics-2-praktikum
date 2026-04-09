#include <iostream>     
#include <cstdlib>      
#include <ctime>
#include <time.h>
#include<vector>


//1. Aufgabe

class CZufall{
    public:
    int wert(int a, int b){int zufall = rand() % (b - a + 1) + a;
    return zufall;};
    void initialisiere(int s){srand(s);};
    void test(int a, int b, int N){
        int array[b-a+1]={0};
        for(int i=0;i<N;++i){
            int zufall=wert(a,b);
            int index=zufall-a;
            array[index]++;
            

        }
         for (int i = 0; i <= (b - a); ++i) {
        std::cout << (a + i) << " kam " << array[i] << " mal vor.\n";
    }
    };
    void test_falsch(int a, int b, int N){  int array[b-a+1]={0};
        for(int i=0;i<N;++i){
            initialisiere(time(NULL));// Dadurch soll jede Ziehung theoretisch eine neue Zufallsfolge starten.
            int zufall=wert(a,b);
            int index=zufall-a;
            array[index]++;
            

        }
         for (int i = 0; i <= (b - a); ++i) {
        std::cout << (a + i) << " kam " << array[i] << " mal vor.\n";
    }
    };

};

//2. Aufgabe: Poker

struct Karte {
    int wert;
    int farbe;
};

std::vector<Karte> erzeugeKartendeck() {
    std::vector<Karte> deck;
    deck.reserve(52); 

    for (int w = 0; w < 13; ++w) {       
        for (int f = 0; f < 4; ++f) {     
            Karte k;
            k.wert = w;
            k.farbe = f;
            deck.push_back(k);
        }
    }

    return deck;
}

void mischeDeck(std::vector<Karte>& deck) {
    int n = deck.size();
    for (int i = n - 1; i > 0; --i) {
        int j = rand() % (i + 1);
        std::swap(deck[i], deck[j]); 
    }
}




void simulierePoker(int versuche, int& paareOut, int& drillingeOut) {
    paareOut = 0;
    drillingeOut = 0;

    std::vector<Karte> basisDeck = erzeugeKartendeck();

    for (int lauf = 0; lauf < versuche; ++lauf) {
        // Neues gemischtes Deck erzeugen
        std::vector<Karte> deck = basisDeck;
        mischeDeck(deck);

        // Zähler zurücksetzen
        int zaehler[13] = {0};

        // Sieben Karten ziehen und Werte zählen
        for (int i = 0; i < 7; ++i) {
            ++zaehler[ deck[i].wert ];
        }

        // Merker, ob Paar oder Drilling gefunden wurde
        int paarAnzahl = 0;
bool hatDrilling = false;

for (int i = 0; i < 13; ++i) {
    paarAnzahl += zaehler[i] / 2;  // 2 Karten = 1 Paar, 4 Karten = 2 Paare
    if (zaehler[i] >= 3)
        hatDrilling = true;
}

bool hatMindZweiPaare = (paarAnzahl >= 2);

if (hatMindZweiPaare) ++paareOut;   // jetzt wirklich ≥ 2 Paare
if (hatDrilling)       ++drillingeOut;

    }
}











int main(){
    CZufall member;
    member.initialisiere(10);
    member.test(3,7,10000);
    std::cout<<"-----------------------------------"<<std::endl;
    member.initialisiere(10);
    member.test(3,7,10000);//Der Zufallszahlengenerator in C++ ist deterministisch.
 //Das heißt: Wenn man ihn mit dem gleichen Startwert (Seed) initialisiert, erzeugt er immer genau dieselbe Folge von Zufallszahlen.


//----------------------------------------------

    member.initialisiere(15);
    member.test(3,7,10000);
    std::cout<<"-----------------------------------"<<std::endl;
    member.initialisiere(12);
    member.test(3,7,10000);
//--------------------------------------------------


    member.initialisiere(time(NULL));
    member.test(3,7,10000);
    std::cout<<"-----------------------------------"<<std::endl;
    member.initialisiere(time(NULL));
    member.test(3,7,10000);//time(NULL) gibt die aktuelle Zeit in Sekunden,Wenn  Programm schnell genug läuft liefert gleiches Ausgabe*/

//---------------------------------------------------

member.test_falsch(3,7,10000);

//----------------------------------------------

srand(static_cast<unsigned>(time(nullptr)));   // einmal seeden

    const int N = 100'000;
    int paare = 0, drillinge = 0;

    simulierePoker(N, paare, drillinge);          

    std::cout << "min. ein Paar:     "
              << static_cast<double>(paare) / N * 100 << " %\n";
    std::cout << "min. ein Drilling: "
              << static_cast<double>(drillinge) / N * 100 << " %\n";

    

    system("pause");
    return 0;
}