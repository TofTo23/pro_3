#include <cmath>
#include <complex>
#include <pybind11/pybind11.h>
#include <matplot/matplot.h>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

const double fs = 4096.0;

using namespace matplot;
using namespace std;


void sinus(double T, double F, double A) {

    vector<double> t = linspace(0, T, (T * fs));

    vector<double> y(t.size());
    for (size_t i = 0; i < t.size(); ++i) {
        y[i] = A * sin(2 * pi * F * t[i]);
    }

    plot(t, y);
    title("Sinus");
    xlabel("t(s)");
    ylabel("y");
    show();
}

void cosinus(double T, double F, double A) {

    vector<double> t = linspace(0, T, (T * fs));

    vector<double> y(t.size());
    for (size_t i = 0; i < t.size(); ++i) {
        y[i] = A * cos(2 * pi * F * t[i]);
    }

    plot(t, y);
    title("Cosinus");
    xlabel("t(s)");
    ylabel("y");
    show();
}

void squarewave(double T, double F, double A) {

    vector<double> t = linspace(0, T, (T * fs));

    vector<double> y(t.size());
    for (size_t i = 0; i < t.size(); ++i) {
        y[i] = A * sin(2 * pi * F * t[i]) > 0 ? A : -A;
    }

    plot(t, y);
    title("Square wave");
    xlabel("t(s)");
    ylabel("y");
    show();
}

void sawtooth(double T, double F, double A) {

    vector<double> t = linspace(0, T, (T * fs));

    vector<double> y(t.size());
    for (size_t i = 0; i < t.size(); ++i) {
        y[i] = fmod(F * t[i] + 0.5, 1) * (2 * A) - A;
    }

    plot(t, y);
    title("Sawtooth wave");
    xlabel("t(s)");
    ylabel("y");
    show();
}

namespace py = pybind11;

PYBIND11_MODULE(_core, m) {
    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: scikit_build_example

        .. autosummary::
           :toctree: _generate

           add
           subtract
    )pbdoc";

    m.def("sinus", &sinus, R"pbdoc(
       Wyswietlanie sinusa, argumenty: 
    )pbdoc");

    m.def("cosinus", &cosinus, R"pbdoc(
       display sinus
    )pbdoc");

    m.def("squarewave", &squarewave, R"pbdoc(
       display sinus
    )pbdoc");

    m.def("sawtooth", &sawtooth, R"pbdoc(
       display sinus
    )pbdoc");

    m.def("subtract", [](int i, int j) { return i + j; }, R"pbdoc(
        Subtract two numbers

        Some other explanation about the subtract function.
    )pbdoc");

    m.attr("__version__") = "0.0.1";

    m.attr("__version__") = "dev";
}
//--------------------------------------------------------------------------------------
#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <matplot/matplot.h>
#include <AudioFile.h>

using namespace matplot;
using namespace std;

const double fs = 4096.0;  

char menu();
vector<double> audio(double* T);
vector<double> sinus(vector<double> t, double F, double A);
vector<double> cosinus(vector<double> t, double F, double A);
vector<double> squarewave(vector<double> t, double F, double A);
vector<double> sawtooth(vector<double> t, double F, double A);
void DFT_IDFT(vector<double> t, vector<double> y);
void usuwanie(vector<double> t, vector<double> y);

int main() {

    double T = 0, F, A; //time, fri, amplituda
    vector<double> y;
    vector<double> t;

    char wybor;
    do {
        wybor = menu();
        if (wybor != 'a' && wybor != 'e' && wybor != 'd' && wybor != 'u') {
            cout << "Podaj czas [s]: ";
            cin >> T;
            cout << "Podaj czestotliwosc [Hz]: ";
            cin >> F;
            cout << "Podaj amplitude: ";
            cin >> A;
            t = linspace(0, T, (T * fs));
        }
        

        switch(wybor) {
            case 'a':
                y = audio(&T);
                t = linspace(0, T, (T * fs));
                break;
            case 's':
                y = sinus(t, F, A);
                break;
            case 'c':
                y = cosinus(t, F, A);
                break;
            case 'q':
                y = squarewave(t, F, A);
                break;
            case 't':
                y = sawtooth(t, F, A);
                break;
            case 'd':
                DFT_IDFT(t, y);
                break;
            case 'u':
                usuwanie(t, y);
                break;
            case 'e':
                exit(0);
                break;
        }
    } while (wybor != 'e');
    

    return 0;
}

char menu() {
    cout << "    MENU:" << endl;
    cout << "s - Sinus" << endl;
    cout << "c - Cosinus" << endl;
    cout << "q - Square wave" << endl;
    cout << "t - Sawtooth signal" << endl;
    cout << "a - Audio" << endl;
    cout << "d - DFT & IDFT" << endl;
    cout << "u - Usuwanie" << endl;
    cout << "e - Exit" << endl;

    char wyb;
    cin >> wyb;

    return wyb;
}

vector<double> sinus(vector<double> t, double F, double A) {

    vector<double> y(t.size());
    for (size_t i = 0; i < t.size(); ++i) {
        y[i] = A * sin(2 * pi * F * t[i]);
    }

    auto w1 = figure();
    w1->size(800, 600);
    plot(t, y);
    title("Sinus");
    xlabel("t(s)");
    ylabel("y");
    show();

    return y;
}

vector<double> cosinus(vector<double> t, double F, double A) {

    vector<double> y(t.size());
    for (size_t i = 0; i < t.size(); ++i) {
        y[i] = A * cos(2 * pi * F * t[i]);
    }

    auto w1 = figure();
    w1->size(800, 600);
    plot(t, y);
    title("Cosinus");
    xlabel("t(s)");
    ylabel("y");
    show();

    return y;
}

vector<double> squarewave(vector<double> t, double F, double A) {

    vector<double> y(t.size());
    for (size_t i = 0; i < t.size(); ++i) {
        y[i] = A * sin(2 * pi * F * t[i]) > 0 ? A : - A;
    }

    auto w1 = figure();
    w1->size(800, 600);
    plot(t, y);
    title("Square wave");
    xlabel("t(s)");
    ylabel("y");
    show();

    return y;
}

vector<double> sawtooth(vector<double> t, double F, double A) {

    vector<double> y(t.size());
    for (size_t i = 0; i < t.size(); ++i) {
        y[i] = fmod(F * t[i] + 0.5 , 1) * (2 * A) - A;
    }

    auto w1 = figure();
    w1->size(800, 600);
    plot(t, y);
    title("Sawtooth wave");
    xlabel("t(s)");
    ylabel("y");
    show();

    return y;
}


vector<double> audio(double* T) {
    AudioFile<double> audioFile;
    string filePath = "/Users/Lenovo/Documents/000Moje/XPP/PRO_3/taco1.wav";  // Taco Hem
    //string filePath = "/Users/Lenovo/Documents/000Moje/XPP/PRO_3/rising-frq.wav";  // plik o zmieniajacej sie (rosnacej) frq

    if (!audioFile.load(filePath)) {
        cerr << "Nie można wczytać pliku audio: " << filePath << endl;
    }

    // info about T, couse gotta creat vector t
    *T = audioFile.getLengthInSeconds();

    vector<double> t;
    t = linspace(0, *T, (*T * fs));

    vector<double> USamples;
    USamples.reserve(fs);

    for (int i = 0; i < fs; ++i) {
        USamples.push_back(audioFile.samples[0][i]); // [0] - tylko jeden kanał 
    }

    // --- wykres dla fs sampli
    auto w1 = figure();
    w1->size(800, 600);
    plot(t, USamples);
    title("Audio");
    xlabel("Próbki");
    ylabel("Amplituda");
    show();

    return USamples;
}

void DFT_IDFT(vector<double> t, vector<double> y) {
    
    vector<complex<double>> dft(y.size());
    for (size_t i = 0; i < y.size(); ++i) {
        complex<double> sum(0.0, 0.0);
        for (size_t j = 0; j < y.size(); ++j) {
            double kat = 2 * pi * i * j / y.size();
            sum += y[j] * polar(1.0, -kat); //zespolak 
        }
        dft[i] = sum;
    }

    vector<double> Adft(y.size());
    for(int i = 0; i < y.size(); ++i) {
        Adft[i] = abs(dft[i]) ;
    }

    vector<double> frq(y.size());
    for (int i = 0; i < y.size(); ++i) {
        frq[i] = i * (fs / y.size());
    }

    double n = y.size();
    vector<double> Idft(y.size());
    for (size_t i = 0; i < y.size(); ++i) {
        complex<double> sum(0.0, 0.0);
        for (size_t j = 0; j < y.size(); ++j) {
            double kat = 2 * pi * i * j / y.size();
            sum += dft[j] * polar(1.0, kat);
        }
        Idft[i] = sum.real() / n;
    }
    
    auto w2 = figure();
    w2->size(800, 600);
    plot(frq, Adft);
    title("DFT");
    xlabel("Częstotliwość");
    ylabel("Amplituda");

    auto w3 = figure();
    w3->size(800, 600);
    plot(t, Idft);
    title("IDFT");
    xlabel("Czas");
    ylabel("Amplituda");
    show();
}

void usuwanie(vector<double> t, vector<double> y) {
    double uf = 0;
    cout << "Podaj czestotlowosc do usuniecia: ";
    cin >> uf;

    vector<complex<double>> dft(y.size());
    for (size_t i = 0; i < y.size(); ++i) {
        complex<double> sum(0.0, 0.0);
        for (size_t j = 0; j < y.size(); ++j) {
            double kat = 2 * pi * i * j / y.size();
            sum += y[j] * polar(1.0, -kat); 
        }
        dft[i] = sum;
    }

    size_t syllabus = static_cast<size_t>(uf * dft.size() / fs);

    for (int i = 0; i < syllabus ; ++i) {
        dft[i] = 0.0;
    }

    vector<double> cutted(dft.size());
    double n = y.size();
    for (size_t i = 0; i < y.size(); ++i) {
        complex<double> sum(0.0, 0.0);
        for (size_t j = 0; j < y.size(); ++j) {
            double kat = 2 * pi * i * j / y.size();
            sum += dft[j] * polar(1.0, kat);
        }
        cutted[i] = sum.real() / n;
    }

    auto w4 = figure();
    w4->size(800, 600);
    plot(t, y);
    title("Oryginalny sygnał");
    xlabel("Czas");
    ylabel("Amplituda");
    
    auto w5 = figure();
    w5->size(800, 600);
    plot(t, cutted);
    title("Sygnał po usunięciu niskich częstotliwości");
    xlabel("Czas");
    ylabel("Amplituda");
    
    show();
}




