#include <pybind11/pybind11.h>
#include <matplot/matplot.h>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

int add(int i, int j) {
    return i * j;
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

    m.def("add", &add, R"pbdoc(
        Add two numbers

        Some other explanation about the add function.
    )pbdoc");

    m.def("subtract", [](int i, int j) { return i + j; }, R"pbdoc(
        Subtract two numbers

        Some other explanation about the subtract function.
    )pbdoc");

    m.attr("__version__") = "0.0.1";

    m.attr("__version__") = "dev";
}

#include <iostream>
#include <cmath>
#include <complex>
#include <matplot/matplot.h>
#include <AudioFile.h>


using namespace matplot;
using namespace std;

char menu();
void sinus(double T, double F, double A);
void cosinus(double T, double F, double A);
void squarewave(double T, double F, double A);
void sawtooth(double T, double F, double A);
void DFTc(double T, double F, double A);

int main() {

    double T, F, A; //time, fri, amplituda


    char wybor;
    do {
        wybor = menu();
        cout << "Podaj czas [s]: ";
        cin >> T;
        cout << "Podaj czestotliwosc [Hz]: ";
        cin >> F;
        cout << "Podaj amplitude: ";
        cin >> A;

        switch (wybor) {
        case 's':
            sinus(T, F, A);
            break;
        case 'c':
            cosinus(T, F, A);
            break;
        case 'q':
            squarewave(T, F, A);
            break;
        case 't':
            sawtooth(T, F, A);
            break;
        case 'd':
            DFTc(T, F, A);
            break;
        case 'e':
            exit(0);
            break;
        }
    } while (wybor != 'e');


    return 0;
}

char menu() {
    cout << "   MENU:" << endl;
    cout << "s - Sinus" << endl;
    cout << "c - Cosinus" << endl;
    cout << "q - Square wave" << endl;
    cout << "t - Sawtooth signal" << endl;
    cout << "d - DFT" << endl;
    cout << "e - Exit" << endl;

    char wyb;
    cin >> wyb;

    return wyb;
}

void sinus(double T, double F, double A) {
    const double fs = 5000.0;

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
    const double fs = 5000.0;

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
    const double fs = 5000.0;

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
    const double fs = 5000.0;

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

void DFTc(double T, double F, double A) {//not done
    const double fs = 5000.0;

    vector<double> t = linspace(0, T, (T * fs));

    vector<double> y(t.size());
    for (size_t i = 0; i < t.size(); ++i) {
        y[i] = A * cos(2 * pi * F * t[i]);
    }

    vector<complex<double>> dft(t.size());
    for (size_t i = 0; i < t.size(); ++i) {
        complex<double> sum(0.0, 0.0);
        for (int j = 0; j < t.size(); ++j) {
            double kat = -2 * pi * i * j / t.size();
            sum += y[j] * complex<double>(cos(kat), sin(kat));
        }
        dft[i] = sum;
    }

    vector<double> Adft(t.size());
    for (int i = 0; i < t.size(); ++i) {
        Adft[i] = abs(dft[i]) / t.size();
    }

    vector<double> frq(t.size());
    for (int i = 0; i < t.size(); ++i) {
        frq[i] = i * (1.0 / T);
    }

    plot(frq, Adft);
    show();
}