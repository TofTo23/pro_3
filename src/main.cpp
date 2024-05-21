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


