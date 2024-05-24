#include <cmath>
#include <complex>
#include <pybind11/pybind11.h>
#include <matplot/matplot.h>
#include <AudioFile.h>
#include <vector>

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

const double fs = 4096.0;

using namespace matplot;
using namespace std;


void display_sinus(double T, double F, double A) {

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

void display_cosinus(double T, double F, double A) {

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

void display_squarewave(double T, double F, double A) {

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

void display_sawtooth(double T, double F, double A) {

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

void display_audio() {
    AudioFile<double> audioFile;
    string filePath = "taco1.wav";  // Taco Hemingway

    if (!audioFile.load(filePath)) {
        cerr << "Nie można wczytać pliku audio: " << filePath << endl;
    }

    double T = audioFile.getLengthInSeconds();

    vector<double> t;
    t = linspace(0, T, (T * fs));

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
}

void DFT_IDFT(char wybor, double F) {

    vector<double> t = linspace(0, 1, fs);
    vector<double> y(t.size()); 
    switch (wybor) {
    case 's':
        for (size_t i = 0; i < t.size(); ++i) {
            y[i] = sin(2 * pi * F * t[i]);
        }
        break;
    case 'c':
        for (size_t i = 0; i < t.size(); ++i) {
            y[i] = cos(2 * pi * F * t[i]);
        }
        break;
    case 'q':
        for (size_t i = 0; i < t.size(); ++i) {
            y[i] = sin(2 * pi * F * t[i]) > 0 ? 1 : -1;
        }
        break;
    case 't':
        for (size_t i = 0; i < t.size(); ++i) {
            y[i] = fmod(F * t[i] + 0.5, 1) * 2 - 1;
        }
        break;
    default: 
        exit(0);
    }

    double n = y.size();
    vector<complex<double>> dft(n);
    for (size_t i = 0; i < n; ++i) {
        complex<double> sum(0.0, 0.0);
        for (size_t j = 0; j < n; ++j) {
            double kat = 2 * pi * i * j / n;
            sum += y[j] * polar(1.0, -kat); //zespolak 
        }
        dft[i] = sum;
    }

    vector<double> Adft(n);
    for (int i = 0; i < n; ++i) {
        Adft[i] = abs(dft[i]);
    }

    vector<double> frq(n);
    for (int i = 0; i < n; ++i) {
        frq[i] = i * (fs / n);
    }

    
    //double n = y.size();
    vector<double> Idft(n);
    for (size_t i = 0; i < n; ++i) {
        complex<double> sum(0.0, 0.0);
        for (size_t j = 0; j < n; ++j) {
            double kat = 2 * pi * i * j / n;
            sum += dft[j] * polar(1.0, kat);
        }
        Idft[i] = sum.real() / fs;
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

void usuwanie(char wybor, double F, double uf) {

    vector<double> t = linspace(0, 1, fs);
    vector<double> y(t.size());
    AudioFile<double> audioFile;
    string filePath = "taco1.wav";  // Taco Hemingway


    switch (wybor) {
    case 's':
        for (size_t i = 0; i < t.size(); ++i) {
            y[i] = sin(2 * pi * F * t[i]);
        }
        break;
    case 'c':
        for (size_t i = 0; i < t.size(); ++i) {
            y[i] = cos(2 * pi * F * t[i]);
        }
        break;
    case 'q':
        for (size_t i = 0; i < t.size(); ++i) {
            y[i] = sin(2 * pi * F * t[i]) > 0 ? 1 : -1;
        }
        break;
    case 't':
        for (size_t i = 0; i < t.size(); ++i) {
            y[i] = fmod(F * t[i] + 0.5, 1) * 2 - 1;
        }
        break;
    case 'a':
        if (!audioFile.load(filePath)) {
            cerr << "Nie można wczytać pliku audio: " << filePath << endl;
        }
        t = linspace(0, audioFile.getLengthInSeconds(), (audioFile.getLengthInSeconds() * fs));
        y.reserve(fs);
        for (int i = 0; i < fs; ++i) {
            y.push_back(audioFile.samples[0][i]); // [0] - tylko jeden kanał 
        }
        break;
    default:
        exit(0);
    }
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

    for (int i = 0; i < syllabus; ++i) {
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

void Audio_DFT_IDFT() {
    AudioFile<double> audioFile;
    string filePath = "taco1.wav";  // Taco Hemingway
    if (!audioFile.load(filePath)) {
        cerr << "Nie można wczytać pliku audio: " << filePath << endl;
    }
    double T = audioFile.getLengthInSeconds(); // czas trwania audio

    vector<double> t;
    t = linspace(0, T, (T * fs));

    vector<double> y;
    y.reserve(fs);
    for (int i = 0; i < fs; ++i) {
        y.push_back(audioFile.samples[0][i]); // [0] - tylko jeden kanał 
    }

    double n = y.size();
    vector<complex<double>> dft(n);
    for (size_t i = 0; i < n; ++i) {
        complex<double> sum(0.0, 0.0);
        for (size_t j = 0; j < n; ++j) {
            double kat = 2 * pi * i * j / n;
            sum += y[j] * polar(1.0, -kat); //zespolak 
        }
        dft[i] = sum;
    }

    vector<double> Adft(n);
    for (int i = 0; i < n; ++i) {
        Adft[i] = abs(dft[i]);
    }

    vector<double> frq(n);
    for (int i = 0; i < n; ++i) {
        frq[i] = i * (fs / n);
    }

    vector<double> Idft(n);
    for (size_t i = 0; i < n; ++i) {
        complex<double> sum(0.0, 0.0);
        for (size_t j = 0; j < n; ++j) {
            double kat = 2 * pi * i * j / n;
            sum += dft[j] * polar(1.0, kat);
        }
        Idft[i] = sum.real() / fs;
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

void Audio_usuwanie(double uf) {
    AudioFile<double> audioFile;
    string filePath = "taco1.wav";  // Taco Hemingway
    if (!audioFile.load(filePath)) {
        cerr << "Nie można wczytać pliku audio: " << filePath << endl;
    }
    double T = audioFile.getLengthInSeconds(); // czas trwania audio

    vector<double> t;
    t = linspace(0, T, (T * fs));

    vector<double> y;
    y.reserve(fs);
    for (int i = 0; i < fs; ++i) {
        y.push_back(audioFile.samples[0][i]); // [0] - tylko jeden kanał 
    }

    double n = y.size();
    vector<complex<double>> dft(n);
    for (size_t i = 0; i < n; ++i) {
        complex<double> sum(0.0, 0.0);
        for (size_t j = 0; j < n; ++j) {
            double kat = 2 * pi * i * j / n;
            sum += y[j] * polar(1.0, -kat); //zespolak 
        }
        dft[i] = sum;
    }

    size_t syllabus = static_cast<size_t>(uf * dft.size() / fs);

    for (int i = 0; i < syllabus; ++i) {
        dft[i] = 0.0;
    }

    vector<double> cutted(dft.size());
    for (size_t i = 0; i < n; ++i) {
        complex<double> sum(0.0, 0.0);
        for (size_t j = 0; j < n; ++j) {
            double kat = 2 * pi * i * j / n;
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


namespace py = pybind11;

PYBIND11_MODULE(_core, m) {
    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: scikit_build_example

        .. autosummary::
           :toctree: _generate

    "s - Sinus"
    "c - Cosinus"
    "q - Square wave"
    "t - Sawtooth signal"


    )pbdoc";

    m.def("display_sinus", &display_sinus, R"pbdoc(
       Wyswietlanie sinusa, argumenty: czas, czestotliowsc, amplituda
    )pbdoc");

    m.def("display_cosinus", &display_cosinus, R"pbdoc(
       Wyswietlanie cosinusa, argumenty: czas, czestotliowsc, amplituda
    )pbdoc");

    m.def("display_squarewave", &display_squarewave, R"pbdoc(
       Wyswietlanie fali prostokatnej, argumenty: czas, czestotliowsc, amplituda
    )pbdoc");

    m.def("display_sawtooth", &display_sawtooth, R"pbdoc(
       Wyswietlanie fali trojkatnej, argumenty: czas, czestotliowsc, amplituda
    )pbdoc");

    m.def("display_audio", &display_audio, R"pbdoc(
       Wyswietlanie audio
    )pbdoc");

    m.def("DFT_IDFT", &DFT_IDFT, R"pbdoc(
       Operacja DFT oraz IDFT, argumenty: char, czestotliwosc
    )pbdoc");

    m.def("usuwanie", &usuwanie, R"pbdoc(
       Usuwanie niskich częstotliwosci za pomoca DFT
    )pbdoc");

    m.def("Audio_DFT_IDFT", &Audio_DFT_IDFT, R"pbdoc(
       Operacja DFT oraz IDFT na pliku dzwiekowym
    )pbdoc");

    m.def("Audio_usuwanie", &Audio_usuwanie, R"pbdoc(
       Usuwanie niskich częstotliwosci za pomoca DFT z pliku dzwiekowego
    )pbdoc");
    

    m.attr("__version__") = "0.0.1";

    m.attr("__version__") = "dev";
}
