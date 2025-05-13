#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>

using namespace std;

struct Harmonic {
    double amplitude;
    double frequency;
    double phase;
};

// Функция для генерации сигнала
vector<double> generateSignal(const vector<Harmonic>& harmonics,
    double samplingRate,
    double duration,
    int& numPoints) {
    double deltaT = 1.0 / samplingRate;
    numPoints = static_cast<int>(duration * samplingRate);
    vector<double> signal(numPoints, 0.0);

    for (int i = 0; i < numPoints; ++i) {
        double t = i * deltaT;
        for (const auto& h : harmonics) {
            signal[i] += h.amplitude * sin(2 * M_PI * h.frequency * t + h.phase);
        }
    }

    return signal;
}

// Функция для квантования сигнала
vector<int> quantizeSignal(const vector<double>& signal, int levels, double maxAmplitude) {
    vector<int> quantized(signal.size(), 0);
    double step = (2 * maxAmplitude) / levels;

    for (size_t i = 0; i < signal.size(); ++i) {
        quantized[i] = static_cast<int>((signal[i] + maxAmplitude) / step);
    }

    return quantized;
}

// Функция для записи данных в файл
void writeToFile(const string& filename,
    const vector<double>& time,
    const vector<double>& signal,
    const vector<int>& quantized = vector<int>()) {
    ofstream out(filename);
    if (!out.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    out << "Time,Signal";
    if (!quantized.empty()) {
        out << ",Quantized";
    }
    out << "\n";

    for (size_t i = 0; i < time.size(); ++i) {
        out << time[i] << "," << signal[i];
        if (!quantized.empty() && i < quantized.size()) {
            out << "," << quantized[i];
        }
        out << "\n";
    }

    out.close();
}

// Функция для чтения гармоник из файла
vector<Harmonic> readHarmonicsFromFile(const string& filename) {
    vector<Harmonic> harmonics;
    ifstream in(filename);
    if (!in.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return harmonics;
    }

    Harmonic h;
    while (in >> h.amplitude >> h.frequency >> h.phase) {
        harmonics.push_back(h);
    }

    in.close();
    return harmonics;
}

int main() {
    vector<Harmonic> harmonics;
    char choice;

    cout << "Read harmonics from file? (y/n): ";
    cin >> choice;

    if (choice == 'y' || choice == 'Y') {
        string filename;
        cout << "Enter filename: ";
        cin >> filename;
        harmonics = readHarmonicsFromFile(filename);
    }
    else {
        int numHarmonics;
        cout << "Enter number of harmonics: ";
        cin >> numHarmonics;

        harmonics.resize(numHarmonics);
        for (int i = 0; i < numHarmonics; ++i) {
            cout << "Harmonic " << i + 1 << ":\n";
            cout << "  Amplitude: ";
            cin >> harmonics[i].amplitude;
            cout << "  Frequency (Hz): ";
            cin >> harmonics[i].frequency;
            cout << "  Phase (rad): ";
            cin >> harmonics[i].phase;
        }
    }

    double samplingRate, duration;
    cout << "Enter sampling rate (Hz): ";
    cin >> samplingRate;
    cout << "Enter duration (s): ";
    cin >> duration;

    int numPoints;
    auto signal = generateSignal(harmonics, samplingRate, duration, numPoints);

    // Генерируем временную ось
    vector<double> time(numPoints);
    double deltaT = 1.0 / samplingRate;
    for (int i = 0; i < numPoints; ++i) {
        time[i] = i * deltaT;
    }

    // Находим максимальную амплитуду для квантования
    double maxAmplitude = 0;
    for (const auto& h : harmonics) {
        maxAmplitude += h.amplitude;
    }

    // Квантование сигнала
    int quantizationLevels;
    cout << "Enter number of quantization levels: ";
    cin >> quantizationLevels;

    auto quantized = quantizeSignal(signal, quantizationLevels, maxAmplitude);

    // Запись результатов в файл
    writeToFile("signal_data.csv", time, signal, quantized);
    cout << "Data written to signal_data.csv\n";

    return 0;
}