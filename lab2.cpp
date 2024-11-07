#include <iostream>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <cmath>

using namespace std;

// Абстрактный класс для матриц
template <typename T>
class Matrix {
public:
    virtual unsigned getRows() const = 0;
    virtual unsigned getCols() const = 0;
    virtual T& at(unsigned row, unsigned col) = 0;
    virtual const T& at(unsigned row, unsigned col) const = 0;
    virtual void print() const = 0;
    virtual void importFromFile(const string& filename) = 0;
    virtual void exportToFile(const string& filename) const = 0;
    virtual ~Matrix() = default;
};

// Класс плотной матрицы
template <typename T = double>
class MatrixDense : public Matrix<T> {
private:
    unsigned _m, _n;
    T* data;

public:
    MatrixDense(unsigned m, unsigned n) : _m(m), _n(n) {
        data = new T[m * n]();
    }

    ~MatrixDense() {
        delete[] data;
    }

    unsigned getRows() const override {
        return _m;
    }

    unsigned getCols() const override {
        return _n;
    }

    T& at(unsigned row, unsigned col) override {
        if (row >= _m || col >= _n) throw out_of_range("Index out of range");
        return data[row * _n + col];
    }

    const T& at(unsigned row, unsigned col) const override {
        if (row >= _m || col >= _n) throw out_of_range("Index out of range");
        return data[row * _n + col];
    }

    void print() const override {
        for (unsigned i = 0; i < _m; ++i) {
            for (unsigned j = 0; j < _n; ++j) {
                cout << at(i, j) << " ";
            }
            cout << endl;
        }
    }

    void importFromFile(const string& filename) override {
        ifstream file(filename);
        if (!file) throw runtime_error("Failed to open file");

        string className;
        file >> className;
        if (className != "MatrixDense") throw runtime_error("Incorrect matrix type");

        file >> _m >> _n;
        delete[] data;
        data = new T[_m * _n];

        for (unsigned i = 0; i < _m; ++i) {
            for (unsigned j = 0; j < _n; ++j) {
                file >> at(i, j);
            }
        }
    }

    void exportToFile(const string& filename) const override {
        ofstream file(filename);
        if (!file) throw runtime_error("Failed to open file");

        file << "MatrixDense\n";
        file << _m << " " << _n << endl;
        for (unsigned i = 0; i < _m; ++i) {
            for (unsigned j = 0; j < _n; ++j) {
                file << at(i, j) << " ";
            }
            file << endl;
        }
    }

    // Операции над матрицами
    MatrixDense<T> operator+(const MatrixDense<T>& other) const {
        if (_m != other._m || _n != other._n) throw runtime_error("Matrix dimensions must match");
        MatrixDense<T> result(_m, _n);
        for (unsigned i = 0; i < _m; ++i) {
            for (unsigned j = 0; j < _n; ++j) {
                result.at(i, j) = at(i, j) + other.at(i, j);
            }
        }
        return result;
    }

    MatrixDense<T> operator-(const MatrixDense<T>& other) const {
        if (_m != other._m || _n != other._n) throw runtime_error("Matrix dimensions must match");
        MatrixDense<T> result(_m, _n);
        for (unsigned i = 0; i < _m; ++i) {
            for (unsigned j = 0; j < _n; ++j) {
                result.at(i, j) = at(i, j) - other.at(i, j);
            }
        }
        return result;
    }

    MatrixDense<T> operator*(const MatrixDense<T>& other) const {
        if (_n != other._m) throw runtime_error("Matrix dimensions do not allow multiplication");
        MatrixDense<T> result(_m, other._n);
        for (unsigned i = 0; i < _m; ++i) {
            for (unsigned j = 0; j < other._n; ++j) {
                result.at(i, j) = 0;
                for (unsigned k = 0; k < _n; ++k) {
                    result.at(i, j) += at(i, k) * other.at(k, j);
                }
            }
        }
        return result;
    }

    MatrixDense<T> elementWiseMultiply(const MatrixDense<T>& other) const {
        if (_m != other._m || _n != other._n) throw runtime_error("Matrix dimensions must match");
        MatrixDense<T> result(_m, _n);
        for (unsigned i = 0; i < _m; ++i) {
            for (unsigned j = 0; j < _n; ++j) {
                result.at(i, j) = at(i, j) * other.at(i, j);
            }
        }
        return result;
    }

    MatrixDense<T> transpose() const {
        MatrixDense<T> result(_n, _m);
        for (unsigned i = 0; i < _m; ++i) {
            for (unsigned j = 0; j < _n; ++j) {
                result.at(j, i) = at(i, j);
            }
        }
        return result;
    }
};

// Класс диагональной матрицы
template <typename T = double>
class MatrixDiagonal : public Matrix<T> {
private:
    unsigned _n;
    T* data;

public:
    MatrixDiagonal(unsigned n) : _n(n) {
        data = new T[n]();
    }

    ~MatrixDiagonal() {
        delete[] data;
    }

    unsigned getRows() const override {
        return _n;
    }

    unsigned getCols() const override {
        return _n;
    }

    T& at(unsigned row, unsigned col) override {
        if (row >= _n || col >= _n) throw out_of_range("Index out of range");
        // Retourne la diagonale ou 0 si ce n'est pas une diagonale
        return (row == col) ? data[row] : *(new T(0)); // Retour temporaire 0
    }

    const T& at(unsigned row, unsigned col) const override {
        if (row >= _n || col >= _n) throw out_of_range("Index out of range");
        // Retourner 0 si ce n'est pas une diagonale
        static const T zero = T(0);  // Valeur zéro statique
        return (row == col) ? data[row] : zero;
    }

    void print() const override {
        for (unsigned i = 0; i < _n; ++i) {
            for (unsigned j = 0; j < _n; ++j) {
                cout << at(i, j) << " ";
            }
            cout << endl;
        }
    }

    void importFromFile(const string& filename) override {
        ifstream file(filename);
        if (!file) throw runtime_error("Failed to open file");

        string className;
        file >> className;
        if (className != "MatrixDiagonal") throw runtime_error("Incorrect matrix type");

        file >> _n;
        delete[] data;
        data = new T[_n];

        for (unsigned i = 0; i < _n; ++i) {
            file >> data[i];
        }
    }

    void exportToFile(const string& filename) const override {
        ofstream file(filename);
        if (!file) throw runtime_error("Failed to open file");

        file << "MatrixDiagonal\n";
        file << _n << endl;
        for (unsigned i = 0; i < _n; ++i) {
            file << data[i] << " ";
        }
        file << endl;
    }

    MatrixDiagonal<T> operator+(const MatrixDiagonal<T>& other) const {
        if (_n != other._n) throw runtime_error("Matrix dimensions must match");
        MatrixDiagonal<T> result(_n);
        for (unsigned i = 0; i < _n; ++i) {
            result.data[i] = data[i] + other.data[i];
        }
        return result;
    }

    MatrixDiagonal<T> operator-(const MatrixDiagonal<T>& other) const {
        if (_n != other._n) throw runtime_error("Matrix dimensions must match");
        MatrixDiagonal<T> result(_n);
        for (unsigned i = 0; i < _n; ++i) {
            result.data[i] = data[i] - other.data[i];
        }
        return result;
    }

    MatrixDiagonal<T> elementWiseMultiply(const MatrixDiagonal<T>& other) const {
        if (_n != other._n) throw runtime_error("Matrix dimensions must match");
        MatrixDiagonal<T> result(_n);
        for (unsigned i = 0; i < _n; ++i) {
            result.data[i] = data[i] * other.data[i];
        }
        return result;
    }
};

// Основная программа
int main() {
    try {
        // Создание и тестирование матриц
        MatrixDense<double> m1(3, 3);
        m1.at(0, 0) = 1.0;
        m1.at(1, 1) = 2.0;
        m1.at(2, 2) = 3.0;
        cout << "Matrix m1 (Dense):\n";
        m1.print();

        // Создание диагональной матрицы
        MatrixDense<double> m1Transposed = m1.transpose();
        cout << "Transposed Matrix m1:\n";
        m1Transposed.print();

        // Создание диагональной матрицы
        MatrixDiagonal<double> m2(3);
        m2.at(0, 0) = 4.0;
        m2.at(1, 1) = 5.0;
        m2.at(2, 2) = 6.0;
        cout << "Matrix m2 (Diagonal):\n";
        m2.print();

        // Тест
        MatrixDense<double> m3 = m1.elementWiseMultiply(m1);
        cout << "Element-wise multiplication of m1 with itself:\n";
        m3.print();

        MatrixDiagonal<double> m4 = m2 + m2;
        cout << "Sum of m2 with itself (Diagonal):\n";
        m4.print();

        m1.exportToFile("matrix_dense.txt");
        m2.exportToFile("matrix_diagonal.txt");

    }
    catch (const exception& ex) {
        cerr << "Error: " << ex.what() << endl;
    }

    return 0;
}
