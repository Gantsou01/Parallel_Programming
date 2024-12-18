#include <iostream>
#include <vector>
#include <stdexcept>
#include <random>

// === Classe de base ===
template <typename T = double>
class MatrixBase {
public:
    virtual ~MatrixBase() = default;
    virtual T& operator()(unsigned i, unsigned j) = 0;
    virtual T operator()(unsigned i, unsigned j) const = 0;
    virtual void display() const = 0;
};

// === Matrice dense ===
template <typename T = double>
class MatrixDense : public MatrixBase<T> {
private:
    unsigned _m, _n;
    T* _data;

public:
    MatrixDense(unsigned m, unsigned n) : _m(m), _n(n), _data(new T[m * n]()) {}
    ~MatrixDense() { delete[] _data; }

    T& operator()(unsigned i, unsigned j) override { return _data[i * _n + j]; }
    T operator()(unsigned i, unsigned j) const override { return _data[i * _n + j]; }

    MatrixDense operator+(const MatrixDense& other) const {
        if (_m != other._m || _n != other._n) throw std::invalid_argument("Taille incompatible pour addition.");
        MatrixDense result(_m, _n);
        for (unsigned i = 0; i < _m; ++i)
            for (unsigned j = 0; j < _n; ++j)
                result(i, j) = (*this)(i, j) + other(i, j);
        return result;
    }

    MatrixDense operator*(const MatrixDense& other) const {
        if (_n != other._m) throw std::invalid_argument("Taille incompatible pour multiplication.");
        MatrixDense result(_m, other._n);
        for (unsigned i = 0; i < _m; ++i)
            for (unsigned j = 0; j < other._n; ++j)
                for (unsigned k = 0; k < _n; ++k)
                    result(i, j) += (*this)(i, k) * other(k, j);
        return result;
    }

    MatrixDense operator*(T scalar) const {
        MatrixDense result(_m, _n);
        for (unsigned i = 0; i < _m; ++i)
            for (unsigned j = 0; j < _n; ++j)
                result(i, j) = (*this)(i, j) * scalar;
        return result;
    }

    MatrixDense transpose() const {
        MatrixDense result(_n, _m);
        for (unsigned i = 0; i < _m; ++i)
            for (unsigned j = 0; j < _n; ++j)
                result(j, i) = (*this)(i, j);
        return result;
    }

    void display() const override {
        for (unsigned i = 0; i < _m; ++i) {
            for (unsigned j = 0; j < _n; ++j) {
                std::cout << (*this)(i, j) << " ";
            }
            std::cout << "\n";
        }
    }
};


// === Matrice diagonale ===
template <typename T = double>
class MatrixDiagonal : public MatrixBase<T> {
private:
    unsigned _m, _n;
    std::vector<int> _diag_indices;
    T* _data;

public:
    MatrixDiagonal(unsigned m, unsigned n, const std::vector<int>& diag_indices)
        : _m(m), _n(n), _diag_indices(diag_indices), _data(new T[diag_indices.size() * m]()) {}
    ~MatrixDiagonal() { delete[] _data; }

    T& operator()(unsigned i, unsigned j) override {
        for (size_t k = 0; k < _diag_indices.size(); ++k) {
            if (j == i + _diag_indices[k]) return _data[k * _m + i];
        }
        static T zero = T();
        return zero;
    }

    T operator()(unsigned i, unsigned j) const override {
        for (size_t k = 0; k < _diag_indices.size(); ++k) {
            if (j == i + _diag_indices[k]) return _data[k * _m + i];
        }
        return T();
    }

    MatrixDiagonal operator+(const MatrixDiagonal& other) const {
        if (_m != other._m || _n != other._n) throw std::invalid_argument("Taille incompatible pour addition.");
        if (_diag_indices != other._diag_indices) throw std::invalid_argument("Indices diagonaux incompatibles.");
        MatrixDiagonal result(_m, _n, _diag_indices);
        for (size_t k = 0; k < _diag_indices.size(); ++k)
            for (unsigned i = 0; i < _m; ++i)
                result._data[k * _m + i] = _data[k * _m + i] + other._data[k * _m + i];
        return result;
    }

    MatrixDiagonal operator*(const MatrixDiagonal& other) const {
        if (_m != other._m || _n != other._n) throw std::invalid_argument("Taille incompatible pour multiplication.");
        if (_diag_indices != other._diag_indices) throw std::invalid_argument("Indices diagonaux incompatibles.");
        MatrixDiagonal result(_m, _n, _diag_indices);
        for (size_t k = 0; k < _diag_indices.size(); ++k)
            for (unsigned i = 0; i < _m; ++i)
                result._data[k * _m + i] = _data[k * _m + i] * other._data[k * _m + i];
        return result;
    }

    MatrixDiagonal transpose() const {
        std::vector<int> neg_diag_indices(_diag_indices.size());
        for (size_t i = 0; i < _diag_indices.size(); ++i) neg_diag_indices[i] = -_diag_indices[i];
        MatrixDiagonal result(_m, _n, neg_diag_indices);
        for (size_t k = 0; k < _diag_indices.size(); ++k)
            for (unsigned i = 0; i < _m; ++i)
                result._data[k * _m + i] = _data[k * _m + i];
        return result;
    }

    MatrixDiagonal operator*(T scalar) const {
        MatrixDiagonal result(_m, _n, _diag_indices);
        for (size_t k = 0; k < _diag_indices.size(); ++k)
            for (unsigned i = 0; i < _m; ++i)
                result._data[k * _m + i] = _data[k * _m + i] * scalar;
        return result;
    }

    void display() const override {
        for (unsigned i = 0; i < _m; ++i) {
            for (unsigned j = 0; j < _n; ++j) {
                std::cout << (*this)(i, j) << " ";
            }
            std::cout << "\n";
        }
    }
};

// === Matrice bloc ===
template <typename T = double>
class MatrixBlock : public MatrixBase<T> {
private:
    unsigned _rows, _cols;                    // Nombre de blocs
    unsigned _block_m, _block_n;              // Taille de chaque bloc
    std::vector<MatrixDense<T>*> _blocks;     // Pointeurs vers des blocs denses

public:
    MatrixBlock(unsigned rows, unsigned cols, unsigned block_m, unsigned block_n)
        : _rows(rows), _cols(cols), _block_m(block_m), _block_n(block_n),
        _blocks(rows* cols, nullptr) {}

    ~MatrixBlock() {
        for (auto block : _blocks) delete block;
    }

    MatrixDense<T>& block(unsigned i, unsigned j) {
        if (!_blocks[i * _cols + j]) throw std::invalid_argument("Bloc inexistant.");
        return *_blocks[i * _cols + j];
    }

    void set_block(unsigned i, unsigned j, MatrixDense<T>* block) {
        _blocks[i * _cols + j] = block;
    }

    T& operator()(unsigned i, unsigned j) override {
        unsigned block_i = i / _block_m;
        unsigned block_j = j / _block_n;
        unsigned local_i = i % _block_m;
        unsigned local_j = j % _block_n;
        return block(block_i, block_j)(local_i, local_j);
    }

    T operator()(unsigned i, unsigned j) const override {
        unsigned block_i = i / _block_m;
        unsigned block_j = j / _block_n;
        unsigned local_i = i % _block_m;
        unsigned local_j = j % _block_n;
        return block(block_i, block_j)(local_i, local_j);
    }

    void display() const override {
        for (unsigned i = 0; i < _rows * _block_m; ++i) {
            for (unsigned j = 0; j < _cols * _block_n; ++j) {
                std::cout << (*this)(i, j) << " ";
            }
            std::cout << "\n";
        }
    }
};

// === Fonction principale ===
int main() {
    // === Matrices Denses ===
    std::cout << "=== Matrices Denses ===\n";
    MatrixDense<double> A(3, 3), B(3, 3);
    for (unsigned i = 0; i < 3; ++i)
        for (unsigned j = 0; j < 3; ++j) {
            A(i, j) = i + j;
            B(i, j) = i * j;
        }

    std::cout << "Matrix A: \n";
    A.display();
    std::cout << "Matrix B: \n";
    B.display();

    // Addition
    auto result_add = A + B;
    std::cout << "A + B: \n";
    result_add.display();

    // Multiplication par éléments
    auto result_elem_mult = A * B;
    std::cout << "A * B (élément par élément): \n";
    result_elem_mult.display();

    // Multiplication matricielle
    auto result_mat_mult = A * B.transpose();
    std::cout << "A * B (matr. complète): \n";
    result_mat_mult.display();

    // Transposition
    auto result_transpose = A.transpose();
    std::cout << "Transpose de A: \n";
    result_transpose.display();

    // Multiplication par un scalaire
    auto result_scalar_mult = A * 2.0;
    std::cout << "A * 2.0: \n";
    result_scalar_mult.display();

    // === Matrices Diagonales ===
    std::cout << "=== Matrices Diagonales ===\n";
    MatrixDiagonal<double> C(3, 3, { 0 });
    MatrixDiagonal<double> D(3, 3, { 0 });
    C(0, 0) = 5; C(1, 1) = 10; C(2, 2) = 15;
    D(0, 0) = 2; D(1, 1) = 4; D(2, 2) = 8;

    std::cout << "Matrix C: \n";
    C.display();
    std::cout << "Matrix D: \n";
    D.display();

    // Addition
    auto result_diag_add = C + D;
    std::cout << "C + D: \n";
    result_diag_add.display();

    // Multiplication par éléments
    auto result_diag_elem_mult = C * D;
    std::cout << "C * D (élément par élément): \n";
    result_diag_elem_mult.display();

    // Transposition
    auto result_diag_transpose = C.transpose();
    std::cout << "Transpose de C: \n";
    result_diag_transpose.display();

    // Multiplication par un scalaire
    auto result_diag_scalar_mult = C * 3.0;
    std::cout << "C * 3.0: \n";
    result_diag_scalar_mult.display();

    return 0;
}
