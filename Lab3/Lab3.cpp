#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <execution> // Pour la parallélisation en C++17

template <typename T>
class Vector {
private:
    std::vector<T> data;  // Contient les éléments du vecteur

public:
    // Constructeur
    Vector(std::initializer_list<T> init) : data(init) {}

    // Accesseur pour récupérer un élément
    T operator[](std::size_t index) const {
        return data[index];
    }

    // Retourne la taille du vecteur
    std::size_t size() const {
        return data.size();
    }

    // Surcharge de l'opérateur << pour afficher le vecteur
    friend std::ostream& operator<<(std::ostream& os, const Vector<T>& v) {
        for (const auto& elem : v.data) {
            os << elem << " ";
        }
        return os;
    }

    // Calcul de la somme des éléments (séquentielle)
    T sum_sequential() const {
        return std::accumulate(data.begin(), data.end(), T(0));
    }

    // Calcul de la norme Euclidienne (séquentielle)
    T euclidean_norm_sequential() const {
        T sum_of_squares = std::accumulate(data.begin(), data.end(), T(0),
            [](T sum, T val) { return sum + val * val; });
        return std::sqrt(sum_of_squares);
    }

    // Calcul de la norme Manhattan (séquentielle)
    T manhattan_norm_sequential() const {
        return std::accumulate(data.begin(), data.end(), T(0),
            [](T sum, T val) { return sum + std::abs(val); });
    }

    // Calcul de la somme des éléments (parallèle)
    T sum_parallel() const {
        return std::reduce(std::execution::par, data.begin(), data.end(), T(0));
    }

    // Calcul de la norme Euclidienne (parallèle)
    T euclidean_norm_parallel() const {
        T sum_of_squares = std::reduce(std::execution::par, data.begin(), data.end(), T(0),
            [](T sum, T val) { return sum + val * val; });
        return std::sqrt(sum_of_squares);
    }

    // Calcul de la norme Manhattan (parallèle)
    T manhattan_norm_parallel() const {
        return std::reduce(std::execution::par, data.begin(), data.end(), T(0),
            [](T sum, T val) { return sum + std::abs(val); });
    }

    // Trouver la valeur maximale (séquentielle)
    T find_max_sequential() const {
        return *std::max_element(data.begin(), data.end());
    }
};

int main() {
    // Création d'un vecteur d'exemple
    Vector<float> v = { 1.0, 2.0, 3.0, 4.0, 5.0 };

    // Affichage du vecteur
    std::cout << "Vecteur: " << v << std::endl;

    // Appels aux différentes méthodes
    std::cout << "Somme séquentielle: " << v.sum_sequential() << std::endl;
    std::cout << "Norme Euclidienne séquentielle: " << v.euclidean_norm_sequential() << std::endl;
    std::cout << "Norme Manhattan séquentielle: " << v.manhattan_norm_sequential() << std::endl;

    std::cout << "Somme parallèle: " << v.sum_parallel() << std::endl;
    std::cout << "Norme Euclidienne parallèle: " << v.euclidean_norm_parallel() << std::endl;
    std::cout << "Norme Manhattan parallèle: " << v.manhattan_norm_parallel() << std::endl;

    std::cout << "Valeur maximale séquentielle: " << v.find_max_sequential() << std::endl;

    return 0;
}
