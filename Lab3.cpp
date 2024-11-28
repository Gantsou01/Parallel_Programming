#include <iostream>
#include <vector>
#include <stdexcept>
#include <random>
#include <chrono>
#include <thread>
#include <numeric>
#include <cmath>
#include <algorithm>

template <typename T>
class Vector {
private:
    T* data;
    size_t n;
    bool is_initialized;

public:
    // Constructeur avec la taille
    Vector(size_t size) : n(size), is_initialized(false) {
        data = new T[n];
        std::cout << "Vector of size " << n << " created.\n";  // Debug
    }

    // Destructeur
    ~Vector() {
        if (data) {
            delete[] data;
            std::cout << "Memory released.\n";  // Debug
        }
    }

    // Initialisation avec une constante
    void initialize_with_constant(T value) {
        if (is_initialized) {
            throw std::runtime_error("Vector already initialized");
        }
        for (size_t i = 0; i < n; ++i) {
            data[i] = value;
        }
        is_initialized = true;
        std::cout << "Vector initialized with constant " << value << ".\n";  // Debug
    }

    // Initialisation avec des nombres aléatoires dans un intervalle
    void initialize_with_random(T lower, T upper) {
        if (is_initialized) {
            throw std::runtime_error("Vector already initialized");
        }
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<T> dis(lower, upper);

        for (size_t i = 0; i < n; ++i) {
            data[i] = dis(gen);
        }
        is_initialized = true;
        std::cout << "Vector initialized with random values between " << lower << " and " << upper << ".\n";  // Debug
    }

    // Vérifier si le vecteur est initialisé
    void check_initialized() const {
        if (!is_initialized) {
            throw std::runtime_error("Vector is not initialized");
        }
    }

    // Trouver l'élément minimal et son indice
    std::pair<T, size_t> find_min() {
        check_initialized();
        T min_value = data[0];
        size_t min_index = 0;
        for (size_t i = 1; i < n; ++i) {
            if (data[i] < min_value) {
                min_value = data[i];
                min_index = i;
            }
        }
        return { min_value, min_index };
    }

    // Trouver l'élément maximal et son indice
    std::pair<T, size_t> find_max() {
        check_initialized();
        T max_value = data[0];
        size_t max_index = 0;
        for (size_t i = 1; i < n; ++i) {
            if (data[i] > max_value) {
                max_value = data[i];
                max_index = i;
            }
        }
        return { max_value, max_index };
    }

    // Calculer la moyenne
    T mean() {
        check_initialized();
        T sum = std::accumulate(data, data + n, T(0));
        return sum / static_cast<T>(n);
    }

    // Somme des éléments
    T sum() {
        check_initialized();
        return std::accumulate(data, data + n, T(0));
    }

    // Norme Euclidienne (L2)
    T euclidean_norm() {
        check_initialized();
        T sum = std::accumulate(data, data + n, T(0), [](T acc, T val) { return acc + val * val; });
        return std::sqrt(sum);
    }

    // Norme de Manhattan (L1)
    T manhattan_norm() {
        check_initialized();
        T sum = std::accumulate(data, data + n, T(0), [](T acc, T val) { return acc + std::abs(val); });
        return sum;
    }

    // Produit scalaire avec un autre vecteur
    T dot_product(const Vector<T>& other) {
        check_initialized();
        other.check_initialized();
        if (n != other.n) {
            throw std::invalid_argument("Vectors must have the same size");
        }
        T product = T(0);
        for (size_t i = 0; i < n; ++i) {
            product += data[i] * other.data[i];
        }
        return product;
    }

    // Exécution parallèle
    template <typename Func>
    double parallel_execution(Func func, size_t num_threads) {
        check_initialized();  // Vérifie si le vecteur est initialisé
        auto start = std::chrono::high_resolution_clock::now();  // Mesure du temps de début

        // Diviser le travail entre les threads
        std::vector<std::thread> threads;
        size_t block_size = n / num_threads;

        // Tableau pour stocker les résultats partiels
        std::vector<T> results(num_threads, T(0));

        // Création des threads
        for (size_t t = 0; t < num_threads; ++t) {
            threads.push_back(std::thread([&, t] {
                size_t start_idx = t * block_size;
                size_t end_idx = (t == num_threads - 1) ? n : (t + 1) * block_size;

                // Effectuer la fonction pour chaque segment
                results[t] = func(start_idx, end_idx);
                }));
        }

        // Attendre la fin des threads
        for (auto& t : threads) {
            t.join();
        }

        // Combiner les résultats
        T final_result = std::accumulate(results.begin(), results.end(), T(0));

        auto end = std::chrono::high_resolution_clock::now();  // Mesure du temps de fin
        std::chrono::duration<double> duration = end - start;  // Calcul de la durée

        return duration.count();  // Retourne le temps en secondes
    }
};

// Fonction main en dehors de la classe
int main() {
    size_t n = 1000000;  // Taille du vecteur
    Vector<int> vec(n);

    try {
        // Initialisation avec des valeurs aléatoires
        vec.initialize_with_random(-1000, 1000);

        // Trouver l'élément minimal
        auto min_pair = vec.find_min();
        std::cout << "Min: " << min_pair.first << " at index " << min_pair.second << std::endl;

        // Exécution parallèle pour la somme des éléments
        double time = vec.parallel_execution([&](size_t start, size_t end) -> int {
            int partial_sum = 0;
            for (size_t i = start; i < end; ++i) {
                partial_sum += vec.data[i];
            }
            return partial_sum;
            }, 4);  // Utilisation de 4 threads

        std::cout << "Parallel execution time: " << time << " seconds" << std::endl;

    }
    catch (const std::exception& e) {
        std::cout << "Exception: " << e.what() << std::endl;
    }

    return 0;
}
