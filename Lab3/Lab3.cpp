#include <iostream>
#include <vector>
#include <fstream>
#include <stdexcept>
#include <thread>
#include <mutex>
#include <cmath>
#include <random>
#include <chrono>
#include <functional>
#include <numeric>

std::mutex output_mutex;

// Type de classe template Vector
template <typename T>
class Vector {
private:
    T* data;
    size_t n;
    bool is_initialized;

    void check_initialized() const {
        if (!is_initialized) {
            throw std::runtime_error("Vector is not initialized.");
        }
    }

public:
    // Constructeur
    explicit Vector(size_t size) : n(size), is_initialized(false) {
        if (n > 0) {
            data = new T[n]();
        }
        else {
            throw std::invalid_argument("Size must be greater than 0.");
        }
    }

    // Destructeur
    ~Vector() {
        delete[] data;
    }

    // Initialisation avec une constante
    void initialize_with_constant(T value) {
        for (size_t i = 0; i < n; ++i) {
            data[i] = value;
        }
        is_initialized = true;
    }

    // Initialisation avec des nombres aléatoires
    void initialize_with_random(T min_val, T max_val) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<T> dis(min_val, max_val);

        for (size_t i = 0; i < n; ++i) {
            data[i] = dis(gen);
        }
        is_initialized = true;
    }

    // Exporter vers un fichier
    void export_to_file(const std::string& filename) const {
        check_initialized();
        std::ofstream file(filename);
        if (!file) {
            throw std::runtime_error("Could not open file for writing.");
        }
        for (size_t i = 0; i < n; ++i) {
            file << data[i] << " ";
        }
        file.close();
    }

    // Importer depuis un fichier
    void import_from_file(const std::string& filename) {
        std::ifstream file(filename);
        if (!file) {
            throw std::runtime_error("Could not open file for reading.");
        }
        for (size_t i = 0; i < n && file >> data[i]; ++i) {}
        is_initialized = true;
    }

    // Trouver l'élément minimum (en parallèle)
    std::pair<T, size_t> find_min_parallel() const {
        check_initialized();
        auto start = std::chrono::high_resolution_clock::now();
        T min_val = data[0];
        size_t min_index = 0;

        auto worker = [this, &min_val, &min_index](size_t start_idx, size_t end_idx) {
            for (size_t i = start_idx; i < end_idx; ++i) {
                if (data[i] < min_val) {
                    min_val = data[i];
                    min_index = i;
                }
            }
            };

        size_t num_threads = std::thread::hardware_concurrency();
        size_t chunk_size = n / num_threads;
        std::vector<std::thread> threads;
        for (size_t i = 0; i < num_threads; ++i) {
            size_t start_idx = i * chunk_size;
            size_t end_idx = (i == num_threads - 1) ? n : (i + 1) * chunk_size;
            threads.emplace_back(worker, start_idx, end_idx);
        }
        for (auto& t : threads) {
            t.join();
        }

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        {
            std::lock_guard<std::mutex> lock(output_mutex);
            std::cout << "Parallel Find Min Time: " << duration.count() << " seconds\n";
        }
        return { min_val, min_index };
    }

    // Trouver l'élément minimum (en séquentiel)
    std::pair<T, size_t> find_min_sequential() const {
        check_initialized();
        auto start = std::chrono::high_resolution_clock::now();
        T min_val = data[0];
        size_t min_index = 0;
        for (size_t i = 1; i < n; ++i) {
            if (data[i] < min_val) {
                min_val = data[i];
                min_index = i;
            }
        }
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        {
            std::lock_guard<std::mutex> lock(output_mutex);
            std::cout << "Sequential Find Min Time: " << duration.count() << " seconds\n";
        }
        return { min_val, min_index };
    }

    // Trouver l'élément maximum (en parallèle)
    std::pair<T, size_t> find_max_parallel() const {
        check_initialized();
        auto start = std::chrono::high_resolution_clock::now();
        T max_val = data[0];
        size_t max_index = 0;

        auto worker = [this, &max_val, &max_index](size_t start_idx, size_t end_idx) {
            for (size_t i = start_idx; i < end_idx; ++i) {
                if (data[i] > max_val) {
                    max_val = data[i];
                    max_index = i;
                }
            }
            };

        size_t num_threads = std::thread::hardware_concurrency();
        size_t chunk_size = n / num_threads;
        std::vector<std::thread> threads;
        for (size_t i = 0; i < num_threads; ++i) {
            size_t start_idx = i * chunk_size;
            size_t end_idx = (i == num_threads - 1) ? n : (i + 1) * chunk_size;
            threads.emplace_back(worker, start_idx, end_idx);
        }
        for (auto& t : threads) {
            t.join();
        }

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        {
            std::lock_guard<std::mutex> lock(output_mutex);
            std::cout << "Parallel Find Max Time: " << duration.count() << " seconds\n";
        }
        return { max_val, max_index };
    }

    // Trouver l'élément maximum (en séquentiel)
    std::pair<T, size_t> find_max_sequential() const {
        check_initialized();
        auto start = std::chrono::high_resolution_clock::now();
        T max_val = data[0];
        size_t max_index = 0;
        for (size_t i = 1; i < n; ++i) {
            if (data[i] > max_val) {
                max_val = data[i];
                max_index = i;
            }
        }
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        {
            std::lock_guard<std::mutex> lock(output_mutex);
            std::cout << "Sequential Find Max Time: " << duration.count() << " seconds\n";
        }
        return { max_val, max_index };
    }

    // Calculer la somme des éléments (en parallèle)
    T sum_parallel() const {
        check_initialized();
        auto start = std::chrono::high_resolution_clock::now();
        T result = 0;

        auto worker = [this, &result](size_t start_idx, size_t end_idx) {
            T partial_sum = 0;
            for (size_t i = start_idx; i < end_idx; ++i) {
                partial_sum += data[i];
            }
            result += partial_sum;
            };

        size_t num_threads = std::thread::hardware_concurrency();
        size_t chunk_size = n / num_threads;
        std::vector<std::thread> threads;
        for (size_t i = 0; i < num_threads; ++i) {
            size_t start_idx = i * chunk_size;
            size_t end_idx = (i == num_threads - 1) ? n : (i + 1) * chunk_size;
            threads.emplace_back(worker, start_idx, end_idx);
        }
        for (auto& t : threads) {
            t.join();
        }

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        {
            std::lock_guard<std::mutex> lock(output_mutex);
            std::cout << "Parallel Sum Time: " << duration.count() << " seconds\n";
        }
        return result;
    }

    // Calculer la somme des éléments (en séquentiel)
    T sum_sequential() const {
        check_initialized();
        auto start = std::chrono::high_resolution_clock::now();
        T result = 0;
        for (size_t i = 0; i < n; ++i) {
            result += data[i];
        }
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        {
            std::lock_guard<std::mutex> lock(output_mutex);
            std::cout << "Sequential Sum Time: " << duration.count() << " seconds\n";
        }
        return result;
    }

    // Calculer la norme Euclidienne (en parallèle)
    T euclidean_norm_parallel() const {
        check_initialized();
        auto start = std::chrono::high_resolution_clock::now();
        T result = 0;

        auto worker = [this, &result](size_t start_idx, size_t end_idx) {
            T partial_sum = 0;
            for (size_t i = start_idx; i < end_idx; ++i) {
                partial_sum += data[i] * data[i];
            }
            result += partial_sum;
            };

        size_t num_threads = std::thread::hardware_concurrency();
        size_t chunk_size = n / num_threads;
        std::vector<std::thread> threads;
        for (size_t i = 0; i < num_threads; ++i) {
            size_t start_idx = i * chunk_size;
            size_t end_idx = (i == num_threads - 1) ? n : (i + 1) * chunk_size;
            threads.emplace_back(worker, start_idx, end_idx);
        }
        for (auto& t : threads) {
            t.join();
        }

        result = std::sqrt(result);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        {
            std::lock_guard<std::mutex> lock(output_mutex);
            std::cout << "Parallel Euclidean Norm Time: " << duration.count() << " seconds\n";
        }
        return result;
    }

    // Calculer la norme Euclidienne (en séquentiel)
    T euclidean_norm_sequential() const {
        check_initialized();
        auto start = std::chrono::high_resolution_clock::now();
        T result = 0;
        for (size_t i = 0; i < n; ++i) {
            result += data[i] * data[i];
        }
        result = std::sqrt(result);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        {
            std::lock_guard<std::mutex> lock(output_mutex);
            std::cout << "Sequential Euclidean Norm Time: " << duration.count() << " seconds\n";
        }
        return result;
    }

    // Calculer la norme Manhattan (en parallèle)
    T manhattan_norm_parallel() const {
        check_initialized();
        auto start = std::chrono::high_resolution_clock::now();
        T result = 0;

        auto worker = [this, &result](size_t start_idx, size_t end_idx) {
            T partial_sum = 0;
            for (size_t i = start_idx; i < end_idx; ++i) {
                partial_sum += std::abs(data[i]);
            }
            result += partial_sum;
            };

        size_t num_threads = std::thread::hardware_concurrency();
        size_t chunk_size = n / num_threads;
        std::vector<std::thread> threads;
        for (size_t i = 0; i < num_threads; ++i) {
            size_t start_idx = i * chunk_size;
            size_t end_idx = (i == num_threads - 1) ? n : (i + 1) * chunk_size;
            threads.emplace_back(worker, start_idx, end_idx);
        }
        for (auto& t : threads) {
            t.join();
        }

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        {
            std::lock_guard<std::mutex> lock(output_mutex);
            std::cout << "Parallel Manhattan Norm Time: " << duration.count() << " seconds\n";
        }
        return result;
    }

    // Calculer la norme Manhattan (en séquentiel)
    T manhattan_norm_sequential() const {
        check_initialized();
        auto start = std::chrono::high_resolution_clock::now();
        T result = 0;
        for (size_t i = 0; i < n; ++i) {
            result += std::abs(data[i]);
        }
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        {
            std::lock_guard<std::mutex> lock(output_mutex);
            std::cout << "Sequential Manhattan Norm Time: " << duration.count() << " seconds\n";
        }
        return result;
    }
};

int main() {
    try {
        Vector<double> v(1000000);
        v.initialize_with_random(0.0, 10.0);

        std::ofstream output_file("vector_performance_results.txt");
        if (!output_file) {
            throw std::runtime_error("Could not open file for writing.");
        }

        auto min_parallel = v.find_min_parallel();
        output_file << "Parallel Min Time: " << min_parallel.first << " seconds\n";

        auto min_sequential = v.find_min_sequential();
        output_file << "Sequential Min Time: " << min_sequential.first << " seconds\n";

        auto max_parallel = v.find_max_parallel();
        output_file << "Parallel Max Time: " << max_parallel.first << " seconds\n";

        auto max_sequential = v.find_max_sequential();
        output_file << "Sequential Max Time: " << max_sequential.first << " seconds\n";

        double sum_parallel = v.sum_parallel();
        output_file << "Parallel Sum Time: " << sum_parallel << " seconds\n";

        double sum_sequential = v.sum_sequential();
        output_file << "Sequential Sum Time: " << sum_sequential << " seconds\n";

        double norm_parallel = v.euclidean_norm_parallel();
        output_file << "Parallel Euclidean Norm Time: " << norm_parallel << " seconds\n";

        double norm_sequential = v.euclidean_norm_sequential();
        output_file << "Sequential Euclidean Norm Time: " << norm_sequential << " seconds\n";

        double manhattan_parallel = v.manhattan_norm_parallel();
        output_file << "Parallel Manhattan Norm Time: " << manhattan_parallel << " seconds\n";

        double manhattan_sequential = v.manhattan_norm_sequential();
        output_file << "Sequential Manhattan Norm Time: " << manhattan_sequential << " seconds\n";

        output_file.close();
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }
}
