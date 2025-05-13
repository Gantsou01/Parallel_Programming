#include <mpi.h>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdio>
#include <omp.h>

#define _CRT_SECURE_NO_WARNINGS

// Задание 15: Программа «I am!»
void task15() {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::cout << "I am " << rank << " process from " << size << " processes!" << std::endl;
}

// Задание 16: Программа «На первый-второй рассчитайся!»
void task16() {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        std::cout << size << " processes." << std::endl;
    }

    if (rank % 2 == 0) {
        std::cout << "I am " << rank << " process: FIRST!" << std::endl;
    }
    else {
        std::cout << "I am " << rank << " process: SECOND!" << std::endl;
    }
}

// Задание 17: Простые блокирующие обмены
void task17() {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        int message = 45;
        MPI_Send(&message, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
    }
    else if (rank == 1) {
        int message;
        MPI_Recv(&message, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        std::cout << "receive message '" << message << "'" << std::endl;
    }
}

// Задание 18: Схема «эстафетная палочка»
void task18() {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int message;
    if (rank == 0) {
        message = 0;
        MPI_Send(&message, 1, MPI_INT, (rank + 1) % size, 0, MPI_COMM_WORLD);
        MPI_Recv(&message, 1, MPI_INT, size - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        std::cout << "[" << rank << "] : receive message '" << message << "'" << std::endl;
    }
    else {
        MPI_Recv(&message, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        message += 1;
        MPI_Send(&message, 1, MPI_INT, (rank + 1) % size, 0, MPI_COMM_WORLD);
        std::cout << "[" << rank << "] : receive message '" << message - 1 << "'" << std::endl;
    }
}

// Задание 19: Схема «мастер-рабочие»
void task19() {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        for (int i = 1; i < size; i++) {
            int message;
            MPI_Recv(&message, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            std::cout << "receive message '" << message << "'" << std::endl;
        }
    }
    else {
        MPI_Send(&rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
}

// Задание 20: Простые неблокирующие обмены
void task20() {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        int message = 45;
        MPI_Request request;
        MPI_Isend(&message, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, &request);
        MPI_Wait(&request, MPI_STATUS_IGNORE);
    }
    else if (rank == 1) {
        int message;
        MPI_Request request;
        MPI_Irecv(&message, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &request);
        MPI_Wait(&request, MPI_STATUS_IGNORE);
        std::cout << "receive message '" << message << "'" << std::endl;
    }
}

// Задание 21: Схема «сдвиг по кольцу»
void task21() {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int send_msg = rank;
    int recv_msg;
    MPI_Status status;

    MPI_Sendrecv(&send_msg, 1, MPI_INT, (rank + 1) % size, 0,
        &recv_msg, 1, MPI_INT, (rank - 1 + size) % size, 0,
        MPI_COMM_WORLD, &status);

    std::cout << "[" << rank << "] : receive message '" << recv_msg << "'" << std::endl;
}

// Задание 22: Схема «каждый каждому»
void task22() {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::vector<int> messages(size, -1);
    std::vector<MPI_Request> requests(size * 2);

    // Отправка сообщений всем процессам
    for (int i = 0, idx = 0; i < size; i++) {
        if (i != rank) {
            MPI_Isend(&rank, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &requests[idx++]);
        }
    }

    // Прием сообщений от всех процессов
    for (int i = 0, idx = size - 1; i < size; i++) {
        if (i != rank) {
            MPI_Irecv(&messages[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &requests[idx++]);
        }
    }

    MPI_Waitall(2 * (size - 1), requests.data(), MPI_STATUSES_IGNORE);

    // Вывод полученных сообщений
    for (int i = 0; i < size; i++) {
        if (i != rank) {
            std::cout << "[" << rank << "] : receive message '" << messages[i] << "' from " << i << std::endl;
        }
    }
}

// Задание 23: Широковещательная рассылка данных
void task23() {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int n = 0;
    std::string str;

    if (rank == 0) {
        std::cout << "Enter string length and the string: ";
        std::cin >> n;
        std::cin >> str;
        str.resize(n);
    }

    // Рассылка длины строки
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank != 0) {
        str.resize(n);
    }

    // Рассылка самой строки
    MPI_Bcast(&str[0], n, MPI_CHAR, 0, MPI_COMM_WORLD);

    // Подсчет символов
    int local_counts[26] = { 0 };
    const int chunk_size = n / size;
    const int start = rank * chunk_size;
    const int end = (rank == size - 1) ? n : start + chunk_size;

    for (int i = start; i < end; i++) {
        if (str[i] >= 'a' && str[i] <= 'z') {
            local_counts[str[i] - 'a']++;
        }
    }

    // Сбор результатов на процессе 0
    int global_counts[26] = { 0 };
    MPI_Reduce(local_counts, global_counts, 26, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    // Вывод результатов
    if (rank == 0) {
        for (int i = 0; i < 26; i++) {
            if (global_counts[i] > 0) {
                std::cout << static_cast<char>('a' + i) << " = " << global_counts[i] << std::endl;
            }
        }
    }
}

// Задание 24: Операции редукции (вычисление π)
void task24(int N = 1000000000) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (N <= 0) {
        if (rank == 0) {
            std::cerr << "Error: N must be positive" << std::endl;
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Рассылка N всем процессам
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Распределение работы между процессами
    const int chunk_size = N / size;
    const int start = rank * chunk_size;
    const int end = (rank == size - 1) ? N : start + chunk_size;

    double local_sum = 0.0;
    for (int i = start; i < end; i++) {
        const double x = (i + 0.5) / N;
        local_sum += 4.0 / (1.0 + x * x);
    }

    double global_sum;
    MPI_Reduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    // Вывод результата
    if (rank == 0) {
        const double pi = global_sum / N;
        std::cout << "Calculated pi: " << pi << std::endl;
    }
}

// Задание 25: Умножение матриц
void task25() {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int n = 0;
    std::vector<double> A, B, C;

    if (rank == 0) {
        std::cout << "Enter matrix size n: ";
        std::cin >> n;
        A.resize(n * n);
        B.resize(n * n);
        C.resize(n * n);

        std::cout << "Enter matrix A (row-wise): ";
        for (int i = 0; i < n * n; i++) std::cin >> A[i];

        std::cout << "Enter matrix B (row-wise): ";
        for (int i = 0; i < n * n; i++) std::cin >> B[i];
    }

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Распределение данных
    const int elements_per_proc = n * n / size;
    std::vector<double> local_A(elements_per_proc), local_B(elements_per_proc), local_C(elements_per_proc, 0);
    MPI_Scatter(A.data(), elements_per_proc, MPI_DOUBLE, local_A.data(), elements_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(B.data(), elements_per_proc, MPI_DOUBLE, local_B.data(), elements_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Локальное умножение
    for (int i = 0; i < n / size; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                local_C[i * n + j] += local_A[i * n + k] * local_B[k * n + j];
            }
        }
    }

    // Сбор результатов
    MPI_Gather(local_C.data(), elements_per_proc, MPI_DOUBLE, C.data(), elements_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        std::cout << "Result matrix C:\n";
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                std::cout << C[i * n + j] << " ";
            }
            std::cout << "\n";
        }
    }
}

// Задание 26: Группы и коммуникаторы
void task26(const std::string& message = "A") {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::string msg;
    if (rank == 0) {
        msg = message;
    }

    // Создание группы для процессов с четными номерами
    MPI_Group world_group, even_group;
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);

    std::vector<int> even_ranks;
    for (int i = 0; i < size; i += 2) {
        even_ranks.push_back(i);
    }

    MPI_Group_incl(world_group, static_cast<int>(even_ranks.size()), even_ranks.data(), &even_group);

    // Создание нового коммуникатора
    MPI_Comm even_comm;
    MPI_Comm_create(MPI_COMM_WORLD, even_group, &even_comm);

    // Рассылка сообщения в новой группе
    int new_rank = MPI_UNDEFINED, new_size = 0;
    std::string received_msg;
    if (even_comm != MPI_COMM_NULL) {
        MPI_Comm_rank(even_comm, &new_rank);
        MPI_Comm_size(even_comm, &new_size);

        // Рассылка длины сообщения
        int msg_len = static_cast<int>(msg.size());
        MPI_Bcast(&msg_len, 1, MPI_INT, 0, even_comm);

        // Рассылка самого сообщения
        received_msg.resize(msg_len);
        if (rank == 0) received_msg = msg;
        MPI_Bcast(&received_msg[0], msg_len, MPI_CHAR, 0, even_comm);
    }

    // Вывод информации
    if (even_comm != MPI_COMM_NULL) {
        std::cout << "MPI_COMM_WORLD: " << rank << " from " << size
            << ". New comm: " << new_rank << " from " << new_size
            << ". Message = " << received_msg << std::endl;
    }
    else {
        std::cout << "MPI_COMM_WORLD: " << rank << " from " << size
            << ". New comm: no from no. Message = no" << std::endl;
    }

    // Освобождение ресурсов
    if (even_comm != MPI_COMM_NULL) {
        MPI_Comm_free(&even_comm);
    }
    MPI_Group_free(&even_group);
    MPI_Group_free(&world_group);
}

// Задание 27: Динамическое создание процессов
void task27(int num_processes = 3) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        if (num_processes <= 0) {
            std::cerr << "Error: Number of processes must be positive" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
            return;
        }

        const char* command = "mpiexec";
        char args[256];
        snprintf(args, sizeof(args), "-n %d ConsoleApplication2MPI.exe 99", num_processes);
        char* argv[] = { const_cast<char*>(command), args, nullptr };

        MPI_Comm intercomm;
        int error_codes[1] = { 0 };

        int result = MPI_Comm_spawn(
            command,
            argv,
            num_processes,
            MPI_INFO_NULL,
            0,
            MPI_COMM_SELF,
            &intercomm,
            error_codes
        );

        if (result != MPI_SUCCESS) {
            std::cerr << "MPI_Comm_spawn failed with error code: " << result << std::endl;
            MPI_Abort(MPI_COMM_WORLD, result);
        }

        std::cout << "I am 0 process from " << size << " processes! My parent is none.\n";
    }
    else {
        std::cout << "I am " << rank << " process from " << size
            << " processes! My parent is none.\n";
    }
}

// Задание 28: Односторонние коммуникации (вычисление π)
void task28(int N = 1000000000) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (N <= 0) {
        if (rank == 0) {
            std::cerr << "Error: N must be positive" << std::endl;
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    MPI_Win win;
    double pi_part = 0.0, pi_total = 0.0;

    MPI_Win_create(&pi_part, sizeof(double), sizeof(double), MPI_INFO_NULL,
        MPI_COMM_WORLD, &win);

    const int chunk_size = N / size;
    const int start = rank * chunk_size;
    const int end = (rank == size - 1) ? N : start + chunk_size;

    for (int i = start; i < end; i++) {
        const double x = (i + 0.5) / N;
        pi_part += 4.0 / (1.0 + x * x);
    }
    pi_part /= N;

    MPI_Win_fence(0, win);
    if (rank != 0) {
        MPI_Accumulate(&pi_part, 1, MPI_DOUBLE, 0, 0, 1, MPI_DOUBLE, MPI_SUM, win);
    }
    MPI_Win_fence(0, win);

    if (rank == 0) {
        pi_total = pi_part;
        std::cout << "Calculated pi: " << pi_total << std::endl;
    }

    MPI_Win_free(&win);
}

// Задание 29: Исследование масштабируемости
void task29() {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        std::cout << "Task 29: Scalability study template\n";
        std::cout << "| № | Program      | Parameter N | Threads | Time (sec) |\n";
        std::cout << "|---|--------------|-------------|---------|------------|\n";
        // Здесь можно добавить реальные измерения времени
    }
}

// Задание 30: Проект с поддержкой MPI+OpenMP
void task30() {
#ifdef _OPENMP
#pragma omp parallel
    {
        const int thread_num = omp_get_thread_num();
        const int threads_count = omp_get_num_threads();
#pragma omp critical
        {
            std::cout << "OpenMP thread " << thread_num << " of " << threads_count << "\n";
        }
    }
#endif

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::cout << "MPI process " << rank << " of " << size << " (with OpenMP support)\n";
}

// Задание 31: Программа «I am» (MPI+OpenMP)
void task31(int num_threads = 3) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (num_threads <= 0) {
        if (rank == 0) {
            std::cerr << "Error: number of threads must be positive" << std::endl;
        }
        return;
    }

#ifdef _OPENMP
#pragma omp parallel num_threads(num_threads)
    {
        const int thread_num = omp_get_thread_num();
        const int total_hybrid_threads = num_threads * size;
#pragma omp critical
        {
            std::cout << "I am " << thread_num << " thread from " << rank
                << " process. Number of hybrid threads = " << total_hybrid_threads << std::endl;
        }
    }
#else
    std::cout << "I am " << rank << " process (OpenMP not enabled)" << std::endl;
#endif
}

// Задание 32: Программа «Число π» (MPI+OpenMP)
void task32(int N = 1000000000) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (N <= 0) {
        if (rank == 0) {
            std::cerr << "Error: N must be positive" << std::endl;
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    const int chunk_size = N / size;
    const int start = rank * chunk_size;
    const int end = (rank == size - 1) ? N : start + chunk_size;

    double local_sum = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:local_sum)
#endif
    for (int i = start; i < end; i++) {
        const double x = (i + 0.5) / N;
        local_sum += 4.0 / (1.0 + x * x);
    }

    double global_sum;
    MPI_Reduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        const double pi = global_sum / N;
        std::cout << "Calculated pi (MPI+OpenMP): " << pi << std::endl;
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc < 2) {
        if (rank == 0) {
            std::cout << "Usage: mpirun -n <procs> ./program <task> [args]\n";
            std::cout << "Available tasks: 15-32\n";
        }
        MPI_Finalize();
        return 1;
    }

    const int task_num = std::stoi(argv[1]);

    switch (task_num) {
    case 15: task15(); break;
    case 16: task16(); break;
    case 17: task17(); break;
    case 18: task18(); break;
    case 19: task19(); break;
    case 20: task20(); break;
    case 21: task21(); break;
    case 22: task22(); break;
    case 23: task23(); break;
    case 24: {
        const int N = (argc > 2) ? std::stoi(argv[2]) : 1000000000;
        task24(N);
        break;
    }
    case 25: task25(); break;
    case 26: {
        const std::string msg = (argc > 2) ? argv[2] : "A";
        task26(msg);
        break;
    }
    case 27: {
        const int n = (argc > 2) ? std::stoi(argv[2]) : 3;
        task27(n);
        break;
    }
    case 28: {
        const int N = (argc > 2) ? std::stoi(argv[2]) : 1000000000;
        task28(N);
        break;
    }
    case 29: task29(); break;
    case 30: task30(); break;
    case 31: {
        const int num_threads = (argc > 2) ? std::stoi(argv[2]) : 3;
        task31(num_threads);
        break;
    }
    case 32: {
        const int N = (argc > 2) ? std::stoi(argv[2]) : 1000000000;
        task32(N);
        break;
    }
    default:
        if (rank == 0) {
            std::cout << "Unknown task number: " << task_num << std::endl;
        }
    }

    MPI_Finalize();
    return 0;
}