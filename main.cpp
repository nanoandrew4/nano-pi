#include <iostream>
#include <chrono>
#include <thread>

#include <gmp.h>

#define PRECISION 100000
#define ITERATIONS (PRECISION / 14 + 1)

struct sumData {
    int workingThreads = 4;
    mpf_t *parts[5];
    mpf_t pi;
};

void chudnovsky();

void calcA(mpz_t *);

void calcB(mpz_t *);

void calcC(mpz_t *);

void calcD(mpz_t *);

void calcE(mpz_t *);

void sum(sumData *, long initPos);

int main() {

    auto init = std::chrono::high_resolution_clock::now();
    chudnovsky();
    auto final = std::chrono::high_resolution_clock::now();

    std::cout << std::endl << "Calculation took: " << std::chrono::duration_cast<std::chrono::milliseconds>(final - init).count()
              << "ms"
              << std::endl;

    return 0;
}

void chudnovsky() {
    mpf_set_default_prec(PRECISION);

    mpz_t *sets[5];

    for (int i = 0; i < 5; i++)
        sets[i] = new mpz_t[ITERATIONS];

    std::thread threads[4];
    threads[0] = std::thread(calcB, sets[1]);
    threads[1] = std::thread(calcC, sets[2]);
    threads[2] = std::thread(calcD, sets[3]);
    threads[3] = std::thread(calcE, sets[4]);

    calcA(sets[0]);

    for (auto &thread : threads)
        thread.join();

    std::cout << "All values calculated, starting type conversion" << std::endl;

    sumData combinedSet;

    for (int i = 0; i < 5; i++) {
        combinedSet.parts[i] = new mpf_t[ITERATIONS];
        for (int j = 0; j < ITERATIONS; j++) {
            mpf_init(combinedSet.parts[i][j]);
            mpf_set_z(combinedSet.parts[i][j], sets[i][j]);
            mpz_clear(sets[i][j]);
        }
        delete sets[i];
    }

    std::cout << "Finished type conversion, proceeding to perform summation" << std::endl;

    mpf_init(combinedSet.pi);

    std::thread sumThread[3];
    for (int i = 0; i < 3; i++)
        sumThread[i] = std::thread(sum, &combinedSet, (ITERATIONS * i) / (double) combinedSet.workingThreads);

    sum(&combinedSet, (long)(ITERATIONS * (3 / 4.0)));

    for (std::thread &i : sumThread)
        i.join();

    std::cout << "Finished summation" << std::endl;

    mpf_t constant;
    mpf_init_set_ui(constant, 10005);
    mpf_sqrt(constant, constant);
    mpf_mul_ui(constant, constant, 426880);

    mpf_div(combinedSet.pi, constant, combinedSet.pi);

    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < ITERATIONS; j++)
            mpf_clear(combinedSet.parts[i][j]);
        delete combinedSet.parts[i];
    }

    mpf_out_str(NULL, 10, PRECISION, combinedSet.pi);
}

void sum(sumData *sumD, long initPos) {
    long endPos = (ITERATIONS / sumD->workingThreads) + initPos;

    for (long i = initPos; i < endPos; i++) {
        mpf_div(sumD->parts[0][i], sumD->parts[0][i], sumD->parts[2][i]);
        mpf_div(sumD->parts[1][i], sumD->parts[1][i], sumD->parts[3][i]);
        mpf_mul(sumD->parts[0][i], sumD->parts[0][i], sumD->parts[1][i]);
        mpf_div(sumD->parts[0][i], sumD->parts[0][i], sumD->parts[4][i]);
        mpf_add(sumD->pi, sumD->pi, sumD->parts[0][i]);
    }
}

void calcA(mpz_t *lData) {
    mpz_init_set_ui(lData[0], 1);

    for (int i = 1; i < ITERATIONS; i++) {
        mpz_init(lData[i]);
        mpz_mul_ui(lData[i], lData[i - 1], 6u + i);
        mpz_mul_ui(lData[i], lData[i], 5u + i);
        mpz_mul_ui(lData[i], lData[i], 4u + i);
        mpz_mul_ui(lData[i], lData[i], 3u + i);
        mpz_mul_ui(lData[i], lData[i], 2u + i);
        mpz_mul_ui(lData[i], lData[i], 1u + i);
    }

    std::cout << "A calculation finished" << std::endl;
}

void calcB(mpz_t *lData) {
    mpz_t var, con;

    mpz_init(var);
    mpz_init_set_ui(con, 13591409);

    for (int i = 0; i < ITERATIONS; i++) {
        mpz_init(lData[i]);
        mpz_add(lData[i], var, con);
        mpz_add_ui(var, var, 545140134u);
    }

    std::cout << "B calculation finished" << std::endl;
}

void calcC(mpz_t *lData) {
    mpz_init_set_ui(lData[0], 1);

    for (int i = 1; i < ITERATIONS; i++) {
        mpz_init(lData[i]);
        mpz_mul_ui(lData[i], lData[i - 1], 3u + i);
        mpz_mul_ui(lData[i], lData[i], 2u + i);
        mpz_mul_ui(lData[i], lData[i], 1u + i);
    }

    std::cout << "C calculation finished" << std::endl;
}

void calcD(mpz_t *lData) {
    mpz_t k;
    mpz_init_set_ui(k, 1);

    mpz_init_set_ui(lData[0], 1);
    mpz_init_set_ui(lData[1], 1);
    for (ulong i = 2; i < ITERATIONS; i++) {
        mpz_init(lData[i]);
        mpz_mul_ui(k, k, i);
        mpz_mul(lData[i], k, k);
        mpz_mul(lData[i], lData[i], k);
    }

    std::cout << "D calculation finished" << std::endl;
}

void calcE(mpz_t *lData) {

    mpz_init_set_si(lData[0], 1);
    for (int i = 1; i < ITERATIONS; i++) {
        mpz_init(lData[i]);
        mpz_mul_si(lData[i], lData[i - 1], -262537412640768000);
    }

    std::cout << "E calculation finished" << std::endl;
}
