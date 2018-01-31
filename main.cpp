#include <iostream>
#include <chrono>
#include <thread>

#include <gmp.h>

#define PRECISION 1000
#define ITERATIONS (PRECISION / 14 + 1)

struct sumData {
    int initPos;
    int workingThreads = 4;
    mpf_t *parts[5];
    mpf_t sum;
};

void opt();

void calcA(mpz_t *);

void calcB(mpz_t *);

void calcC(mpz_t *);

void calcD(mpz_t *);

void calcE(mpz_t *);

void sum(sumData *);

void prev();

int main() {

    auto init = std::chrono::high_resolution_clock::now();
    opt();
    auto final = std::chrono::high_resolution_clock::now();

    std::cout << "Calculation took: " << std::chrono::duration_cast<std::chrono::milliseconds>(final - init).count()
              << "ms"
              << std::endl;

    return 0;
}

void opt() {
    mpf_set_default_prec(PRECISION);

    mpz_t **sets = new mpz_t *[5];

    for (int i = 0; i < 5; i++)
        sets[i] = new mpz_t[ITERATIONS];

    std::thread threads[4];
    threads[0] = std::thread(calcB, sets[1]);
    threads[1] = std::thread(calcC, sets[2]);
    threads[2] = std::thread(calcD, sets[3]);
    threads[3] = std::thread(calcE, sets[4]);

    calcA(sets[0]);

    for (int i = 0; i < 4; i++)
        threads[i].join();

    std::cout << "All values calculated, starting type conversion and combination" << std::endl;

    mpz_out_str(NULL, 10, sets[0][1]);

    sumData combinedSet;

    mpf_init(combinedSet.sum);

    for (int i = 0; i < 5; i++) {
        combinedSet.parts[i] = new mpf_t[ITERATIONS];
        for (int j = 0; j < ITERATIONS; j++)
            mpf_set_z(combinedSet.parts[i][j], sets[i][j]);
        mpz_clears(*sets[i]);
        delete sets[i];
    }

    std::cout << "Finished type conversion and cleanup, proceeding to perform combinations" << std::endl;

    sum(&combinedSet);

//    pthread_t sumThreads[combinedSet.workingThreads];

    for (int i = 0; i < 5; i++) {
        mpf_clears(*combinedSet.parts[i]);
        delete combinedSet.parts[i];
    }

    mpf_out_str(NULL, 10, PRECISION, combinedSet.sum);
}

void sum(sumData *sumD) {
//    sumData* sumD = (sumData *) arg;
    int initPos = sumD->initPos;
    int endPos = (ITERATIONS / sumD->workingThreads) + initPos;

    for (int i = initPos; i < endPos; i++) {
        mpf_div(sumD->parts[0][i], sumD->parts[0][i], sumD->parts[2][i]);
        mpf_div(sumD->parts[1][i], sumD->parts[1][i], sumD->parts[3][i]);
        mpf_mul(sumD->parts[0][i], sumD->parts[0][i], sumD->parts[1][i]);
        mpf_div(sumD->parts[0][i], sumD->parts[0][i], sumD->parts[4][i]);
        mpf_add(sumD->sum, sumD->sum, sumD->parts[0][i]);
    }

    std::cout << "Finished summation" << std::endl;
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

void prev() {
    mpf_set_default_prec(PRECISION);

    mpf_t pi, sum, top, bottom, Z;
    mpf_inits(pi, sum, top, bottom, Z, 0);
    mpf_sqrt_ui(Z, 10005);
    mpf_mul_ui(Z, Z, 426880);

    mpz_t A, B, C, D, E;
    mpz_inits(A, B, C, D, E, 0);

    for (ulong i = 0; i < ITERATIONS; i++) {
        mpz_fac_ui(A, 6 * i);
        mpz_set_ui(B, 545140134 * i + 13591409);
        mpz_fac_ui(C, 3 * i);
        mpz_fac_ui(D, i);
        mpz_pow_ui(D, D, 3);
        mpz_set_si(E, -262537412640768000);
        mpz_pow_ui(E, E, i);

        mpz_mul(A, A, B);
        mpf_set_z(top, A);

        mpz_mul(C, C, D);
        mpz_mul(C, C, E);
        mpf_set_z(bottom, C);

        mpf_div(top, top, bottom);
        mpf_add(sum, sum, top);
    }

    mpf_div(pi, Z, sum);

    mpf_out_str(NULL, 10, PRECISION, pi);
    std::cout << std::endl;

    mpf_clears(pi, sum, top, bottom, Z, nullptr);
    mpz_clears(A, B, C, D, E, nullptr);
}