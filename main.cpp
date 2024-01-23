#include <iostream>
#include <gmp.h>
#include <vector>

#define BITSTRENGTH 14
#define DEBUG true

// Function to generate a random prime number
void generate_prime(mpz_t num, int bit_strength, gmp_randstate_t gmpRandState) 
{
    mpz_urandomb(num, gmpRandState, bit_strength);
    mpz_setbit(num, 0);
    mpz_nextprime(num, num);
}

// Function to generate a random secret within the range [0, prime)
void generate_secret(mpz_t secret, mpz_t prime, gmp_randstate_t gmpRandState) 
{
    mpz_urandomm(secret, gmpRandState, prime);
}

// Function to generate coefficients for the polynomial and set the last coefficient as the secret
void generate_coefficients(mpz_t* coefficients, mpz_t prime, int& change, mpz_t secret, gmp_randstate_t gmpRandState) 
{
    // Generate random coefficients for the polynomial
    for (int i = 0; i < change - 1; i++) {
        mpz_urandomm(coefficients[i], gmpRandState, prime);
    }

    // Set the last coefficient as the secret
    mpz_set(coefficients[change - 1], secret);
}

// Function to compute shares for a given user
void compute_shares(mpz_t x, mpz_t y, mpz_t* coefficients, int change) 
{
    // Initialize y to 0
    mpz_init_set_ui(y, 0);

    mpz_t temp;
    mpz_init(temp);

    // Compute the polynomial expression for the given x
    for (int j = 0; j < change; j++) {
        mpz_pow_ui(temp, x, j);     // x^j
        mpz_mul(temp, temp, coefficients[j]);  // coefficients[j] * x^j
        mpz_add(y, y, temp);        // Sum up the terms
    }

    mpz_clear(temp);
}

// Function to compute Lagrange coefficients for interpolation
void compute_lagrange_coefficients(std::vector<mpz_t>& alphas, mpz_t* x, int k, mpz_t prime) 
{
    // Compute the Lagrange coefficients for interpolation
    for (int i = 0; i < k; i++) {
        mpz_init(alphas[i]);
        mpz_init_set_ui(alphas[i], 1);

        for (int j = 0; j < k; j++) {
            if (j != i) {
                mpz_t temp;
                mpz_init(temp);

                // Compute (x[j] - x[i])^-1 mod prime
                mpz_sub(temp, x[j], x[i]);
                mpz_invert(temp, temp, prime);

                // Update the Lagrange coefficient
                mpz_mul(alphas[i], alphas[i], temp);
                mpz_clear(temp);
            }
        }
    }
}

// Function to reconstruct the secret using Lagrange coefficients and shares
void reconstruct_secret(mpz_t reconstructedSecret, std::vector<mpz_t>& alphas, mpz_t* shares, int k, mpz_t p) 
{
    // Initialize the reconstructed secret to 0
    mpz_init_set_ui(reconstructedSecret, 0);

    mpz_t temp;
    mpz_init(temp);

    // Reconstruct the secret using Lagrange interpolation
    for (int i = 0; i < k; i++) {
        mpz_mul(temp, alphas[i], shares[i]);
        mpz_add(reconstructedSecret, reconstructedSecret, temp);
    }

    // Take the result modulo p
    mpz_mod(reconstructedSecret, reconstructedSecret, p);
    mpz_clear(temp);
}

int main() 
{
    int n = 4;  // Numbers of users (max)
    int k = 3;  // Threshold: minimal number of users => secret

    mpz_t p;            // Prime number
    mpz_t S;            // Secret
    mpz_t Sr;           // Reconstruction of the Secret

    std::vector<mpz_t> a(k);       // Coefficients of polynomial
    std::vector<mpz_t> alphas(k);  // Lagrangian polynomials in zero

    std::vector<mpz_t> x(n);  // Login users
    std::vector<mpz_t> y(n);  // Shares of users

    // Initialize GMP random state
    gmp_randstate_t gmpRandState;
    gmp_randinit_default(gmpRandState);
    gmp_randseed_ui(gmpRandState, static_cast<unsigned long>(time(NULL)));

    /* 
     * This function creates the shares computation. The basic algorithm is...
     * 1. Initialize Prime Number: we work into Z/pZ
     * 2. Initialize Secret Number: S
     * 3. Compute a random polynomial of order k-1
     * 4. Shares computation for each user (xi, yi) for i in [1,n]
     * 5. Reconstruct the secret with k users or more
     */

    /*
     * Step 1: Initialize Prime Number: we work into Z/pZ
     */

    mpz_init(p);
    generate_prime(p, BITSTRENGTH, gmpRandState);

    if (DEBUG) 
    {
        char p_str[1000]; mpz_get_str(p_str, 10, p);
        std::cout << "Random Prime 'p' = " << p_str << std::endl;
    }

    /*
     * Step 2: Initialize Secret Number
     */

    mpz_init(S);
    generate_secret(S, p, gmpRandState);

    if (DEBUG) 
    {
        char S_str[1000]; mpz_get_str(S_str, 10, S);
        std::cout << "Secret number 'S' = " << S_str << std::endl;
    }

    /*
     * Step 3: Initialize Coefficient of polynomial
     */

    // Initialize coefficients of the polynomial
    for (int i = 0; i < k - 1; i++) {
        mpz_init(a[i]);
    }
    mpz_init_set(a[k - 1], S);

    if (DEBUG) 
    {
        char a1_str[1000]; mpz_get_str(a1_str, 10, a[0]);
        char a2_str[1000]; mpz_get_str(a2_str, 10, a[1]);
        char S_str[1000];  mpz_get_str(S_str, 10, S);
        std::cout << "Polynom 'P(X)' = " << a2_str << "X^2 + " << a1_str << "X + " << S_str << std::endl;
    }

    /*
     * Step 4: Shares computation for each user (xi, yi)
     */

    // Initialize login users
    for (int i = 0; i < n; i++) {
        mpz_init(x[i]);
        mpz_set_ui(x[i], i + 1);
    }

    // Compute shares for each user
    for (int i = 0; i < n; i++) {
        mpz_init(y[i]);
        compute_shares(x[i], y[i], a.data(), k);
    }

    if (DEBUG) 
    {
        char x1_str[1000]; mpz_get_str(x1_str, 10, x[0]);
        char x2_str[1000]; mpz_get_str(x2_str, 10, x[1]);
        char x3_str[1000]; mpz_get_str(x3_str, 10, x[2]);
        char x4_str[1000]; mpz_get_str(x4_str, 10, x[3]);

        char y1_str[1000]; mpz_get_str(y1_str, 10, y[0]);
        char y2_str[1000]; mpz_get_str(y2_str, 10, y[1]);
        char y3_str[1000]; mpz_get_str(y3_str, 10, y[2]);
        char y4_str[1000]; mpz_get_str(y4_str, 10, y[3]);

        std::cout << "Login and share of each users : " << "( x1="<< x1_str << " ; y1=" << y1_str << " ) , "  << "( x2="<< x2_str << " ; y2=" << y2_str << " ) , "  << "( x3="<< x3_str << " ; y3=" << y3_str << " ) , "  << "( x4="<< x4_str << " , y4=" << y4_str << " )" << std::endl;
    }

    /*
     * Step 5: Sample for reconstruct the secret with 3 users (x1, x2, x3)
     */

    // Compute Lagrange coefficients
    compute_lagrange_coefficients(alphas, x.data(), k, p);

    // Reconstruct the secret
    reconstruct_secret(Sr, alphas, y.data(), k, p);

    if (DEBUG) 
    {
        char Sr_str[1000];
        mpz_get_str(Sr_str, 10, Sr);
        std::cout << "Reconstruction of the secret : S = " << Sr_str << std::endl;
    }

    // Clean up the GMP integers
    mpz_clear(S);
    mpz_clear(p);
    mpz_clear(Sr);

    for (int i = 0; i < n; i++) {
        mpz_clear(x[i]);
        mpz_clear(y[i]);
    }

    for (int i = 0; i < k; i++) {
        mpz_clear(a[i]);
        mpz_clear(alphas[i]);
    }

    gmp_randclear(gmpRandState);

    return 0;
}