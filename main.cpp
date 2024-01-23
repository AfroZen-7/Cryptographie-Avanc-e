//
//  TP7_Shamir's Secret Sharing Scheme
//  

#include <stdio.h>
#include <iostream>
#include <gmp.h>
#include <gmpxx.h>
#include <cstdlib>
#include <ctime>
#include <unistd.h>

#define BUFFER_SIZE 1024 //Pour faire la fonction de génération aléatoire de prime number
#define BITSTRENGTH  14              /* size of prime number (p) in bits */
#define DEBUG true

void generate_prime(mpz_t num, int bit_strength) {
    // Initialiser l'état aléatoire avec le temps
    gmp_randstate_t gmpRandState;
    gmp_randinit_default(gmpRandState);
    gmp_randseed_ui(gmpRandState, static_cast<unsigned long>(time(NULL)));

    // Génère un nombre aléatoire avec le nombre spécifié de bits
    mpz_urandomb(num, gmpRandState, bit_strength);

    // Assurez-vous que le nombre est impair (dernier bit est 1)
    mpz_setbit(num, 0);

    // Trouve le prochain nombre premier à partir de la valeur générée
    mpz_nextprime(num, num);

    // Libérer la mémoire de l'état aléatoire
    gmp_randclear(gmpRandState);
}

void generate_secret(mpz_t s, int bit_strength, mpz_t p) {
    // Initialiser l'état aléatoire avec le temps et le nombre 4444
    gmp_randstate_t gmpRandState;
    gmp_randinit_default(gmpRandState);
    gmp_randseed_ui(gmpRandState, static_cast<unsigned long>(time(NULL) ^ 4444));

    // Génère un nombre aléatoire avec le nombre spécifié de bits et le prime number spécifié
    mpz_urandomm(s, gmpRandState, p);

    // Libérer la mémoire de l'état aléatoire
    gmp_randclear(gmpRandState);
}

/* Main subroutine */
int main()
{
    /* Declare variables */
    int n = 4;          // Numbers of users (max)
    int k = 3;          // Threshold : minimal number of users => secret

    mpz_t p;            // Prime number
    mpz_t S;            // Secret
    mpz_t Sr;           // Reconstruction of the Secret

    mpz_t a1, a2;       // Coefficients of polynom
    mpz_t alpha1,alpha2,alpha3;  // Lagrangian polynomials in zero

    mpz_t x1,x2,x3,x4;  // Login users
    mpz_t y1,y2,y3,y4;  // Shares of users

    /* This function creates the shares computation. The basic algorithm is...
    *+
    *  1. Initialize Prime Number : we work into Z/pZ
    *  2. Initialize Secret Number : S
    *  3. Compute a random polynom of order k-1
    *  4. Shares computation for each users (xi, yi) for i in [1,n]
    *  5. Reconstruct the secret with k users or more
    *
    */
    
    /*
     *  Step 1: Initialize Prime Number : we work into Z/pZ
     */

    mpz_init(p);
    generate_prime(p, BITSTRENGTH);

    if (DEBUG)
    {
        char p_str[1000]; mpz_get_str(p_str,10,p);
        std::cout << "Random Prime 'p' = " << p_str <<  std::endl;
    }

    /*
     *  Step 2: Initialize Secret Number
     */

    // Générer le nombre secret aléatoirement dans Z/pZ
    mpz_init(S);
    generate_secret(S,BITSTRENGTH,p);
    
    if (DEBUG)
    {
        char S_str[1000]; mpz_get_str(S_str,10,S);
        std::cout << "Secret number 'S' = " << S_str <<  std::endl;
    }

    /*
     *  Step 3: Initialize Coefficient of polynom
     */
    mpz_init(a1); mpz_init_set_str(a1, "3", 0);
    mpz_init(a2); mpz_init_set_str(a2, "10", 0);
    
    //TODO: Delete this part and compute the coeffiecients randomly ( warning: inside Z/pZ )
    
    if (DEBUG)
    {
        char a1_str[1000]; mpz_get_str(a1_str,10,a1);
        char a2_str[1000]; mpz_get_str(a2_str,10,a2);
        char S_str[1000];  mpz_get_str(S_str,10,S);
        std::cout << "Polynom 'P(X)' = " << a2_str << "X^2 + " << a1_str << "X + " << S_str << std::endl;
    }
    
    /*
     *  Step 4: Shares computation for each users (xi, yi)
     */
    mpz_init(x1); mpz_init_set_str(x1, "2", 0);
    mpz_init(x2); mpz_init_set_str(x2, "4", 0);
    mpz_init(x3); mpz_init_set_str(x3, "6", 0);
    mpz_init(x4); mpz_init_set_str(x4, "8", 0);

    mpz_init(y1); mpz_init_set_str(y1, "7", 0);
    mpz_init(y2); mpz_init_set_str(y2, "1", 0);
    mpz_init(y3); mpz_init_set_str(y3, "9", 0);
    mpz_init(y4); mpz_init_set_str(y4, "9", 0);
    
    //TODO: Delete this part and compute the shares of all users with public login
    
    if (DEBUG)
    {
        char x1_str[1000]; mpz_get_str(x1_str,10,x1);
        char x2_str[1000]; mpz_get_str(x2_str,10,x2);
        char x3_str[1000]; mpz_get_str(x3_str,10,x3);
        char x4_str[1000]; mpz_get_str(x4_str,10,x4);

        char y1_str[1000]; mpz_get_str(y1_str,10,y1);
        char y2_str[1000]; mpz_get_str(y2_str,10,y2);
        char y3_str[1000]; mpz_get_str(y3_str,10,y3);
        char y4_str[1000]; mpz_get_str(y4_str,10,y4);
        
        std::cout << "Login and share of each users : " << "( x1="<< x1_str << " ; y1=" << y1_str << " ) , "  << "( x2="<< x2_str << " ; y2=" << y2_str << " ) , "  << "( x3="<< x3_str << " ; y3=" << y3_str << " ) , "  << "( x4="<< x4_str << " , y4=" << y4_str << " )" << std::endl;
    }

    /*
     *  Step 5: Sample for reconstruct the secret with 3 users (x1, x2, x3)
     */
    mpz_init(alpha1); mpz_init_set_str(alpha1, "3", 0);
    mpz_init(alpha2); mpz_init_set_str(alpha2, "8", 0);
    mpz_init(alpha3); mpz_init_set_str(alpha3, "1", 0);

    //TODO: Delete this part and automatically compute the secret with k or more shares
    
    // Compute Secret = sum_{i=1}^{k} alpha_i x y_i
    mpz_init(Sr); mpz_init_set_str(Sr, "0", 0);
    mpz_t temp; mpz_init(temp);
    
    mpz_mul(temp,alpha1,y1);
    mpz_add(Sr, Sr, temp);
    mpz_mul(temp,alpha2,y2);
    mpz_add(Sr, Sr, temp);
    mpz_mul(temp,alpha3,y3);
    mpz_add(Sr, Sr, temp);
    mpz_mod(Sr, Sr, p );
    
    if (DEBUG)
    {
        char Sr_str[1000]; mpz_get_str(Sr_str,10,Sr);
        std::cout << "Reconstruction of the secret : S = " << Sr_str << std::endl;
    }
    
    /* Clean up the GMP integers */
    mpz_clear(y1);mpz_clear(y2);mpz_clear(y3);mpz_clear(y4);
    mpz_clear(x1);mpz_clear(x2);mpz_clear(x3);mpz_clear(x4);
    mpz_clear(alpha1);mpz_clear(alpha2);mpz_clear(alpha3);
    mpz_clear(a1);mpz_clear(a2);
    mpz_clear(temp);
    mpz_clear(Sr);
    mpz_clear(S);
    mpz_clear(p);
}

