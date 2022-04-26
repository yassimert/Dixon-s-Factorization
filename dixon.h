/*
 * dixon.h
 *
 *  Created on: June 8, 2017
 *      Author: Mert Yassı
 */
// DIXON's ALGORITHM

#ifndef DIXON_H_
#define DIXON_H_

typedef unsigned long long int ull_int;

typedef struct arrays
{
	ull_int *factorBase, *exponents;
	ull_int numOfPrimes;
}arrays_t;

void dixon_fact_method(ull_int n, int fbb);
arrays_t sieveOfEra(int n);
ull_int gcd(ull_int a, ull_int b);

#endif /* DIXON_H_ */
