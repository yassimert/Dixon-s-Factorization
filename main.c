/*
 * main.c

 *
 *  	Date: June 8, 2017
 *      Author: Mert Yassı
 */
// DIXON's ALGORITHM

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"dixon.h"

#define FBBOUND 200

int main() {
	printf("DIXON's ALGORITHM\n\n");
	ull_int n = 721516342875857;
	printf("n = %llu\n",n);
	dixon_fact_method(n, FBBOUND);
	return 0;
}

