/*
 * dixon.c
 *
 *   	Date: June 8, 2017
 *      Author: Mert Yassı
 */
// DIXON's ALGORITHM
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"dixon.h"

void dixon_fact_method(ull_int n, int fbb) {
	arrays_t a;
	a = sieveOfEra(fbb);
	ull_int k = ceil(sqrt(n)), size = a.numOfPrimes + 3, r, i, j, li = 0, x = 1, y2 = 1, y, xpy, xmy, p ,q, old_x, old_y;
	ull_int **arr, *lArr, *kArr;
	arr = (ull_int**)malloc(sizeof(ull_int*) * a.numOfPrimes);

	for(i = 0; i < a.numOfPrimes; i++){
		arr[i] = (ull_int*)malloc( sizeof(ull_int) * size);
	}

	a.exponents = (ull_int*)malloc(sizeof(ull_int) * a.numOfPrimes);
	lArr = (ull_int*)malloc(sizeof(ull_int) * size);
	kArr = (ull_int*)malloc(sizeof(ull_int) * size);

	printf("Factorbase primes: ");
	for(i = 0; i < a.numOfPrimes; i++){
		printf("%llu ",a.factorBase[i]);
	}
	printf("\n");

	for(r = 0; r < size;){
		for(i = 0; i < a.numOfPrimes; i++){
			a.exponents[i] = 0; // make all exponents initially 0
		}
		ull_int l = (k * k) % n;
		ull_int ktemp = k;
		ull_int ltemp = l;
		int control = 0; i = 0;
		while(control == 0){
			if(l % a.factorBase[i] == 0 && l != 0){
				while(l % a.factorBase[i] == 0 && l != 0){
					l = l / a.factorBase[i];
					a.exponents[i]++;
				}
			}else{
				if(i < a.numOfPrimes - 1){
					i++;
				}else{
					control = 1;
				}
			}
		}
		if(l == 1){
			kArr[li] = ktemp;
			lArr[li] = ltemp;
			li++;
			printf("[ ");
			for(i=0;i<a.numOfPrimes;i++){
				arr[i][r] = a.exponents[i] % 2;
				printf("%llu ",arr[i][r]);
			}
			printf("]\n");
			r++;
		}
		k++;
	}

	FILE *f = fopen("matrix.txt", "w");
	if (f == NULL){
		printf("Error opening file!\n");
		exit(1);
	}
	ull_int count = 0, count2 = 0;
	fprintf(f,"M:=Matrix(GF(2),%llu,%llu,[", size, a.numOfPrimes);
	for(j = 0; j < r; j++){
		count++;
		for(i = 0; i < a.numOfPrimes; i++){
			count2++;
			fprintf(f,"%llu",arr[i][j]);
			if(count2 < a.numOfPrimes){
				fprintf(f,",");
			}
		}
		if(count < r){
			fprintf(f,",");
		}
		count2 = 0;
	}
	fprintf(f,"]);");
	fprintf(f,"\nN:=Nullspace(M);");
	fprintf(f,"\nB:=Basis(N);");
	fprintf(f,"\nB;");
	fclose(f);
	//int ret = system("~/magma < matrix.txt > result.txt"); // you can factor larger numbers using local Magma
	int ret = system("python3 encoder.py matrix.txt > result.txt"); // will not work if the size of matrix.txt > 50 kB (Online Magma input limit)

	FILE *f2 = fopen("result.txt", "rb");
	fseek(f2, 0, SEEK_END);
	long fsize = ftell(f2);
	fseek(f2, 0, SEEK_SET);
	char *string = malloc(fsize + 1);
	ret = fread(string, fsize, 1, f2);
	fclose(f2);
	string[fsize] = 0;

	int sct = 1; // sct = 0 when using local Magma
	i = 0;
	while(string[i] != '\0'){
		if(string[i] == 44){
			sct++;
		}
		i++;
	}

	char **arr2;
	arr2 = (char**)malloc(sizeof(char*) * sct);
	for(i = 0; i < sct; i++){
		arr2[i] = (char*)malloc(sizeof(char) * size);
	}

	i = 0, j = 0, k = 0;
	while(string[i] != '\0'){
		if(string[i] == '('){
			while(string[i] != ')'){
				if(((string[i] == 48) || (string[i] == 49))){
					arr2[k][j] = string[i];
					j++;
				}
				i++;
			}
			j = 0;
			k++;
		}else{
			i++;
		}
	}

	if(sct == 1){
		printf("Magma could not find NullSpace of M!\n");
	}else{
		printf("NullSpace:\n");
		for(i = 0; i < sct; i++){
			printf(" [ ");
			for(j = 0; j < size; j++){
				printf("%c ",arr2[i][j]);
			}
			printf("]\n");
		}
	}
	int ctrl = 1;
	for(i = 0; i < sct; i++){
		for(j = 0; j < size; j++){
			if(arr2[i][j] == 49){
				old_x = kArr[j];
				old_y = lArr[j];
				x = old_x * x;
				y2 = old_y * y2;
			}
		}
		y = sqrt(y2);
		xpy = x + y;
		xmy = x - y;
		p = gcd(xpy, n);
		q = gcd(xmy, n);

		if((p != 1 && p != n) && (q != 1 && q != n)){
			printf("1ST FACTOR P: %llu\n",p);
			printf("2ND FACTOR Q: %llu\n",q);
			ctrl = 1;
			i = sct;
		}else{
			x = 1;
			y2 = 1;
		}
	}
	if(ctrl != 1){
		printf("The algorithm could not find proper factors!\n");
	}

	free(a.factorBase);
	free(a.exponents);
	free(kArr);
	free(lArr);
	for(i = 0; i < a.numOfPrimes; i++){
		free(arr[i]);
	}
	free(arr);
	for(i = 0; i < sct; i++){
		free(arr2[i]);
	}
	free(arr2);
}

arrays_t sieveOfEra(int n) {
	ull_int *primes, *factorArr, i, j, count = 0, k = 0;
	arrays_t a;

	primes = (ull_int*)malloc(sizeof(ull_int) * n);

	for(i = 2; i < n; i++){
		primes[i] = 1; // make all elements true
	}

	for(i = 2; i < n; i++){
		if(primes[i] == 1){
			for(j = i; i * j < n; j++){
				primes[i * j] = 0; // mark composites
			}
		}
	}

	for(i = 2; i < n; i++){ // loop for determining numOfPrimes
		if(primes[i] == 1){
			count++;
		}
	}

	factorArr = (ull_int*)malloc(sizeof(ull_int) * count);
	for(i = 2; i < n; i++){
		if(primes[i] == 1){
			factorArr[k] = i;
			k++;
		}
	}
	a.factorBase = factorArr;
	a.numOfPrimes = count;
	free(primes);

	return a;
}

ull_int gcd(ull_int a, ull_int b) {
	ull_int i, j, temp;
	if(a > b){
		i = a;
		j = b;
	}
	else{
		i = b;
		j = a;
	}
	while(j != 0) {
		temp = i % j;
		i = j;
		j = temp;
	}
	return i;
}

