//
//  main.c
//  lab4
//
//  Created by barry on 27.01.2018.
//  Copyright © 2018 barry. All rights reserved.
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define N 5
#define M 6
#define EPS 0.0001
void swap(double*a, double *b);
void method_Gauss(void);
void method_Gauss_Zeidel(void);
void read_data(double [N][M]);
int main(void) {
    method_Gauss();
    method_Gauss_Zeidel();
    return 0;
}
void method_Gauss(){
    double a[N][M] = {{0}}, temp;
    int buff, i, j;
    read_data(a);
    int k = 0; // текущая операция
    
    while (k < N - 1){ // ищем строку с максимальным коэф. при x1
        buff = k;
            for (i = k; i < N; i++) {
                if (a[buff][k] < a[i][k]){
                    buff = i;
                }
            }
        printf("Выбираем строку с максимальным коэффициентом и меняем ее с первой.\n");
            for (j = k; j < M; j++) { // меняем местами строку с элементом max со строкой "k"
                temp = a[buff][j];
                a[buff][j] = a[k][j];
                a[k][j] = temp;
            }
        for (i = 0; i < N; i++) {
            printf("\n");
            for (j = 0; j < M; j++) {
                printf("%lf ", a[i][j]);
            }
        }
        printf("\n\n");
        printf("ормируем уравнения относительно коэффициента");
        for (i = k; i < N; i++) { // нормируем уравн относительно k-того эл-та
            temp = a[i][k];
            for (j = k; j < M; j++) {
                if (!a[i][j]) continue;
                a[i][j] /= temp;
            }
        }
        for (i = 0; i < N; i++) {
            printf("\n");
            for (j = 0; j < M; j++) {
                printf("%lf ", a[i][j]);
            }
        }
        printf("\n\n");
        printf("Вычитаем");
        for (i = k + 1; i < N; i++) { // вычитаем k-тую строку из ставшихся строк
            for (j = k; j < M; j++) {
                a[i][j] -= a[k][j];
            }
        }
        for (i = 0; i < N; i++) {
            printf("\n");
            for (j = 0; j < M; j++) {
                printf("%lf ", a[i][j]);
            }
        }
        printf("\n\n");
        k++;
    }
    printf("последние нормирование");
    temp = a[N-1][M-2];
    a[N-1][M-2] /= temp; //последние нормирование
    a[N-1][M-1] /= temp;
    for (i = 0; i < N; i++) {
        printf("\n");
        for (j = 0; j < M; j++) {
            printf("%lf ", a[i][j]);
        }
    }
    printf("\n\n");
    double ot[N] = {0}, x[N] = {0};
    x[N-1] = a[N-1][M-1];
    for (i = N - 1; i >= 0; i--) { // обратная подстановка
        ot[i] = a[i][M-1];
        for (j = M - 2; j>=0; j--) {
            if (j == i) continue;
            ot[i] -= a[i][j] * x[j];
        }
        x[i] = ot[i];
    }
    printf("Ответы:\n");
    for (i = 0 ; i < N; i++) {
        printf("%lf ", x[i]);
    }
    FILE *file = fopen("answear_data.txt", "w");
    char c = 'a';
    for (i = 0 ; i < N; i++) { // Формат вывода должен быть:«a = ... b = ... c = ... d = ... e = ...»
        fprintf(file, " %c = %lf", c+i, x[i]);
    }
    fprintf(file, "\n");
    double er[N] = {0};
    printf("\n\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            er[i] += a[i][j] * x[j];
        }
        er[i] -= a[i][M-1];
        printf("%lf/ ", er[i]);
    }
    fclose(file);
}

void method_Gauss_Zeidel(){
    double a[N][M] = {{0}};
    int count = N;
    read_data(a);
    while (count==N){ // проверка на нули в главной диагонали и исправление этого
        count = 0;
        for (int i = 0; i < N; i++) {
            if (a[i][i] == 0){
                count++;
            }
        }
        if (count == 5) {
            for (int i = rand()%5; i < N; i++) {
                for (int j = 0; j < M; j++) {
                    swap(&a[0][j], &a[i][j]);
                }
            }
        }
    }
    double old_x[N], new_x[N] = {0};
    //
    double max;
    do{
        for (int i = 0; i < N; i++) { //переприсваивание векторов
            old_x[i] = new_x[i];
        }
        for (int i = 0; i < N; i++) {
            double var = 0;
            for (int j = 0; j<i; j++) {
                var += a[i][j] * new_x[j];
            }
            for (int j = i+ 1; j < N; j++) {
                var += a[i][j] * old_x[j];
            }
            new_x[i] = (a[i][M - 1] - var) / a[i][i];
        }
        double error[N] = {0};
        for (int i = 0; i<N; i++) { //найдем массив погрешностей
            error[i] = fabs((new_x[i] - old_x[i])/new_x[i]);
        }
        max = error[0];
        for (int i = 0; i<N; i++) { //найдем максимальное значение погрешности
            if (error[i] > max) {
                max = error[i];
            }
        }
    }while (max > EPS);
    //
    FILE *file = fopen("answear_data.txt", "a");
    char c = 'a';
    for (int i = 0 ; i < N; i++) { // Формат вывода должен быть:«a = ... b = ... c = ... d = ... e = ...»
        fprintf(file, " %c = %lf", c+i, new_x[i]);
    }
    fprintf(file, "\n");
    double er[N] = {0};
    printf("\n\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            er[i] += a[i][j] * new_x[j];
        }
        er[i] -= a[i][M-1];
        printf("%lf\\ ", er[i]);
    }
}
void read_data(double a[N][M]){
    FILE *file = fopen("input_data.txt", "r");
    int i, j;
    for (i = 0; i< N; i++) {
        printf("\n");
        for (j = 0; j < M; j++) {
            fscanf(file, "%lf", &a[i][j]);
        }
    }
    fclose(file);
}
void swap(double*a, double *b){
    double temp = *a;
    *a = *b;
    *b = temp;
}
