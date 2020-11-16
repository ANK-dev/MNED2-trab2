/******************************************************************************
 *        Métodos Numéricos para Equações Diferenciais II -- Trabalho 2       *
 *                          Ariel Nogueira Kovaljski                          *
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "main.h"

int main(void)
{
    double Q_new[nx];   /* Array de Q no tempo n+1 */
    double Q_old[nx];   /* Array de Q no tempo n */

    puts("MNED II - Trabalho 2\n"
         "====================\n"
         "por Ariel Nogueira Kovaljski\n");

    listParameters();

    /* Cálculo de Q via método Forward Time-Backward Space (FTBS) */
    initializeArray(Q_old, nx, A, B, C, D, E);
    initializeArray(Q_new, nx, A, B, C, D, E);
    calculateQ_FTBS(Q_old, Q_new);
    printAndSaveResults(Q_new, nx, FTBS);

    /* Cálculo de Q via método Lax-Friedrichs */
    initializeArray(Q_old, nx, A, B, C, D, E);
    initializeArray(Q_new, nx, A, B, C, D, E);
    calculateQ_LF(Q_old, Q_new);
    printAndSaveResults(Q_new, nx, LF);

    /* Cálculo de Q via método Lax-Wendroff */
    initializeArray(Q_old, nx, A, B, C, D, E);
    initializeArray(Q_new, nx, A, B, C, D, E);
    calculateQ_LW(Q_old, Q_new);
    printAndSaveResults(Q_new, nx, LW);

    /* Cálculo de Q via método Beam-Warmimg */
    initializeArray(Q_old, nx, A, B, C, D, E);
    initializeArray(Q_new, nx, A, B, C, D, E);
    calculateQ_BW(Q_old, Q_new);
    printAndSaveResults(Q_new, nx, BW);

    return 0;
}

void listParameters()
{
    puts("Parametros\n----------");
    puts("Constantes da equacao:");
    printf("Delta_t = %f, Delta_x = %f, u_bar = %3.2e\n\n",
           Delta_t, Delta_t, u_bar);
    puts("Constantes da simulacao:");
    printf("nx = %d, t_final = %f\n\n", nx, t_final);
}

/* Inicializa um array para um valor de entrada */
void initializeArray(double arr[], int len, double a, double b, double c,
                     double d, double e)
{
    /*************************************************************
     *
     *                    Condição de contorno
     *                c(x,0) = exp(-A(x-B)²) + s(x)
     * 
     *           onde s(x) = { E,    se  C <= x <= D
     *                       { 0,    c.c.
     */

    int i;
    double x, s;

    for (i = 0; i < len; ++i) {
        x = i * Delta_x;
        s = (x >= c && x <= d ? e : 0);
        arr[i] = exp( -a * ((x - b)*(x - b)) ) + s;
    }
}

/* Calcula as concentrações na malha ao longo do tempo */
/* Método Forward Time-Backward Space (FTBS) */
void calculateQ_FTBS(double old[], double new[])
{
    int i;
    int progress = 0, progress_count = 0;
    int progress_incr = (t_final/Delta_t) * 5.0 / 100;
    double t = 0;

    puts("Calculando FTBS...");

    do {
        /*************************************************************
         *
         *            |==@==|==@==|==@==|==@==|==@==|==@==|
         *               ^
         *            Para o volume da fronteira esquerda
         *        o índice 0 refere-se ao volume nº 1 da malha
         */
        new[0] = old[0];

        /*************************************************************
         *
         *            |==@==|==@==|==@==|==@==|==@==|==@==|
         *                     ^     ^     ^     ^     ^
         *                  Para os volumes do centro 
         *               e na fronteira direita da malha
         */
        for (i = 1; i < nx; ++i) {
            new[i] = old[i] - u_bar*Delta_t/Delta_x * (
                         old[i] - old[i-1]
                     );
        }

        /* Incrementa o progresso a cada 5% */
        if (progress_count == progress_incr){
            progress_count = 0;
            ++progress;
            printf("\rCalculando... %d%% concluido", progress * 5);
            fflush(stdout);
        }

        /* Atualiza array de valores antigos com os novos para a próxima
           iteração */
        for (i = 0; i < nx; ++i) {
            old[i] = new[i];
        }

        /* Incrementa contador para cada 5% */
        ++progress_count;

    } while ( (t += Delta_t) <= t_final );
}

/* Método Lax-Friedrichs */
void calculateQ_LF(double old[], double new[])
{
    int i;
    int progress = 0, progress_count = 0;
    int progress_incr = (t_final/Delta_t) * 5.0 / 100;
    double t = 0;

    puts("Calculando Lax-Friedrichs...");

    do {
        /*************************************************************
         *
         *            |==@==|==@==|==@==|==@==|==@==|==@==|
         *               ^
         *            Para o volume da fronteira esquerda
         *        o índice 0 refere-se ao volume nº 1 da malha
         */
        new[0] = (
                     old[1] + old[0]
                 )/2 - u_bar*Delta_t/(2*Delta_x) * (
                     old[1] - old[0]
                 );

        /*************************************************************
         *
         *            |==@==|==@==|==@==|==@==|==@==|==@==|
         *                     ^     ^     ^     ^
         *             Para os volumes do centro da malha
         */
        for (i = 1; i < nx - 1; ++i) {
            new[i] = (
                         old[i+1] + old[i-1]
                     )/2 - u_bar*Delta_t/(2*Delta_x) * (
                         old[i+1] - old[i-1]
                     );
        }

        /*************************************************************
         *
         *            |==@==|==@==|==@==|==@==|==@==|==@==|
         *                                             ^
         *             Para o volume da fronteira direita
         *            i possui valor de nx - 1 nesse ponto
         */
        new[i] = (
                     old[i] + old[i-1]
                 )/2 - u_bar*Delta_t/(2*Delta_x) * (
                     old[i] - old[i-1]
                 );

        /* Incrementa o progresso a cada 5% */
        if (progress_count == progress_incr){
            progress_count = 0;
            ++progress;
            printf("\rCalculando... %d%% concluido", progress * 5);
            fflush(stdout);
        }

        /* Atualiza array de valores antigos com os novos para a próxima
           iteração */
        for (i = 0; i < nx; ++i) {
            old[i] = new[i];
        }

        /* Incrementa contador para cada 5% */
        ++progress_count;

    } while ( (t += Delta_t) <= t_final );
}

/* Método Lax-Wendroff */
void calculateQ_LW(double old[], double new[])
{
    int i;
    int progress = 0, progress_count = 0;
    int progress_incr = (t_final/Delta_t) * 5.0 / 100;
    double t = 0;

    puts("Calculando Lax-Wendroff...");

    do {
        /*************************************************************
         *
         *            |==@==|==@==|==@==|==@==|==@==|==@==|
         *               ^
         *            Para o volume da fronteira esquerda
         *        o índice 0 refere-se ao volume nº 1 da malha
         */
        new[0] = old[0] - u_bar*Delta_t/(2*Delta_x) * (
                     old[1] - old[0]
                 ) + (u_bar*u_bar)*(Delta_t*Delta_t)/(2*Delta_x*Delta_x) * (
                     old[1] - old[0]
                 );

        /*************************************************************
         *
         *            |==@==|==@==|==@==|==@==|==@==|==@==|
         *                     ^     ^     ^     ^
         *             Para os volumes do centro da malha
         */
        for (i = 1; i < nx - 1; ++i) {
            new[i] = old[i] - u_bar*Delta_t/(2*Delta_x) * (
                         old[i+1] - old[i-1]
                     ) + (u_bar*u_bar)*(Delta_t*Delta_t)/(2*Delta_x*Delta_x) * (
                         old[i+1] - 2*old[i] + old[i-1]
                     );
        }

        /*************************************************************
         *
         *            |==@==|==@==|==@==|==@==|==@==|==@==|
         *                                             ^
         *             Para o volume da fronteira direita
         *            i possui valor de nx - 1 nesse ponto
         */
        new[i] = old[i] - u_bar*Delta_t/(2*Delta_x) * (
                     old[i] - old[i-1]
                 ) + (u_bar*u_bar)*(Delta_t*Delta_t)/(2*Delta_x*Delta_x) * (
                     old[i-1] - old[i]
                 );

        /* Incrementa o progresso a cada 5% */
        if (progress_count == progress_incr){
            progress_count = 0;
            ++progress;
            printf("\rCalculando... %d%% concluido", progress * 5);
            fflush(stdout);
        }

        /* Atualiza array de valores antigos com os novos para a próxima
           iteração */
        for (i = 0; i < nx; ++i) {
            old[i] = new[i];
        }

        /* Incrementa contador para cada 5% */
        ++progress_count;

    } while ( (t += Delta_t) <= t_final );
}

/* Método Beam-Warming */
void calculateQ_BW(double old[], double new[])
{
    int i;
    int progress = 0, progress_count = 0;
    int progress_incr = (t_final/Delta_t) * 5.0 / 100;
    double t = 0;
    
    puts("Calculando Beam-Warming...");

    do {
        /*************************************************************
         *
         *            |==@==|==@==|==@==|==@==|==@==|==@==|
         *               ^
         *            Para o volume da fronteira esquerda
         *        o índice 0 refere-se ao volume nº 1 da malha.
         *     O método de Lax-Wendroff é utilizado para este caso
         */
        new[0] = old[0] - u_bar*Delta_t/(2*Delta_x) * (
                     old[1] - old[0]
                 ) + (u_bar*u_bar)*(Delta_t*Delta_t)/(2*Delta_x*Delta_x) * (
                     old[1] - old[0]
                 );

        /*************************************************************
         *
         *            |==@==|==@==|==@==|==@==|==@==|==@==|
         *                     ^ 
         *         Para o volume logo após a fronteira esquerda
         */
        new[1] = old[1] - 3 * u_bar*Delta_t/(2*Delta_x) * (
                     old[1] - old[0]
                 ) + (u_bar*u_bar)*(Delta_t*Delta_t)/(2*Delta_x*Delta_x) * (
                     old[1] - old[0]
                 );
                 
        /*************************************************************
         *
         *            |==@==|==@==|==@==|==@==|==@==|==@==|
         *                           ^     ^     ^     ^
         *                  Para os volumes do centro 
         *               e na fronteira direita da malha
         */
        for (i = 2; i < nx; ++i) {
            new[i] = old[i] - u_bar*Delta_t/(2*Delta_x) * (
                         3*old[i] - 4*old[i-1] + old[i-2]
                     ) + (u_bar*u_bar)*(Delta_t*Delta_t)/(2*Delta_x*Delta_x) * (
                         old[i] - 2*old[i-1] + old[i-2]
                     );
        }

        /* Incrementa o progresso a cada 5% */
        if (progress_count == progress_incr){
            progress_count = 0;
            ++progress;
            printf("\rCalculando... %d%% concluido", progress * 5);
            fflush(stdout);
        }

        /* Atualiza array de valores antigos com os novos para a próxima
           iteração */
        for (i = 0; i < nx; ++i) {
            old[i] = new[i];
        }

        /* Incrementa contador para cada 5% */
        ++progress_count;

    } while ( (t += Delta_t) <= t_final );
}

/* Imprime na tela e salva os resultados num arquivo de saída */
void printAndSaveResults(double arr[], int len, int method)
{
    int i;
    char filename[50], file0[50], file1[50];
    FILE *results_file;     /* Ponteiro para o arquivo de resultados */

    /* Prepara nome do arquivo de saída */
    switch (method) {
        case FTBS:
            snprintf(filename, 50, "%s%d%s", "results", method, ".txt");
            break;
        case LF:
            snprintf(filename, 50, "%s%d%s", "results", method, ".txt");
            break;
        case LW:
            snprintf(filename, 50, "%s%d%s", "results", method, ".txt");
            break;
        case BW:
            snprintf(filename, 50, "%s%d%s", "results", method, ".txt");
            break;
    }

    /* Imprime os resultados no console */
    printf("\n\nQ[%d] (tempo final: %.2fs) = [", nx, t_final);
    for (i = 0; i < len - 1; ++i) {
        printf("%f, ", arr[i]);
    }
    printf("%f]\n\n", arr[i]);

    /* Error Handling -- Verifica se é possível criar/escrever o arquivo de
                         resultados */
    snprintf(file0, 50, "%s%s", "./results/", filename);
    snprintf(file1, 50, "%s%s", "./../results/", filename);
    if (   ( results_file = fopen(file0, "w") ) == NULL
        && ( results_file = fopen(file1, "w") ) == NULL) {
        fprintf(stderr, "[ERR] Houve um erro ao escrever o arquivo \"%s\"! "
                        "Os resultados nao foram salvos.\n", filename);
        exit(1);
    }

    /* Adiciona os resultados no arquivo "results.txt" */
    fprintf(results_file,
            "nx=%d\n"
            "Delta_t=%f\n"
            "Delta_x=%f\n"
            "t_final=%f\n"
            "u_bar=%f\n",
            nx, Delta_t, Delta_x, t_final, u_bar);
    fputs("********************\n", results_file);

    for (i = 0; i < len; ++i) {
        fprintf(results_file, "%f,%f\n", i * Delta_x, arr[i]);
    }

    fclose(results_file);   /* Fecha o arquivo */

    printf("[INFO] Os resultados foram salvos no arquivo \"%s\" "
           "no diretorio \"results/\".\n\n", filename);
}