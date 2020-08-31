#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <unistd.h>
#include <omp.h>

//Cantidad total de iteraciones
#define ITERATIONS	10000
//CLSIZE es la cantidad de nodos cercanos que la hormiga tendrá en cuenta
#define CLSIZE		15
//Alpha representa el parámetro de influencia de feromona
#define ALPHA		1
//Beta representa el parámetro de influencia de distancia
#define BETA		2
//TAUO representa los valores de feromonas iniciales en los vértices
#define TAU0		0.00001
//RH0 es la tasa de evaporación de feromonas
#define RHO			0.5
//Q es una constante para actualizar los niveles de feromonas
#define Q			100

/*
* Estructura para almacenar los datos de la instacia TSP
*/
typedef struct {
	char name[100];
	char comment[100];
	char type[100];
	int dimension;
	char wtype[100];
} InstanceData;

/*
* Estructura para guardar el tour y los nodos visitados por una hormiga
*/
typedef struct {
	int* tour;
	int* visited;
} Ant;

float dist(int, int, int, int);
double tourCost(int*);
int hasCandidatesLeft(int, int);
int argMax(int, int);
int NN(int, int);
void moveAntTo(int, int, int);
void updatePheromoneLevel(int*);