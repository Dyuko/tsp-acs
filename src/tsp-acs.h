#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <unistd.h>
#include <pthread.h>
#include <semaphore.h>

#define ITERATIONS	10000
#define CLSIZE		15
#define Q0			0.9
#define ALPHA		1
#define BETA		2
#define TAU0		0.00001
#define RHO			0.5
#define Q			100

typedef struct {
	char name[100];
	char comment[100];
	char type[100];
	int dimension;
	char wtype[100];
} InstanceData;

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
void *routeAnt(void *);