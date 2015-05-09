#include "tsp-acs.h"

int M;

double** distanceMatrix;
double** tau;
int** cl;
InstanceData ins;
Ant ant[1000];

int runParallel;

int main(int argc, char *argv[]) {
	float bestKnownCost = FLT_MAX;
	int* bestKnownTour;

	FILE *file;
	int i, j, k;
	file = fopen(argv[1], "r");
	M = atoi(argv[2]);
	runParallel = atoi(argv[3]);
	fscanf(file, "NAME : %[^\n]s", ins.name);
	fscanf(file, "\nCOMMENT : %[^\n]s", ins.comment);
	fscanf(file, "\nTYPE : %[^\n]s", ins.type);
	fscanf(file, "\nDIMENSION : %d", &ins.dimension);
	fscanf(file, "\nEDGE_WEIGHT_TYPE : %[^\n]s", ins.wtype);
	fscanf(file, "\nNODE_COORD_SECTION");
	if (strcmp(ins.wtype,"EUC_2D")) {
		return 1;
	}
	if (runParallel) {
		printf("Running parallel version\n");
	}
	else {
		printf("Running standard version\n");
	}
	int** coord;
	coord = (int**) malloc(ins.dimension * sizeof(int*));
	bestKnownTour = (int*) malloc(ins.dimension * sizeof(int));
	for(i = 0; i <ins.dimension; i++) {
		coord[i] = (int*) malloc(2 * sizeof(int));
	}
	for (i = 0; i < ins.dimension; i++) {
		fscanf(file, "\n %*[^ ] %d %d", &coord[i][0], &coord[i][1]);
	}
	fclose(file);
	
	distanceMatrix = (double**) malloc(ins.dimension * sizeof(double*));
	for (i = 0; i < ins.dimension; i++) {
		distanceMatrix[i] = (double*) malloc(ins.dimension * sizeof(double));
	}
	
	for (i = 0; i < ins.dimension; i++){
		for ( j = i + 1; j < ins.dimension ; j++) {
			distanceMatrix[i][j] = dist(coord[i][0], coord[i][1], coord[j][0], coord[j][1]);
			distanceMatrix[j][i] = distanceMatrix[i][j];
		}
	}
	free(coord);

	tau = (double**) malloc(ins.dimension*sizeof(double*));
	for (i = 0; i < ins.dimension; i++) {
		tau[i] = (double*) malloc(ins.dimension * sizeof(double));
	}
	for (i = 0; i < ins.dimension; i++) {
		for( j = 0; j < ins.dimension; j++) {
			tau[i][j] = TAU0;
		}
	}
	
	srand(time(NULL));
	for (i = 0; i < M; i++) {
		ant[i].tour = (int*) malloc(ins.dimension * sizeof(int));
		ant[i].tour[0] = rand() % ins.dimension;
		ant[i].visited = (int*) malloc(ins.dimension * sizeof(int));
		for (j = 1; j < ins.dimension; j++) {
			ant[i].visited[j] = 0;
		}
		ant[i].visited[ant[i].tour[0]] = 1;
	}
	
	cl = (int**) malloc(ins.dimension * sizeof(int*));
	int* visited = (int*) malloc(ins.dimension * sizeof(int));
	double nearestDistance;
	int nearestNeighbor;
	for (i = 0; i < ins.dimension; i++) {
		cl[i] = (int*) malloc(CLSIZE * sizeof(int));
	}
	for (i = 0; i < ins.dimension; i++) {
		for(j = 0; j < ins.dimension; j++) {
			visited[j] = 0;
		}
		for(j = 0; j < CLSIZE; j++) {
			nearestDistance = DBL_MAX;
			nearestNeighbor = 0;
			for(k = 0; k < ins.dimension; k++) {
				if(!visited[k] && distanceMatrix[i][k] < nearestDistance && i != k) {
					nearestDistance = distanceMatrix[i][k];
					nearestNeighbor = k;
				}
			}
			visited[nearestNeighbor] = 1;
			cl[i][j] = nearestNeighbor;
		}
	}
	
	pthread_t *thread;
	thread = (pthread_t*) malloc(M * sizeof(pthread_t));
	
	for (i = 0; i < ITERATIONS; i++) {
		if (runParallel) {
			for(j = 0; j < M; j++) {
				int st = pthread_create(&thread[j], NULL, routeAnt, (void *) j);
				if (st) {
					printf("pthread_create() error %d\n",st);
					exit(-1);
				}
			}
		
			for (j = 0; j < M; j++) {
				pthread_join(thread[j],NULL);
			}
			
		}
		else {
			for (j = 0; j < M; j++) {
				for (k = 1; k < ins.dimension; k++) {
					if (hasCandidatesLeft(j,ant[j].tour[k-1])) {
						moveAntTo(j, k, argMax(j, ant[j].tour[k-1]));
					}
					else {
						moveAntTo(j, k, NN(j,k-1));
					}
				}
			}
		}
		
		int bestAnt = 0;
		double lowerCost = DBL_MAX;
		for (j = 0; j < M; j++) {
			if (tourCost(ant[j].tour) < lowerCost) {
				bestAnt = j;
				lowerCost = tourCost(ant[j].tour);
			}
		}

		updatePheromoneLevel(ant[bestAnt].tour);
		
		if (lowerCost < bestKnownCost) {
			bestKnownCost = lowerCost;
			for (k = 0; k < ins.dimension; k++) {
				bestKnownTour[k] = ant[bestAnt].tour[k];
			}
		}

		for (j = 0; j < M; j++) {
			for (k = 0; k < ins.dimension; k++) {
				ant[j].visited[k] = 0;
			}
			ant[j].tour[0] = rand() % ins.dimension;
			ant[j].visited[ant[j].tour[0]] = 1;
			for (k = 1; k < ins.dimension; k++) {
				ant[j].tour[k] = 0;
			}
		}
	}
	
	printf("Best tour: ");
	for (j = 0; j < ins.dimension; j++) {
		printf("%d ", bestKnownTour[j]);
	}
	printf("; Length = %f\n", bestKnownCost);

	return 0;
}

void *routeAnt(void *id) {
	int j, i = (int)id;
	for(j=1;j<ins.dimension;j++){
		if(hasCandidatesLeft(i,ant[i].tour[j-1]))
			moveAntTo(i,j,argMax(i,ant[i].tour[j-1]));
		else
			moveAntTo(i,j,NN(i,j-1));
	}
}

void updatePheromoneLevel(int* tour) {
	int i;
	for (i = 1; i < ins.dimension; i++) {
		tau[tour[i-1]][tour[i]] = (1-RHO) * tau[tour[i-1]][tour[i]] + RHO * (Q/tourCost(tour));
		tau[tour[i]][tour[i-1]] = (1-RHO) * tau[tour[i]][tour[i-1]] + RHO * (Q/tourCost(tour));
	}
}

int NN(int antIndex, int origin) {
	double nearestDistance = DBL_MAX;
	int nearestNeighbor = 0;
	int i;
	for (i = 0; i < ins.dimension; i++) {
		if(!ant[antIndex].visited[i] && distanceMatrix[origin][i] < nearestDistance && origin != i) {
			nearestDistance = distanceMatrix[origin][i];
			nearestNeighbor = i;
		}
	}
	return nearestNeighbor;
}

void moveAntTo(int antIndex, int slot, int city){
	ant[antIndex].tour[slot] = city;
	ant[antIndex].visited[city] = 1;
	tau[ant[antIndex].tour[slot]][ant[antIndex].tour[slot-1]] = (1-RHO) * tau[ant[antIndex].tour[slot]][ant[antIndex].tour[slot-1]] + RHO * TAU0;
	tau[ant[antIndex].tour[slot-1]][ant[antIndex].tour[slot]] = (1-RHO) * tau[ant[antIndex].tour[slot-1]][ant[antIndex].tour[slot]] + RHO * TAU0;
}

int argMax(int antIndex, int currentCity){
	int i, bestCity;
	float maxu, Tau, Nu;
	maxu = (float) - FLT_MAX;
	bestCity = 0;
	for(i = 0; i < CLSIZE; i++) {
		if(!ant[antIndex].visited[cl[currentCity][i]]){
			Tau = pow(tau[currentCity][cl[currentCity][i]], ALPHA);
			Nu = pow((1/distanceMatrix[currentCity][cl[currentCity][i]]), BETA);
			if (Tau * Nu > maxu) {
				maxu = Tau * Nu;
				bestCity = cl[currentCity][i];
			}
		}
	}
	return bestCity;
}

int hasCandidatesLeft(int antIndex, int currentCity) {
	int i;
	int clvisited[CLSIZE];
	for (i = 0; i < CLSIZE; i++) {
		clvisited[i] = ant[antIndex].visited[cl[currentCity][i]];
	}
	for (i = 0; i < CLSIZE; i++) {
		if (!clvisited[i]) {
			return 1;
		}
	}
	return 0;
}

float dist(int x1, int y1, int x2, int y2) {
        float vdist, hdist;
        if (x1 > x2) vdist = (float)x1 - (float)x2;
        else vdist = (float)x2 - (float)x1;
        if (y1 > y2) hdist = (float)y1 - (float)y2;
        else hdist = (float)y2 - (float)y1;
        return sqrtf((float)(pow(hdist,2) + pow(vdist,2)));
}

double tourCost(int* tour) {
	int i;
	double cost = (double)0;
	for (i = 1; i < ins.dimension; i++) {
		cost += distanceMatrix[tour[i-1]][tour[i]];
	}
	return cost + distanceMatrix[tour[ins.dimension-1]][tour[0]];
}