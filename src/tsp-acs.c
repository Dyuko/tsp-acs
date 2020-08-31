#include "tsp-acs.h"

int M;	//Número de hormigas

double** distanceMatrix;	//Almacena la distancia entre todos los nodos
double** tau;				//tau[i][j] representa el nivel de feromonas en el vértice(i,j)
int** cl;					//Matriz con todos los nodos (filas) y sus CLSIZE nodos más cercanos (columnas) ordenados de forma descendiente
InstanceData ins;			//Instacia TSP
Ant* ant;					//Matriz con el tour y los nodos visitados de cada hormiga

int runParallel;

int main(int argc, char *argv[]) {
	// Timer start
	struct timespec ts_start;
	clock_gettime(CLOCK_MONOTONIC, &ts_start);

	float bestKnownCost = FLT_MAX;	//Almacena el mejor costo del tour entre las hormigas
	int* bestKnownTour;				//Almacena el mejor tour entre las hormigas

	FILE *file;
	int i, j, k;
	
	//argv[1] indica el directorio al archivo que guarda los datos de una instancia de TSP
	//El archivo posee una extensión .tsp, cuyo formato es TSPLIB
	file = fopen(argv[1], "r"); //Se abre el archivo con la instancia de TSP							
	
	M = atoi(argv[2]);			//argv[2] indica el número de hormigas a enrutar
	
	runParallel = atoi(argv[3]);//argv[3] es un booleano que indica la paralización con hilos, 0 o 1
	
	//Se extraen los datos del archivo con la instancia TSP
	fscanf(file, "NAME : %[^\n]s", ins.name);				//ins.name almacena el nombre de la instancia TSP
	fscanf(file, "\nCOMMENT : %[^\n]s", ins.comment);		//ins.comment almacena un comentario adicional
	fscanf(file, "\nTYPE : %[^\n]s", ins.type);				//ins.type almacena el tipo de la instancia, en este caso "TSP"
	fscanf(file, "\nDIMENSION : %d", &ins.dimension);		//ins.dimension almacena el número de nodos
	fscanf(file, "\nEDGE_WEIGHT_TYPE : %[^\n]s", ins.wtype);//ins.wtype almacena cómo los pesos de borde están dadas, en este caso "EUC_2D", distancia euclidiana en 2-D
	fscanf(file, "\nNODE_COORD_SECTION");					//especifica las coordenadas asociada a cada nodo
	
	if (strcmp(ins.wtype,"EUC_2D")) {	//Comprueba que el tipo de los pesos de borde sea "EUC_2D"
		return 1;						//Si no es EUC_2D se cancela la ejecución retornando 1
	}

	if (runParallel) {					//Si se ejecuta en paralelo 
		printf("Running parallel version\n");
	}
	else {								//Si se ejecuta en secuencial 
		printf("Running standard version\n");
	}

	int** coord;	//coord almacenará las coordenadas de los nodos
	coord = (int**) malloc(ins.dimension * sizeof(int*));	//Reserva memoria para las filas, una fila por nodo
	bestKnownTour = (int*) malloc(ins.dimension * sizeof(int));	//Almacena el mejor tour entre las hormigas
	
	//Por cada una de las filas reserva memoria para dos columnas, coordenadas x e y del nodo
	for(i = 0; i <ins.dimension; i++) {
		coord[i] = (int*) malloc(2 * sizeof(int));
	}

	//Se lee las coordenadas desde el archivo y se almacena en la matriz coord 
	for (i = 0; i < ins.dimension; i++) {
		fscanf(file, "\n %*[^ ] %d %d", &coord[i][0], &coord[i][1]);
	}
	fclose(file);	//Cierra el archivo
	
	//distanceMatrix es una matriz con la distancia entre todos los nodos
	distanceMatrix = (double**) malloc(ins.dimension * sizeof(double*));		//Reserva memoria para la cantidad de nodos en fila
	for (i = 0; i < ins.dimension; i++) {
		distanceMatrix[i] = (double*) malloc(ins.dimension * sizeof(double));	//Por cada fila reserva memoria para la cantidad de nodos en columna
	}
	
	//Calcula la distancia entre todos los nodos y lo almacena en distanceMatrix
	for (i = 0; i < ins.dimension; i++){
		for ( j = i + 1; j < ins.dimension ; j++) {
			distanceMatrix[i][j] = dist(coord[i][0], coord[i][1], coord[j][0], coord[j][1]);
			distanceMatrix[j][i] = distanceMatrix[i][j];
		}
	}
	free(coord);	//Se libera la matriz con las coordenadas de los nodos

	//tau es una matriz con el nivel de feromonas en los vértices entre todos los nodos
	tau = (double**) malloc(ins.dimension*sizeof(double*));			//Reserva memoria para la cantidad de nodos en fila
	for (i = 0; i < ins.dimension; i++) {
		tau[i] = (double*) malloc(ins.dimension * sizeof(double));	//Por cada fila reserva memoria para la cantidad de nodos en columna
	}

	//Inicializa todos los vértices con el nivel de feromona inicial TAUO
	for (i = 0; i < ins.dimension; i++) {
		for( j = 0; j < ins.dimension; j++) {
			tau[i][j] = TAU0;
		}
	}
	
	srand(time(NULL));	//Establece la hora como la semilla del rand
	
	ant = (Ant*) malloc(M * sizeof(Ant));	//Reserva las filas para almacenar el tour y nodos visitados de cada hormiga
	
	//Inicializa las hormigas 
	for (i = 0; i < M; i++) {	//Por cada hormiga
		ant[i].tour = (int*) malloc(ins.dimension * sizeof(int));	//Reserva memoria para el tour de la hormiga
		ant[i].tour[0] = rand() % ins.dimension;					//random para establecer la posición inicial de la hormiga
		ant[i].visited = (int*) malloc(ins.dimension * sizeof(int));//Reserva memoria para los nodos visitados de la hormiga
		for (j = 1; j < ins.dimension; j++) {						//Indica que al comienzo la hormiga no visitó ningún nodo
			ant[i].visited[j] = 0;
		}
		ant[i].visited[ant[i].tour[0]] = 1;							//Indica como visitado la posición inicial de la hormiga
	}
	
	cl = (int**) malloc(ins.dimension * sizeof(int*));				//Reserva filas para la cantidad de nodos
	for (i = 0; i < ins.dimension; i++) {
		cl[i] = (int*) malloc(CLSIZE * sizeof(int));				//Por cada fila de cl, reserva CLSIZE columnas
	}
	int* visited = (int*) malloc(ins.dimension * sizeof(int));		//Vector auxiliar de visitados
	double nearestDistance;		//Menor distancia a un nodo
	int nearestNeighbor;		//Nodo con menor distancia

	/*
	* Cada fila de la matriz cl representa un nodo, cada fila posee CLSIZE columnas.
	* Las columnas representan los CLSIZE nodos vecinos más cercanos ordenados de forma ascendente (de izquierda a dercha)
	* En el siguiente for se carga los valores de la matriz cl
	*/
	for (i = 0; i < ins.dimension; i++) {	
		for(j = 0; j < ins.dimension; j++) {	
			visited[j] = 0;
		}
		for(j = 0; j < CLSIZE; j++) {	
			nearestDistance = DBL_MAX;	
			nearestNeighbor = 0;	
			for(k = 0; k < ins.dimension; k++) {	
				//Si el nodo no fue visitado y la distancia a ese nodo es aún menor y no es el nodo mismo
				if(!visited[k] && distanceMatrix[i][k] < nearestDistance && i != k) {
					nearestDistance = distanceMatrix[i][k];	//Actualiza la menor distancia
					nearestNeighbor = k;					//Indica cuál es el nodo con menor distancia
				}
			}
			visited[nearestNeighbor] = 1;	
			cl[i][j] = nearestNeighbor;	
		}
	}
	
	for (i = 0; i < ITERATIONS; i++) {	//Realiza las iteraciones
		if (runParallel) {	//Modo Paralelo
			//HotSpot
			#pragma omp parallel
			{
				#pragma omp single nowait
				{
					if(!i) {	//Para la primera iteración y el primer hilo imprime la cantidad de hilos
						int num_threads = omp_get_num_threads();
						printf("Number of threads: %d\n", num_threads);
					}
				}
				#pragma omp for
				for (j = 0; j < M; j++) {	//Divide la cantidad total de hormigas entre la cantidad de núcleos 
					for (k = 1; k < ins.dimension; k++) {					//Por cada movimiento a un nodo
						if (hasCandidatesLeft(j,ant[j].tour[k-1])) {		//Si hay un nodo entre los CLSIZE más cercanos no visitado
							moveAntTo(j, k, argMax(j, ant[j].tour[k-1]));	//Calcula y mueve la hormiga según la fórmula de probabilidad
						}
						else {												//Si ha visitado todos los CLSIZE nodos más cercanos
							moveAntTo(j, k, NN(j,k-1));						//Mueve la hormiga simplemente al nodo más cercano restante
						}
					}
				}
			}
		}
		else {	//Modo secuencial
			for (j = 0; j < M; j++) {									//Por cada hormiga
				for (k = 1; k < ins.dimension; k++) {					//Por cada movimiento a un nodo
					if (hasCandidatesLeft(j,ant[j].tour[k-1])) {		//Si hay un nodo entre los CLSIZE más cercanos no visitado
						moveAntTo(j, k, argMax(j, ant[j].tour[k-1]));	//Calcula y mueve la hormiga según la fórmula de probabilidad
					}
					else {												//Si ha visitado todos los CLSIZE nodos más cercanos
						moveAntTo(j, k, NN(j,k-1));						//Mueve la hormiga simplemente al nodo más cercano restante
					}
				}
			}
		}
		
		int bestAnt = 0;
		double lowerCost = DBL_MAX;
		//Encuentra el tour con menor costo
		for (j = 0; j < M; j++) {
			if (tourCost(ant[j].tour) < lowerCost) {
				bestAnt = j;
				lowerCost = tourCost(ant[j].tour);
			}
		}
		
		//Actualiza los niveles de feromonas
		updatePheromoneLevel(ant[bestAnt].tour);
		
		//Si en esta iteracion se encuentra un tour con menor costo se guarda
		if (lowerCost < bestKnownCost) {
			bestKnownCost = lowerCost;
			for (k = 0; k < ins.dimension; k++) {
				bestKnownTour[k] = ant[bestAnt].tour[k];
			}
		}

		//Reinicia las hormigas
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
	printf("%d ", bestKnownTour[0]);
	printf("\nLength = %f\n", bestKnownCost);

	free(ant);

	// timer stop
	struct timespec ts_stop;
	clock_gettime(CLOCK_MONOTONIC, &ts_stop);
	double start = (double)ts_start.tv_sec + (double)ts_start.tv_nsec/1000000000.0;
	double stop = (double)ts_stop.tv_sec + (double)ts_stop.tv_nsec/1000000000.0;
	double elapsed = (stop - start);

	// display time
	printf ("Time = %f\n", elapsed);
	return 0;
}

/*
* Actualiza los niveles de feromonas
*/
void updatePheromoneLevel(int* tour) {
	int i;
	for (i = 1; i < ins.dimension; i++) {
		tau[tour[i-1]][tour[i]] = (1-RHO) * tau[tour[i-1]][tour[i]] + RHO * (Q/tourCost(tour));
		tau[tour[i]][tour[i-1]] = (1-RHO) * tau[tour[i]][tour[i-1]] + RHO * (Q/tourCost(tour));
	}
}

/*
* Dada una hormiga y su nodo actual, calcula y retorna el nodo no visitado más cercano
*/
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

/*
* Dada una hormiga, el slot actual en su tour y el siguiente nodo a visitar, mueve la hormiga al destino y actualiza los niveles de feromonas
*/
void moveAntTo(int antIndex, int slot, int city){
	ant[antIndex].tour[slot] = city;
	ant[antIndex].visited[city] = 1;
	tau[ant[antIndex].tour[slot]][ant[antIndex].tour[slot-1]] = (1-RHO) * tau[ant[antIndex].tour[slot]][ant[antIndex].tour[slot-1]] + RHO * TAU0;
	tau[ant[antIndex].tour[slot-1]][ant[antIndex].tour[slot]] = (1-RHO) * tau[ant[antIndex].tour[slot-1]][ant[antIndex].tour[slot]] + RHO * TAU0;
}

/*
* Dada una hormiga y su nodo actual, calcula y retorna el nodo a donde se moverá la hormiga según la fórmula de probabilidad
*/
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

/*
* Dada una hormiga y su nodo actual, verifica si hay algún nodo cercano que todavía no ha visitado
* En caso de encontrar un nodo cercano no visitado retorna 1, en caso contrario 0
*/
int hasCandidatesLeft(int antIndex, int currentCity) {
	int i;
	int clvisited[CLSIZE];
	for (i = 0; i < CLSIZE; i++) {	//Verifica si la hormiga ya visitó los nodos más cercanos de su nodo actual
		clvisited[i] = ant[antIndex].visited[cl[currentCity][i]];
	}
	for (i = 0; i < CLSIZE; i++) {	//Si hay un nodo cercano que no ha visitado entonces retorna 1
		if (!clvisited[i]) {
			return 1;
		}
	}
	return 0;
}

/*
* Calcula la distancia total entre dos coordenadas 
*/
float dist(int x1, int y1, int x2, int y2) {
        float vdist, hdist;
        if (x1 > x2) vdist = (float)x1 - (float)x2;
        else vdist = (float)x2 - (float)x1;
        if (y1 > y2) hdist = (float)y1 - (float)y2;
        else hdist = (float)y2 - (float)y1;
        return sqrtf((float)(pow(hdist,2) + pow(vdist,2)));
}

/*
* Calcula es costo total asociado al tour de una hormiga
*/
double tourCost(int* tour) {
	int i;
	double cost = (double)0;
	for (i = 1; i < ins.dimension; i++) {
		cost += distanceMatrix[tour[i-1]][tour[i]];
	}
	return cost + distanceMatrix[tour[ins.dimension-1]][tour[0]];
}