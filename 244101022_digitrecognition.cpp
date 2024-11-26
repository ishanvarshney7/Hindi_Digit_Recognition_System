#include "StdAfx.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <string.h>
#include <direct.h>  
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <Windows.h>

using namespace std;
#pragma comment(lib, "winmm.lib")
typedef long long ll;
#define P 12
#define Max_CBsize 32
#define DELTA 0.0001      // Convergence threshold for distortion
#define EPS 0.03          // Epsilon value used for codebook splitting
#define FRAME_SIZE 320
#define ITERATIONS 1000

const int NUM_GROUPS = 10;    // Total number of groups
const int NUM_FILES_PER_GROUP = 100; // Files per group
const int MAX_FILENAME_LENGTH = 100; // Max length for each filename

vector<vector<long double>> universe;
vector<vector<int>> assigned_idx;
vector<int> observation;          
vector<vector<long double>> clusters;

int CBsize = 1;            
const int STRIDE = 320;
// Tokhura weights for calculating Tokhura distance
double tokhuraWt[] = {1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};

const int N = 5;   // Number of states
const int M = 32;  // Number of observations
const int T = 150;       // Length of observation sequence
long double threshold = 1e-30;
long double P_prev = -1.0, P_star = 0;


// Model parameters
long double A[N][N];    // State transition matrix
long double B[N][M];    // Observation probability matri
long double pi[N];      // Initial state distribution

// Global arrays for Forward, Backward, Viterbi, and Baum-Welch algorithms
long double** alpha;    // Forward probabilities
long double** beta;     // Backward probabilities
long double** gamma;    // Baum-Welch gamma values
long double** delta;    // Viterbi delta values
int** psi;         // Viterbi psi values
long double*** xi;      // Baum-Welch xi values

int observation_sequence[NUM_FILES_PER_GROUP * NUM_GROUPS][T];

// Function to emove the DC offset (mean) from the data array
void DCshift(double* dataArray, ll elementCount) {
    double sum = 0;
    for (ll i = 0; i < elementCount; i++) {
        sum += dataArray[i];
    }
    double mean = sum / elementCount;

    //printf("DC shift: %lf\n", mean);
    for (ll i = 0; i < elementCount; i++) {
        dataArray[i] -= mean;
    }
}

// Function to normalize the data to the range [-5000, 5000]
void normalization(double* dataArray, ll elementCount) {
    double max = fabs(dataArray[0]);
    for (ll i = 1; i < elementCount; i++) {
        double absValue = fabs(dataArray[i]);
        if (absValue > max)
            max = absValue;
    }
    //printf("Max value before normalization: %lf\n", max);

    double normalizationFactor = 5000.0 / max;
    //printf("Normalization factor: %lf\n", normalizationFactor);

    for (ll i = 0; i < elementCount; i++) {
        dataArray[i] *= normalizationFactor;
    }

    double newMax = fabs(dataArray[0]);
    for (ll i = 1; i < elementCount; i++) {
        double absValue = fabs(dataArray[i]);
        if (absValue > newMax)
            newMax = absValue;
    }
    //printf("Max after normalization: %lf\n\n", newMax);
}

void calculateRi(double* Arr, int samples, double* R) {
    for (int k = 0; k < 13; k++) {
        R[k] = 0.0;
        for (int i = 0; i < samples - k; i++) {
            R[k] += Arr[i] * Arr[i + k];
        }
    }
}

void calculateAi(double* R, double* A) {
    double E[13] = {0.0};
    double K[13] = {0.0};
    double alpha[13][13] = {0.0};

    E[0] = R[0];
    for (int i = 1; i <= 12; i++) {
        double val = 0.0;
        for (int j = 1; j < i; j++)
            val += alpha[j][i - 1] * R[i - j];
        K[i] = (R[i] - val) / E[i - 1];
        alpha[i][i] = K[i];
        if (i > 1) {
            for (int j = 1; j < i; j++) {
                alpha[j][i] = alpha[j][i - 1] - K[i] * alpha[i - j][i - 1];
            }
        }
        E[i] = (1 - K[i] * K[i]) * E[i - 1];
    }
    for (int i = 1; i <= 12; i++) {
        A[i - 1] = alpha[i][12];
    }
}

void concatFileName(char* dest, const char* baseName, const char* suffix) {
    while (*baseName) {
        *dest++ = *baseName++;
    }
    while (*suffix) {
        *dest++ = *suffix++;
    }
    *dest = '\0';
}

void calculateCi(double* A, double* C, double* R) {
    C[0] = log(R[0]);  // Cepstral coefficient C0 is the log of the energy of the signal

    // Recursive calculation of cepstral coefficients
    for (int m = 1; m <= 12; m++) {
        C[m] = A[m - 1];
        for (int k = 1; k < m; k++) {
            C[m] += ((double)k / m) * C[k] * A[m - k - 1];
        }
    }
}

void applyRaisedSineWindow(double* C, int Ma) {
    for (int m = 1; m <= Ma; m++) {
        double windowFactor = 1.0 + (Ma / 2.0) * sin(3.14 * m / Ma);
        C[m] *= windowFactor;
    }
}

long double calcTokhuraDistance(vector<long double>& currentRegion, vector<long double>& currentUniverse) {
    long double distance = 0;
    for (int i = 0; i < P; i++) {
        long double diff = currentRegion[i] - currentUniverse[i];
        distance += tokhuraWt[i] * diff * diff;
    }
    return distance;
}

// Assigns universe vectors to the closest cluster using Tokhura distance
void assignClusters(vector<vector<long double>>& universe, vector<vector<long double>>& clusters) {
    int N = universe.size(); // Number of vectors in the universe
	observation.clear();
    for (int i = 0; i < N; i++) {
        int nearestClusterIndex = 0;
        long double minDistance = LDBL_MAX;
        for (int j = 0; j < CBsize; j++) {
            long double distance = calcTokhuraDistance(clusters[j], universe[i]);
            if (distance < minDistance) {
                minDistance = distance;
                nearestClusterIndex = j;
            }
        }
        assigned_idx[nearestClusterIndex].push_back(i); // Assign vector to the nearest cluster	
		observation.push_back(nearestClusterIndex);
    }
}


// Calculates the average total distortion for the current clustering
long double totalDistortion(vector<vector<long double>>& universe, vector<vector<long double>>& clusters) {
    long double total_dist = 0;
    long double totalVectors = 0;

    for (int i = 0; i < CBsize; i++) {
        int clusterSize = assigned_idx[i].size();
        totalVectors += clusterSize;

        for (int j = 0; j < clusterSize; j++) {
            int vectorIndex = assigned_idx[i][j];
            total_dist += calcTokhuraDistance(clusters[i], universe[vectorIndex]);
        }
    }
    return total_dist / totalVectors; // Return the average distortion
}

// Updates the centroids (clusters) by averaging the assigned vectors
void updateClusters(vector<vector<long double>>& universe, vector<vector<long double>>& clusters) {
    for (int i = 0; i < CBsize; i++) {
        vector<long double> newCentroid(P, 0.0);
        int clusterSize = assigned_idx[i].size();

        for (int j = 0; j < clusterSize; j++) {
            int vectorIndex = assigned_idx[i][j];
            for (int k = 0; k < P; k++) {
                newCentroid[k] += universe[vectorIndex][k];
            }
        }

        // Compute the new centroid by averaging the assigned vectors
        for (int k = 0; k < P; k++) {
            clusters[i][k] = newCentroid[k] / (long double)clusterSize;
        }
    }
}

// Implements Lloyd's K-Means algorithm to iteratively minimize distortion
void lloydsKmeansAlgo(vector<vector<long double>>& universe, vector<vector<long double>>& clusters) {
    assigned_idx.clear();
    assigned_idx.resize(CBsize);

    assignClusters(universe, clusters);
    long double currentDist = totalDistortion(universe, clusters);
    int round = 1;
    long double oldDist = 0;

    printf("Cycle %d has total distortion: %Lf\n", round, currentDist);

    // Repeat until the distortion change is below the threshold (DELTA)
    while (abs(currentDist - oldDist) > DELTA) {
        updateClusters(universe, clusters);
        assigned_idx.clear();
        assigned_idx.resize(CBsize);

        assignClusters(universe, clusters);
        oldDist = currentDist;
        currentDist = totalDistortion(universe, clusters);
        round++;

        printf("Cycle %d has total distortion: %Lf, DELTA: %Lf\n", round, currentDist, abs(currentDist - oldDist));
    }
}

// Computes the centroid of the entire universe
vector<long double> findCentroidOfUniverse(vector<vector<long double>>& universe) {
    vector<long double> centroid(P, 0.0);
    int N = universe.size();

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < P; j++) {
            centroid[j] += universe[i][j];
        }
    }

    for (int j = 0; j < P; j++) {
        centroid[j] /= N; // Compute the average for each dimension
    }

    return centroid;
}

// Updates the codebook by splitting each centroid into two new centroids
void updateCodebook(vector<vector<long double>>& clusters, vector<vector<long double>>& copyCB) {
    for (int i = 0; i < copyCB.size(); i++) {
        vector<long double> positive(P), negative(P);

        for (int j = 0; j < P; j++) {
            positive[j] = (1 + EPS) * copyCB[i][j];
            negative[j] = (1 - EPS) * copyCB[i][j];
        }

        clusters.push_back(positive);
        clusters.push_back(negative);
    }
}

// Implements the Linde-Buzo-Gray (LBG) algorithm to iteratively increase the codebook size
void LBGAlgo(vector<vector<long double>>& universe, vector<vector<long double>>& clusters) {
    while (CBsize <= Max_CBsize) {
        printf("\n\n\nCodebook of size: %d\n", CBsize);
        lloydsKmeansAlgo(universe, clusters);

    //    printf("\nGenerated Codebook:\n");
        for (int i = 0; i < clusters.size(); i++) {
            for (int j = 0; j < clusters[i].size(); j++) {
                //printf("%Lf ", clusters[i][j]);
            }
		//	printf("\n");
        }
		
        if (CBsize >= Max_CBsize) {
            return;
        }

        vector<vector<long double>> copyCB(clusters);
        clusters.clear();
        CBsize *= 2;

        updateCodebook(clusters, copyCB);
    }
}

// Function to initialize dynamic global arrays
void initializeGlobalArrays(int T) {
    alpha = new long double*[T];
    beta = new long double*[T];
    gamma = new long double*[T];
    delta = new long double*[T];
    psi = new int*[T];
    xi = new long double**[T - 1];

    for (int t = 0; t < T; ++t) {
        alpha[t] = new long double[N];
        beta[t] = new long double[N];
        gamma[t] = new long double[N];
        delta[t] = new long double[N];
        psi[t] = new int[N];

        if (t < T - 1) {
            xi[t] = new long double*[N];
            for (int i = 0; i < N; ++i) {
                xi[t][i] = new long double[N];
            }
        }
    }
}

void deallocateGlobalArrays(int T) {
    for (int t = 0; t < T; ++t) {
        delete[] alpha[t];
        delete[] beta[t];
        delete[] gamma[t];
        delete[] delta[t];
        delete[] psi[t];

        if (t < T - 1) {
            for (int i = 0; i < N; ++i) {
                delete[] xi[t][i];
            }
            delete[] xi[t];
        }
    }
    delete[] alpha;
    delete[] beta;
    delete[] gamma;
    delete[] delta;
    delete[] psi;
    delete[] xi;
}

// Function to read model (A, B, and pi) from file
void readModelFromFile(const string& filename) {
    ifstream file(filename);
    string line;

    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    bool isAMatrix = false, isBMatrix = false, isPiMatrix = false;
    int rowA = 0, rowB = 0, rowPi = 0;

    while (getline(file, line)) {
        istringstream iss(line);
        if (line.find("A Matrix") != string::npos) {
            isAMatrix = true;
            isBMatrix = false;
            isPiMatrix = false;
            rowA = 0;
        } else if (line.find("B Matrix") != string::npos) {
            isAMatrix = false;
            isBMatrix = true;
            isPiMatrix = false;
            rowB = 0;
        } else if (line.find("PI Matrix") != string::npos) {
            isAMatrix = false;
            isBMatrix = false;
            isPiMatrix = true;
            rowPi = 0;
        } else if (isAMatrix) {
            for (int j = 0; j < N; ++j) {
                iss >> A[rowA][j];
            }
            rowA++;
        } else if (isBMatrix) {
            for (int j = 0; j < M; ++j) {
                iss >> B[rowB][j];
            }
            rowB++;
        } else if (isPiMatrix) {
            for (int i = 0; i < N; ++i) {
                iss >> pi[i];
            }
        }
    }

    file.close();
}

// Forward Algorithm (Problem 1)
long double forward(const vector<int>& observations) {
    int T = observations.size();

    for (int i = 0; i < N; ++i) {
        alpha[0][i] = pi[i] * B[i][observations[0]];
    }

    for (int t = 1; t < T; ++t) {
        for (int j = 0; j < N; ++j) {
            long double sum = 0.0;
            for (int i = 0; i < N; ++i) {
                sum += alpha[t - 1][i] * A[i][j];
            }
            alpha[t][j] = sum * B[j][observations[t]];
        }
    }

    long double probability = 0.0;
    for (int i = 0; i < N; ++i) {
        probability += alpha[T - 1][i];
    }

    return probability;
}

// Backward Algorithm (Problem 3)
long double backward(const vector<int>& observations) {
    int T = observations.size();

    for (int i = 0; i < N; ++i) {
        beta[T - 1][i] = 1.0;
    }

    for (int t = T - 2; t >= 0; --t) {
        for (int i = 0; i < N; ++i) {
            long double sum = 0.0;
            for (int j = 0; j < N; ++j) {
                sum += A[i][j] * B[j][observations[t + 1]] * beta[t + 1][j];
            }
            beta[t][i] = sum;
        }
    }

    long double probability = 0.0;
    for (int i = 0; i < N; ++i) {
        probability += pi[i] * B[i][observations[0]] * beta[0][i];
    }

    return probability;
}

// Viterbi Algorithm (Problem 2)
long double viterbi(const vector<int>& observations) {
    int T = observations.size();

    for (int i = 0; i < N; ++i) {
        delta[0][i] = pi[i] * B[i][observations[0]];
        psi[0][i] = 0;
    }

    for (int t = 1; t < T; ++t) {
        for (int j = 0; j < N; ++j) {
            long double maxVal = -1;
            int maxArg = 0;
            for (int i = 0; i < N; ++i) {
                long double val = delta[t - 1][i] * A[i][j];
                if (val > maxVal) {
                    maxVal = val;
                    maxArg = i;
                }
            }
            delta[t][j] = maxVal * B[j][observations[t]];
            psi[t][j] = maxArg;
        }
    }

    vector<int> stateSequence(T);
    long double maxProb = 0.0;
    int maxState = 0;
    for (int i = 0; i < N; ++i) {
        if (delta[T - 1][i] > maxProb) {
            maxProb = delta[T - 1][i];
            maxState = i;
        }
    }
    stateSequence[T - 1] = maxState;

    for (int t = T - 2; t >= 0; --t) {
        stateSequence[t] = psi[t + 1][stateSequence[t + 1]];
    }

    return maxProb;
}

// Baum-Welch Algorithm (Problem 3)
void baumWelch(const vector<int>& observations) {
    int T = observations.size();
	int count = 1;
    while(P_prev < P_star && count < ITERATIONS) {
		printf("iteration %d\n", count);
		count++;
        // Step 1: Forward
        forward(observations);

        // Step 2: Backward
        backward(observations);

        // Step 3: Compute gamma and xi
        for (int t = 0; t < T - 1; ++t) {
            long double denom = 0.0;
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) {
                    denom += alpha[t][i] * A[i][j] * B[j][observations[t + 1]] * beta[t + 1][j];
                }
            }
            for (int i = 0; i < N; ++i) {
                gamma[t][i] = 0.0;
                for (int j = 0; j < N; ++j) {
                    xi[t][i][j] = (alpha[t][i] * A[i][j] * B[j][observations[t + 1]] * beta[t + 1][j]) / denom;
                    gamma[t][i] += xi[t][i][j];
                }
            }
        }

        // Last gamma
        long double denom = 0.0;
        for (int i = 0; i < N; ++i) {
            denom += alpha[T - 1][i];
        }
        for (int i = 0; i < N; ++i) {
            gamma[T - 1][i] = alpha[T - 1][i] / denom;
        }

        // Step 4: Re-estimate pi, A, and B
        for (int i = 0; i < N; ++i) {
            pi[i] = gamma[0][i];
        }
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                long double numer = 0.0, denom = 0.0;
                for (int t = 0; t < T - 1; ++t) {
                    numer += xi[t][i][j];
                    denom += gamma[t][i];
                }
                A[i][j] = numer / denom;
            }
        }
        for (int i = 0; i < N; ++i) {
            for (int k = 0; k < M; ++k) {
                long double numer = 0.0, denom = 0.0;
                for (int t = 0; t < T; ++t) {
                    if (observations[t] == k) {
                        numer += gamma[t][i];
                    }
                    denom += gamma[t][i];
                }
                B[i][k] = numer / denom;
            }
		}
		P_prev = P_star;
		P_star = viterbi(observations);
        
    }
}

// Threshold adjustment for B matrix
void applyThresholdToBMatrix() {
    for (int i = 0; i < N; ++i) {
        long double rowSum = 0.0;
        int maxIndex = 0;
        long double maxValue = 0.0;

        // First, apply threshold and calculate the row sum
        for (int j = 0; j < M; ++j) {
            if (B[i][j] < threshold) {
                B[i][j] = threshold;  // Apply threshold
            }
            rowSum += B[i][j];

            // Track the maximum value in the row
            if (B[i][j] > maxValue) {
                maxValue = B[i][j];
                maxIndex = j;
            }
        }

        // If the row sum exceeds 1, adjust the maximum value
        if (rowSum > 1.0) {
            long double excess = rowSum - 1.0;
            B[i][maxIndex] -= excess;  // Subtract excess from the maximum value
        }
    }
}




void main()
{
	
    char fileNames[NUM_GROUPS * NUM_FILES_PER_GROUP][MAX_FILENAME_LENGTH];
    int fileIndex = 0;

    for (int group = 0; group < NUM_GROUPS; ++group) {
        for (int fileNum = 1; fileNum <= NUM_FILES_PER_GROUP; ++fileNum) {
            sprintf(fileNames[fileIndex], "Recordings/244101022_H_%d_%d", group, fileNum);
            fileIndex++;
        }
    }
    int numFiles = sizeof(fileNames) / sizeof(fileNames[0]);  
	
    double* data = NULL;
    ll capacity = 100;

    for (int fileIdx = 0; fileIdx < numFiles; fileIdx++) \
	{
        char inputFileName[100];
        char outputFileName[100];

        strcpy(inputFileName, fileNames[fileIdx]);
        strcat(inputFileName, ".txt");
       
        FILE* file = fopen(inputFileName, "r");
        if (file == NULL) {
            printf("Error: Could not open file %s\n", inputFileName);
            continue;
        }

        char buffer[256];
        for (int i = 0; i < 5; i++) {
            if (fgets(buffer, sizeof(buffer), file) == NULL) {
                printf("Error: Not enough lines in the file to skip.\n");
                fclose(file);
                continue;
            }
        }

        data = (double*)malloc(capacity * sizeof(double));
        if (data == NULL) {
            printf("Error: Memory allocation failed.\n");
            fclose(file);
            continue;
        }

        ll count = 0;
        while (fscanf(file, "%lf", &data[count]) == 1) {
            count++;
            if (count >= capacity) {
                capacity *= 2;
                double* newData = (double*)realloc(data, capacity * sizeof(double));
                if (newData == NULL) {
                    printf("Error: Memory reallocation failed.\n");
                    free(data);
                    fclose(file);
                    continue;
                }
                data = newData;
            }
        }
        fclose(file);

        int frameIndex;

       
		DCshift(data, count);
        normalization(data, count);
		vector<long double> vt;
		frameIndex = 0;
        for (int offset = 0; offset < T; offset++) 
		{
            int currentFrameIndex = frameIndex + offset * STRIDE;
            if (currentFrameIndex < 0 || currentFrameIndex + FRAME_SIZE > count) {
                printf("Frame %d (starting at index %d) is out of bounds, skipping calculations.\n",offset , currentFrameIndex);
                continue;
            }

            double frameData[FRAME_SIZE];
            for (int j = 0; j < FRAME_SIZE; j++) {
                frameData[j] = data[currentFrameIndex + j];
            }

            double R[13] = {0.0};
            double A[12] = {0.0};
			double C[13] = {0.0};
            calculateRi(frameData, FRAME_SIZE, R);
            calculateAi(R, A);
			calculateCi(A, C, R);
			applyRaisedSineWindow(C, 12);
            for (int i = 1; i <= 12; i++)
				vt.push_back(C[i]);
			universe.push_back(vt);
			vt.clear();
		}
        free(data);
    }
	assigned_idx.clear();
    assigned_idx.resize(CBsize);
    


    // Start with one cluster at the centroid of the universe
    clusters.push_back(findCentroidOfUniverse(universe));

    // Execute the LBG algorithm to generate the final codebook
    LBGAlgo(universe, clusters);
	printf("LBG HAS CONVERGED!!!\n");

	int index = 0;
	for(int i = 0; i < NUM_FILES_PER_GROUP * NUM_GROUPS; i++)
	{
		for(int j = 0; j < T; j++)
			observation_sequence[i][j] = observation[index++];
	}
		
	
	vector<int> observations;
	for(int digit = 0; digit < NUM_GROUPS; digit++)
	{
		long double B_final[N][M] = {0}, A_final[N][N] = {0};    
		for(int samples = 0; samples < NUM_FILES_PER_GROUP; samples++)
		{
			observations.clear();
			for(int i = 0; i < T; i++)
				observations.push_back(observation_sequence[NUM_FILES_PER_GROUP * digit + samples][i]);
	
			readModelFromFile("Initial_Model.txt");
			initializeGlobalArrays(T);
			
			P_prev = -1;
			P_star = 0;
			// Baum-Welch Algorithm (Problem 3)
			baumWelch(observations);
			cout<<endl;			
			for (int i = 0; i < N; ++i) {
				for (int j = 0; j < N; ++j) 
					A_final[i][j] +=  A[i][j];
			}
			applyThresholdToBMatrix();
			for (int i = 0; i < N; ++i) {
				for (int j = 0; j < M; ++j) 
					B_final[i][j] += B[i][j];
			}
			deallocateGlobalArrays(T);
		}
	
		char filename[50];
		FILE *outFile;
		int i, j;

	    /* Construct the filename */
		sprintf(filename, "Models/model_digit_%d", digit);

	    /* Open the file in write mode */
		outFile = fopen(filename, "w");
	    if (outFile == NULL) {
		    fprintf(stderr, "Error opening file: %s\n", filename);
			return ; 
		}

		/* Write the A matrix of model to the file */
		fprintf(outFile, "A matrix of model:\n");
		for (i = 0; i < N; ++i) {
			for (j = 0; j < N; ++j) {
				A_final[i][j] = A_final[i][j] / (double)NUM_FILES_PER_GROUP;
				fprintf(outFile, "%35.30lf ",A_final[i][j]);
			}
			fprintf(outFile, "\n");
		}

		/* Write the B matrix of model to the file */
		fprintf(outFile, "\nB matrix of model:\n");
		for (i = 0; i < N; ++i) {
			for (j = 0; j < M; ++j) {
				B_final[i][j] = B_final[i][j] / (double)NUM_FILES_PER_GROUP;
				fprintf(outFile,"%35.30lf ",B_final[i][j]);
			}
			fprintf(outFile, "\n");
		}

		printf("Model for digit %d saved\n", digit);

		/* Close the file */
		fclose(outFile);
	}
	// Open file for writing
    ofstream outFile("Models/cluster.txt");
    if (outFile.is_open()) {
        // Iterate through the 2D vector and write each value to the file
        for (size_t i = 0; i < clusters.size(); ++i) {
            for (size_t j = 0; j < clusters[i].size(); ++j) {
                outFile << clusters[i][j];
                if (j < clusters[0].size() - 1) {
                    outFile << " ";  // Add a space between elements in the same row
                }
            }
            outFile << "\n";  // New line after each row
        }
        outFile.close();
        cout << "Data successfully written to cluster.txt\n";
    } else {
        cerr << "Error opening file for writing.\n";
    }
}