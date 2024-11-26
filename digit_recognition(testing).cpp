// digit_recognition(testing).cpp : Defines the entry point for the console application.

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

vector<vector<long double>> universe;
vector<int> observation;          
vector<vector<long double>> clusters;

int CBsize = 32;            
// Tokhura weights for calculating Tokhura distance
double tokhuraWt[] = {1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};

const int ITERATIONS = 1000;
const int N = 5;   // Number of states
const int M = 32;  // Number of observations
const int T = 150;       // Length of observation sequence
long double threshold = 1e-30;

// Model parameters
long double A[N][N] = {0};    // State transition matrix
long double B[N][M] = {0};    // Observation probability matri
long double pi[N] = {1.0, 0.0, 0.0, 0.0, 0.0};      // Initial state distribution
char nums[10];

const int STRIDE = 320; 

const int NUM_GROUPS = 10;              // Ending index
    
int observation_sequence[1][T];
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
		observation.push_back(nearestClusterIndex);
    }
}

long double forward(vector<int>& observations, const long double A[N][N], const long double B[N][M], const long double pi[N]) {

    long double alphai[T][N] = {0.0};

    // Initialize alpha values
    for (int i = 0; i < N; ++i) {
        alphai[0][i] = pi[i] * B[i][observations[0]];
    }

    // Compute alpha values for each time step
    for (int t = 1; t < T; ++t) {
        for (int j = 0; j < N; ++j) {
            long double sum = 0.0;
            for (int i = 0; i < N; ++i) {
                sum += alphai[t - 1][i] * A[i][j];
            }
            alphai[t][j] = sum * B[j][observations[t]];
        }
    }

    // Calculate the final probability
    long double probability = 0.0;
    for (int i = 0; i < N; ++i) {
        probability += alphai[T - 1][i];
    }

    return probability;
}

void read_matrices_from_file(const char* filename) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        fprintf(stderr, "Error opening matrix file.\n");
        return;
    }

    char line[256];
    int i = 0, j = 0;

    // Read lines until finding "A matrix of model"
    while (fgets(line, sizeof(line), file)) {
        if (strstr(line, "A matrix of model:") != NULL) {
            break;
        }
    }

    // Read matrix A
    for (i = 0; i < N; ++i) {
        for (j = 0; j < N; ++j) {
            if (fscanf(file, "%lf", &A[i][j]) != 1) {
                fprintf(stderr, "Error reading A matrix data.\n");
                fclose(file);
                return;
            }
        }
    }

    // Skip lines until finding "B matrix of model"
    while (fgets(line, sizeof(line), file)) {
        if (strstr(line, "B matrix of model:") != NULL) {
            break;
        }
    }

    // Read matrix B
    for (i = 0; i < N; ++i) {
        for (j = 0; j < M; ++j) {
            if (fscanf(file, "%lf", &B[i][j]) != 1) {
                fprintf(stderr, "Error reading B matrix data.\n");
                fclose(file);
                return;
            }
        }
    }
	fclose(file);
}

vector<vector<long double>> readFileToVector(const char* filename) {
    vector<vector<long double>> data;
    ifstream inFile(filename);

    if (!inFile.is_open()) {
        cerr << "Error opening cluster file: " << filename << "\n";
        return data; // Return an empty vector on failure
    }

    string line;
    while (getline(inFile, line)) {
        istringstream iss(line);
        vector<long double> row;
        long double value;

        // Read values separated by spaces into the row vector
        while (iss >> value) {
            row.push_back(value);
        }

        // Add the row to the data vector

        if (!row.empty()) {
            data.push_back(row);
        }
    }

    inFile.close();
    return data;
}



short waveIn[16025 * 3];  // 3 seconds buffer
double wave[16025*3];

void PlayRecord()
{
    const int NUMPTS = 16025 * 3;   // 3 seconds
    int sampleRate = 16025;  
    HWAVEIN  hWaveIn;
    WAVEFORMATEX pFormat;
    pFormat.wFormatTag=WAVE_FORMAT_PCM;    
    pFormat.nChannels=1;    
	
    pFormat.nSamplesPerSec=sampleRate;      
    pFormat.nAvgBytesPerSec=sampleRate*2;  
    pFormat.nBlockAlign=2;                 
    pFormat.wBitsPerSample=16;             
    pFormat.cbSize=0;
	
    waveInOpen(&hWaveIn, WAVE_MAPPER, &pFormat, 0L, 0L, WAVE_FORMAT_DIRECT);
    WAVEHDR WaveInHdr;
    WaveInHdr.lpData = (LPSTR)waveIn;
    WaveInHdr.dwBufferLength = NUMPTS*2;
    WaveInHdr.dwBytesRecorded=0;
    WaveInHdr.dwUser = 0L;
    WaveInHdr.dwFlags = 0L;
    WaveInHdr.dwLoops = 0L;
	
    waveInPrepareHeader(hWaveIn, &WaveInHdr, sizeof(WAVEHDR));
    HWAVEOUT hWaveOut;

    cout << "Playing..." << endl;
    waveOutOpen(&hWaveOut, WAVE_MAPPER, &pFormat, 0, 0, WAVE_FORMAT_DIRECT);
    waveOutWrite(hWaveOut, &WaveInHdr, sizeof(WaveInHdr)); // Playing the data
    Sleep(3 * 1000); // Sleep for as long as there was recorded
    waveInClose(hWaveIn);
    waveOutClose(hWaveOut);
}

void StartRecord()
{
    const int NUMPTS = 16025 * 3;   // 3 seconds
    int sampleRate = 16025;  
    HWAVEIN hWaveIn;
    MMRESULT result;
    WAVEFORMATEX pFormat;
    pFormat.wFormatTag=WAVE_FORMAT_PCM;    
    pFormat.nChannels=1;                   
    pFormat.nSamplesPerSec=sampleRate;      
    pFormat.nAvgBytesPerSec=sampleRate*2;  
    pFormat.nBlockAlign=2;                 
    pFormat.wBitsPerSample=16;             
    pFormat.cbSize = 0;

    result = waveInOpen(&hWaveIn, WAVE_MAPPER, &pFormat, 0L, 0L, WAVE_FORMAT_DIRECT);
    WAVEHDR WaveInHdr;
    WaveInHdr.lpData = (LPSTR)waveIn;
    WaveInHdr.dwBufferLength = NUMPTS*2;
    WaveInHdr.dwBytesRecorded=0;
    WaveInHdr.dwUser = 0L;
    WaveInHdr.dwFlags = 0L;
    WaveInHdr.dwLoops = 0L;
    waveInPrepareHeader(hWaveIn, &WaveInHdr, sizeof(WAVEHDR));
    result = waveInAddBuffer(hWaveIn, &WaveInHdr, sizeof(WAVEHDR));
    result = waveInStart(hWaveIn);
    cout << "\n\n\nRecording for 3 seconds..." << endl;
    Sleep(1 * 1000);
	cout<<"Speak Now\n";
	Sleep(2 * 1000);
	cout<< "Recording Ended \n";
    waveInClose(hWaveIn);
    //PlayRecord();
}


int live_testing()
{
		universe.clear();       // Clear universe at the start of each test run
		observation.clear(); 
		StartRecord();
		int frameIndex;
		ll count = 16025 * 3;
		
		for (int i = 0; i < 16025*3; ++i) 
			wave[i] = static_cast<double>(waveIn[i]);
	
		DCshift(wave, count);
        normalization(wave, count);

		vector<long double> vt;
		frameIndex = 0;
        for (int offset = 0; offset < 150; offset++) 
		{
            int currentFrameIndex = frameIndex + offset * STRIDE;
            if (currentFrameIndex < 0 || currentFrameIndex + FRAME_SIZE > count) {
                printf("Frame %d (starting at index %d) is out of bounds, skipping calculations.\n",offset , currentFrameIndex);
                continue;
            }

            double frameData[FRAME_SIZE];
            for (int j = 0; j < FRAME_SIZE; j++) {
                frameData[j] = wave[currentFrameIndex + j];
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
		
	assignClusters(universe, clusters);
	
	int index = 0;
	for(int j = 0; j < T; j++)
		observation_sequence[0][j] = observation[index++];

	vector<int> observations;
	long double max_prob = -1.0; 
	int predicted_digit = -1; 

	observations.clear();
	for(int i = 0; i < T; i++)
		observations.push_back(observation_sequence[0][i]);
	for(int digit = 0; digit < NUM_GROUPS; digit++)
	{
		char filename[100];
		sprintf(filename, "Models/model_digit_%d", digit);
		read_matrices_from_file(filename);
			
		long double prob = forward(observations, A, B, pi);	
		//cout << "Digit " << digit  << ", probability is: " << prob << endl;
		if (prob > max_prob) 
		{
			max_prob = prob;
			predicted_digit = digit;
		}
    }
   
	return predicted_digit;
}




int main() {
	
	printf("Namaste!!!\n");
	printf("Kripya apna Phone Number ki ek ek krke Bole\n");
	printf("\n\n\n");
	Sleep(2 * 1000);
	clusters = readFileToVector("Models/cluster.txt");

	for(int i = 0; i < 10; i++)
	{	
		int x = live_testing();
		nums[i] = static_cast<char>(x + '0');
		cout << "\nApne " << nums[i] <<" bola hai" << endl ;
	}

	printf("Apka Phone Number Hai: ");
	for(int i = 0; i < 10; i++)
		printf("%c",nums[i]);
	cout<<endl;

	return 0;
}
