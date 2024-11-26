// Record.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <Windows.h>
#include <mmsystem.h>
#include <direct.h>  // For creating directories
#pragma comment(lib, "winmm.lib")
using namespace std;

#define frame_size 320
#define LENGTH_WAV 16025 * 3
#define MAX_OCCURRENCE 3  // Occurrence limit

short int waveIn[LENGTH_WAV];  // Array where the sound sample will be stored

void PlayRecord() {
    const int NUMPTS = LENGTH_WAV;
    int sampleRate = 16025;

    HWAVEOUT hWaveOut;
    WAVEFORMATEX pFormat;
    pFormat.wFormatTag = WAVE_FORMAT_PCM;
    pFormat.nChannels = 1;
    pFormat.nSamplesPerSec = sampleRate;
    pFormat.nAvgBytesPerSec = sampleRate * 2;
    pFormat.nBlockAlign = 2;
    pFormat.wBitsPerSample = 16;
    pFormat.cbSize = 0;

    if (waveOutOpen(&hWaveOut, WAVE_MAPPER, &pFormat, 0L, 0L, WAVE_FORMAT_DIRECT) != MMSYSERR_NOERROR) {
        printf("Failed to open waveform output device.\n");
        return;
    }

    WAVEHDR WaveOutHdr;
    WaveOutHdr.lpData = (LPSTR)waveIn;
    WaveOutHdr.dwBufferLength = NUMPTS * 2;
    WaveOutHdr.dwBytesRecorded = 0;
    WaveOutHdr.dwUser = 0L;
    WaveOutHdr.dwFlags = 0L;
    WaveOutHdr.dwLoops = 0L;
    waveOutPrepareHeader(hWaveOut, &WaveOutHdr, sizeof(WAVEHDR));

    printf("Playing...\n");
    waveOutWrite(hWaveOut, &WaveOutHdr, sizeof(WaveOutHdr));

    Sleep(3 * 1000);  // Sleep for duration of playback

    waveOutClose(hWaveOut);
}

void StartRecord() {
    const int NUMPTS = LENGTH_WAV;
    int sampleRate = 16025;

    HWAVEIN hWaveIn;
    MMRESULT result;

    WAVEFORMATEX pFormat;
    pFormat.wFormatTag = WAVE_FORMAT_PCM;
    pFormat.nChannels = 1;
    pFormat.nSamplesPerSec = sampleRate;
    pFormat.nAvgBytesPerSec = sampleRate * 2;
    pFormat.nBlockAlign = 2;
    pFormat.wBitsPerSample = 16;
    pFormat.cbSize = 0;

    result = waveInOpen(&hWaveIn, WAVE_MAPPER, &pFormat, 0L, 0L, WAVE_FORMAT_DIRECT);

    if (result != MMSYSERR_NOERROR) {
        printf("Failed to open waveform input device.\n");
        return;
    }

    WAVEHDR WaveInHdr;
    WaveInHdr.lpData = (LPSTR)waveIn;
    WaveInHdr.dwBufferLength = NUMPTS * 2;
    WaveInHdr.dwBytesRecorded = 0;
    WaveInHdr.dwUser = 0L;
    WaveInHdr.dwFlags = 0L;
    WaveInHdr.dwLoops = 0L;
    waveInPrepareHeader(hWaveIn, &WaveInHdr, sizeof(WAVEHDR));

    result = waveInAddBuffer(hWaveIn, &WaveInHdr, sizeof(WAVEHDR));
    if (result != MMSYSERR_NOERROR) {
        printf("Failed to add buffer to waveform input device.\n");
        waveInClose(hWaveIn);
        return;
    }

    result = waveInStart(hWaveIn);
    if (result != MMSYSERR_NOERROR) {
        printf("Failed to start waveform input device.\n");
        waveInClose(hWaveIn);
        return;
    }

    printf("\n\n\nRecording for 3 seconds...\n");
    Sleep(1 * 1000);  // Wait until finished recording
	cout<<"Speak Now\n";
	Sleep(2 * 1000);
    waveInClose(hWaveIn);
	PlayRecord();
}


void save_data(int digit, int occurrence) {
    char filepath[1000];
    sprintf(filepath, "Hindi/txt/244101022_H_%d_%d.txt", digit, occurrence);  // File path for text file

    FILE* file = fopen(filepath, "w+");
    if (!file) {
        printf("Error opening file\n");
        exit(1);
    }

    for (int i = 0; i < LENGTH_WAV; i++) {
        fprintf(file, "%d\n", waveIn[i]);
    }

    fclose(file);
    printf("Saved as %s\n", filepath);
}


int main() {
    int occurrences[10] = {0};  // Array to store occurrences for digits 0-9
    int digit;

    // Create directory structure: 'Hindi', 'Hindi/wav', 'Hindi/txt'
    _mkdir("Hindi");
    _mkdir("Hindi/txt");

    printf("Enter a digit (0-9): ");
    scanf("%d", &digit);

    if (digit < 0 || digit > 9) {
        printf("Invalid input! Please enter a digit between 0 and 9.\n");
        return 1;
    }

    for (int i = 1; i <= 100; i++) {
        StartRecord();
        save_data(digit, i);
        printf("Recorded and saved occurrence %d for digit %d\n", i, digit);
    }

    printf("All occurrences for digit %d have been recorded and saved.\n", digit);
	return 0;
}