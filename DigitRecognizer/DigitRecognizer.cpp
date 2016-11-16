// DigitRecognizer.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <typeinfo>
#include <vector>
#include <math.h>
#include <iomanip>
using namespace std;

int NORM_CONSTANT = 10000;
double THRESHOLD_CONSTANT = 10;
int FRAME_SIZE = 320;

vector<float> v, temp, noise;
ofstream logfile;

// Given a frame of 320 samples, return ci (find Ri,ai,ci)
vector<double> find_cepstral(vector<double> frame){
	vector<double> R, c;
	vector<vector<double>> a(100, vector<double> (100, 0));
	vector<double> E(100, 0);
	for (int i = 0; i <= 12; i++){
		R.push_back(1);
	}
	for (int i = 0; i <= 18; i++){
		c.push_back(0);
	}

	// Find auto-correlation co-efficients
	for (int m = 0; m <= 12; m++){
		R[m] = 0;
		for (int n = 0; n < (FRAME_SIZE - 1 - m); n++){
			R[m] += frame[n] * frame[n + m];
		}
	}


	double K[100];
	int p = 12;
	E[0] = R[0];
	for (int i = 1; i < (p + 1); i++){
		//DURBIN(i, R, a, E);
		//the algorithm starts here
		E[0] = R[0];
		for (int j = 1; j <= i; j++){
			if (j == 1){
				//cout << "Rn[1] "<<Rn[1]<<endl;
				K[1] = R[1] / E[0];
				//cout <<"k[1] " <<k[1]<<endl;
			}
			else{
				double term2 = 0;
				for (int k = 1; k <= j - 1; k++){
					term2 = term2 + a[k][j - 1] * R[j - k];
				}
				//cout <<"term2 "<< term2<<endl;
				K[j] = (R[j] - term2) / E[j - 1];
				//cout <<"k["<<i<<"]: " <<k[i]<<endl;
			}
			a[j][j] = K[j];
			//cout <<"a["<<i<<"]["<<i<<"] "<<a[i][i]<<endl;
			for (int k = 1; k <= j - 1; k++){
				a[k][j] = (a[k][j - 1] - K[j] * a[j - k][j - 1]);
				//cout <<"a["<<j<<"]["<<i<<"] "<<a[i][i]<<endl;
			}
			E[j] = (1 - pow(K[j], 2))*E[j - 1];
			//cout << "E["<<i<<"] " <<E[i]<<endl;
		}
		//the algorithm ends here
	}

	for (int i = 1; i < (p + 1); i++){
		// LPC co-efficients a[i][p]
		// cout << a[i][p] << endl;
	}

	// Computing cepstral co-efficients
	c[0] = log(double(R[0] * R[0])) / log(double(2));
	for (int m = 1; m < c.size(); m++){
		if (m <= p){
			c[m] = a[m][p];
			for (int k = 1; k <= m - 1; k++){
				c[m] += (k * c[k] * a[m - k][p]) / m;
			}
		}
		else {
			c[m] = 0;
			for (int k = 1; k <= m - 1; k++){
				c[m] += (k * c[k] * a[m - k][p]) / m;
			}
		}
	}

	// print output c[i]
	//cout << endl << "Cepstral co-efficients : " << endl;
	for (int i = 0; i < c.size(); i++){
		//cout << std::fixed << std::setprecision(18) << c[i] << endl;
	}


	for (int i = 0; i < c.size(); i++){
		//myfile << std::fixed << std::setprecision(18) << c[i] << ",";
	}
	//myfile << "\n";

	return c;
}

// returns float energy_avg per frame(320 samples per frame) of v.size() samples
double calculate_energy(vector<double> v){
	long long energy = 0;
	for (int i = 0; i < v.size(); i++){
		energy += v[i] * v[i];
	}
	return (double)((double)(energy * FRAME_SIZE) / (double)v.size());
}

// returns float zcr_avg per frame(320 samples per frame) of v.size() samples
double calculate_zcr(vector<double> v){
	long long zcr = 0;
	int y = 0;
	if (v[0] > 0){
		y = 1;
	}
	else if (v[0] < 0){
		y = -1;
	}
	for (int i = 0; i < v.size(); i++){
		if (v[i] > 0 && y > 0){
			y = 1;
		}
		else if (v[i] > 0 && y == 0){
			y = 1;
		}
		else if (v[i] > 0 && y < 0){
			zcr++;
			y = 1;
		}
		else if (v[i] == 0 && y > 0){
			y = 0;
		}
		else if (v[i] == 0 && y == 0){
			y = 0;
		}
		else if (v[i] == 0 && y < 0){
			y = 0;
		}
		else if (v[i] < 0 && y > 0){
			zcr++;
			y = -1;
		}
		else if (v[i] < 0 && y == 0){
			y = -1;
		}
		else if (v[i] < 0 && y < 0){
			y = -1;
		}
	}
	return (double)((double)(zcr * FRAME_SIZE) / (double)v.size());
}

double find_DC_coefficient(vector<double> voice){
	double sum = 0;
	for (int i = 0; i < voice.size(); i++){
		sum += voice[i];
	}
	return sum / (double)(voice.size());
}

double find_normalization_coefficient(vector<double> voice){
	double max = 0;
	for (int i = 0; i < voice.size(); i++){
		if (max < voice[i]){
			max = voice[i];
		}
	}
	return (double)NORM_CONSTANT / max;
}

vector<vector<double>> frame_generator(string filename, string noise_filename, int noise_frame_count){
	vector<vector<double>> frames;
	vector<double> voice, noise;
	double voice_dc_coeff, voice_norm_coeff;
	double noise_energy_avg, noise_zcr_avg;
	double energy_threshold = 0, zcr_threshold = 0;
	ifstream frame_stream, noise_stream;
	frame_stream.open(filename);
	noise_stream.open(noise_filename);

	// Check for errors
	if (!frame_stream){
		cout << endl << "FATAL ERROR : frame_generator :: frame_stream :: file does not exist :: " + filename << endl;
		return frames;
	}
	if (!noise_stream){
		cout << endl << "LOG ERROR : frame_generator :: noise_stream :: file does not exist :: " + noise_filename << endl;
	}

	// Generate frames
	// Read data
	int data;
	if (!noise_stream){
		// x frames of 320 samples selected as noise
		for (int i = 0; i < noise_frame_count * FRAME_SIZE; i++){
			frame_stream >> data;
			noise.push_back(data);
		}
		while (frame_stream >> data){
			voice.push_back(data);
		}
	}
	else {
		while (noise_stream >> data){
			noise.push_back(data);
		}
		while (frame_stream >> data){
			voice.push_back(data);
		}
	}
	cout << "Noise samples count = " << noise.size() << " || Voice samples count = " << voice.size() << endl;
	// Find DC coefficient and Normalization constant
	voice_dc_coeff = find_DC_coefficient(voice);
	voice_norm_coeff = find_normalization_coefficient(voice);
	cout << endl << "DC coefficient = " << voice_dc_coeff << endl;
	cout << "Normalization coefficient = " << voice_norm_coeff << " || Max value ~~ " << (int)NORM_CONSTANT / voice_norm_coeff << endl;
	for (int i = 0; i < noise.size(); i++){
		noise[i] = (double)((noise[i] - voice_dc_coeff) * voice_norm_coeff);
	}
	for (int i = 0; i < voice.size(); i++){
		voice[i] = (double)((voice[i] - voice_dc_coeff) * voice_norm_coeff);
	}

	// Find average energy and zcr of noise
	noise_energy_avg = calculate_energy(noise);
	noise_zcr_avg = calculate_zcr(noise);
	energy_threshold = THRESHOLD_CONSTANT * noise_energy_avg;
	zcr_threshold = THRESHOLD_CONSTANT * noise_zcr_avg;
	cout << endl << "Noise average energy per frame : " << (int)noise_energy_avg << endl;
	cout << endl << "Noise average ZCR per frame : " << noise_zcr_avg << endl;
	cout << endl << "Energy threshold : " << (int)energy_threshold << endl;
	cout << endl << "ZCR threshold : " << zcr_threshold << endl;

	// Generate frames
	vector<double> current_frame;
	int frame_count = 1;
	for (int i = 0; i < voice.size(); i++){
		if (current_frame.size() >= FRAME_SIZE){
			if (calculate_energy(current_frame) > energy_threshold){
				frames.push_back(current_frame);
				cout << "Frame no : " << frame_count << "  Energy over threshold. Energy = " << calculate_energy(current_frame) << " ; ZCR = " << calculate_zcr(current_frame) << endl;
			}
			else {
				if (calculate_zcr(current_frame) > zcr_threshold){
					frames.push_back(current_frame);
					cout << "Frame no : " << frame_count << "  ZCR over threshold. Energy = " << calculate_energy(current_frame) << " ; ZCR = " << calculate_zcr(current_frame) << endl;
				}
				else{
					cout << "Frame no : " << frame_count << "  Low energy and ZCR. Energy = " << calculate_energy(current_frame) << " ; ZCR = " << calculate_zcr(current_frame) << endl;
				}
			}
			frame_count++;
			current_frame.clear();
		}
		else {
			current_frame.push_back(voice[i]);
		}
	}

	return frames;
}

vector<vector<double>> hamming_window_applier(vector<vector<double>> frames){
	vector<vector<double>> new_frames;
	double w_m;
	int M = FRAME_SIZE, m;
	vector<double> temp;
	for (int i = 0; i < frames.size(); i++){
		for (int j = 0; j < frames[i].size(); j++){
			m = j % M;
			w_m = 0.54 - 0.46*cos(2 * M_PI * m / M);
			temp.push_back(frames[i][j] * w_m);
		}
		new_frames.push_back(temp);
		temp.clear();
	}
	return new_frames;
}

vector<vector<double>> cepstral_generator(vector<vector<double>> frames){
	vector<vector<double>> c, hammed_frames;

	// Apply hamming window
	cout << endl << "Applying hamming window over all voice frames . . ." << endl;
 	hammed_frames = hamming_window_applier(frames);

	cout << endl << "Generating cepstral coefficients . . ." << endl;
	for (int i = 0; i < hammed_frames.size(); i++){
		c.push_back(find_cepstral(hammed_frames[i]));
	}
	cout << "Generated cepstral coefficients." << endl << endl;
	return c;
}

int write_2d_csv_log(vector<vector<double>> data, string filename){
	ofstream logstream;
	logstream.open(filename);
	logstream.precision(30);
	if (logstream){
		for (int i = 0; i < data.size(); i++){
			for (int j = 0; j < data[i].size(); j++){
				logstream << std::fixed << std::setprecision(30) << data[i][j] << ", ";
			}
			logstream << "\n";
		}
		logstream.close();
		return 1;
	}
	else {
		return 0;
	}
}

int recognize(string filename){
	vector<vector<double>> frames, cepstral;
	string noise_filename = "noise.txt";
	int noise_frame_count = 3; // If noise_filename is not found, few initial frames of original file will be used as noise reference
	frames = frame_generator(filename, noise_filename, noise_frame_count);
	if (frames.size() < 5){
		cout << endl << "RECOGNITION FAILED :: recognize :: Minimum frame limits reached :: frame_generator output size is too low" << endl << endl;
		return 0;
	}
	cepstral = cepstral_generator(frames);
	cout << endl << "Writing cepstral coefficient to recognition_cepstral.csv . . ." << endl;
	write_2d_csv_log(cepstral, "recognition_cepstral.csv");


	return 0;
}

int _tmain(int argc, _TCHAR* argv[])
{
	cout << endl << "\t*** DIGIT RECOGNIZER ***\t" << endl;
	int t = 1;
	string recog_filename;
	while (t){
		cout << endl << "Enter 1 to RECOGNIZE !\tEnter 2 to TRAIN !\tEnter 0 to EXIT !" << endl;
		cin >> t;
		switch (t)
		{
		case 1:
			cout << endl << "Enter filename : " << endl;
			cin >> recog_filename;
			recognize(recog_filename);
			break;
		case 2:
			//train();
			break;
		case 0:
			//exit();
			break;
		default:
			t = 1;
			break;
		}
	}
	/*
	string filename1, filename2;
	char data[100];
	int x;


	cout << endl << "SAVING OUTPUT . . ." << endl;
	myfile.open(filename1);
	for (int i = 0; i <= 18; i++){
	myfile << "c[" << i << "],";
	}
	myfile << "\n";

	cout << endl << "READING VOICE DATA . . . " << endl;
	temp.clear();
	ofstream out("out.txt");
	int sample_count = 1;
	for (int i = 0; i < v.size(); i++){
	int M = 320;
	int m = i % 320;
	float w_m = 0.54 - 0.46*cos(2 * M_PI * m / M);
	temp.push_back(v[i] * w_m);
	out << v[i] * w_m << endl;
	sample_count++;
	if (sample_count == 320){
	cout << endl << "Passing frame . . . " << endl;
	find_cepstral(temp);
	temp.clear();
	sample_count = 1;
	// ----------------------------------- //
	//break;
	}
	}

	cout << endl << "OUTPUT SAVED IN CSV" << endl;
	myfile.close();


	*/

	cout << endl << endl;
	system("pause");
	return 0;
}

