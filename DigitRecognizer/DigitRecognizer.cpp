// DigitRecognizer.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <typeinfo>
#include <vector>
#include <math.h>
#include <random>
#include <iomanip>
#include <ios>
using namespace std;

#define MAX_DISTORTION 9999
int NORM_CONSTANT = 10000;
double THRESHOLD_CONSTANT = 10;
int FRAME_SIZE = 320;

//vector<float> v, temp, noise;
ofstream logstream;
// N = no of states in the model
int HMM_N = 5;
// M = no of distinct observation symbols per state i.e. discrete alphabet size
int HMM_M = 32;
int HMM_T = 100;

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

int write_1d_csv_log(vector<int> data, string filename){
	ofstream logstream;
	logstream.open(filename);
	logstream.precision(30);
	if (logstream){
		for (int i = 0; i < data.size(); i++){
			logstream << std::fixed << std::setprecision(30) << data[i] << ", ";
		}
		logstream.close();
		return 1;
	}
	else {
		return 0;
	}
}

//find euclidean distance between "a" and "b" vector with a.size() = b.size() = dimension
double distance(vector<double> a, vector<double> b, int dimension){
	double distance = 0;
	for (int i = 0; i < dimension; i++){
		distance += (a[i] - b[i])*(a[i] - b[i]);
	}
	return sqrt(distance);
}

//codebook = stores all the M centroids :: vector<vector<double>>
//training_set = contains all L cepstral vectors :: vector<vector<double>>
//dimension = dimension of cepstral vector inputs, generally 12 :: int
//prev_loop_distortion = avg_dist in previous loop :: double
//threshold = Stop when average distance falls below preset threshold :: double
//max_loop - loop_count = remaining iteration count after which KMeans will stop :: int
//loop_count = current loop number :: int
vector<vector<double>> KMeans(vector<vector<double>> codebook, vector<vector<double>> training_set, int dimension, double prev_loop_distortion, double threshold, int max_loop, int loop_count){
	int M = codebook.size();
	int L = training_set.size();
	double min_distance = 9999;
	double avg_dist = 0;
	vector<vector<double>> new_codebook = codebook;
	vector<vector<vector<double>>> cluster;

	//initialize clusters
	for (int i = 0; i < M; i++){
		vector<double> temp_c;
		vector<vector<double>> temp_i;
		temp_c.push_back((double)i);
		temp_i.push_back(temp_c);
		cluster.push_back(temp_i);
	}

	//find min distance and assign clusters
	for (int i = 0; i < L; i++){
		int min_codeword = 0;
		min_distance = 9999;
		for (int j = 0; j < M; j++){
			double curr_distance = distance(training_set[i], new_codebook[j], dimension);
			if (curr_distance < min_distance){
				min_distance = curr_distance;
				min_codeword = j;
			}
		}

		//assign cluster for training_set[i]
		vector<vector<double>> temp_i;
		for (int i = 0; i < cluster[min_codeword].size(); i++){
			temp_i.push_back(cluster[min_codeword][i]);
		}
		temp_i.push_back(training_set[i]);
		cluster[min_codeword] = temp_i;
	}

	//store total minimum distance for each cluster
	double total_distance = 0;
	//for each cluster find new centroid
	for (int i = 0; i < M; i++){
		vector<double> new_centroid = new_codebook[i];
		min_distance = 9999;
		//find min distance
		for (int j = 1; j < cluster[i].size(); j++){
			double curr_distance = 0;
			for (int k = 1; k < cluster[i].size(); k++){
				curr_distance += distance(cluster[i][j], cluster[i][k], dimension);
			}
			if (curr_distance < min_distance){
				min_distance = curr_distance;
				new_centroid = cluster[i][j];
			}
		}
		//assign new centroid
		new_codebook[i] = new_centroid;
		if (min_distance == 9999){
			min_distance = 0;
		}
		total_distance += min_distance;
	}

	avg_dist = total_distance / (double)training_set.size();

	// Generate Log file
	// log file should contain the various intermediate outputs like the no. of vectors that are allotted
	// to each cluster as you go through the algorithm. Also the distortion for the current codebook at that time.
	ofstream outfile;
	outfile.open("log_" + to_string(codebook.size()) + "_" + to_string(loop_count) + ".txt");
	outfile.precision(18);
	outfile << "Distortion for the current codebook : " << avg_dist << endl;
	for (int i = 0; i < M; i++){
		outfile << "Cluster " << (i + 1) << " elements : " << endl;
		for (int j = 1; j < cluster[i].size(); j++){
			for (int k = 1; k < cluster[i][j].size(); k++){
				outfile << cluster[i][j][k] << ",";
			}
			outfile << endl;
		}
		outfile << endl;

	}
	outfile.close();


	//print output
	ofstream graphfile;
	graphfile.precision(18);
	//std::cout.precision(18);
	graphfile.open("graph.csv", ios::app);
	cout << "KMeans loop -> " << loop_count << "| average distance = " << (float)avg_dist << " | threshold = " << (int)threshold << endl;
	graphfile << codebook.size() << "," << (float)avg_dist << "\n";
	graphfile.close();

	// if avg_dist is less than threshold or same as previous loop, then return
	if (avg_dist < threshold){
		cout << "AVERAGE DISTORTION IS LESS THAN THRESHOLD." << endl;
		return new_codebook;
	}
	else if (avg_dist == prev_loop_distortion){
		cout << "AVERAGE DISTORTION IS SAME AS PREVIOUS LOOP." << endl;
		return new_codebook;
	}
	else if (loop_count >= max_loop){
		cout << "MAX LOOP COUNT EXCEEDED!" << endl;
		return new_codebook;
	}
	else {
		//iterate KMeans with newly generated codebook
		new_codebook = KMeans(new_codebook, training_set, dimension, avg_dist, threshold, max_loop, loop_count + 1);
		return new_codebook;
	}
}

vector<double> split_vector(vector<double> v, double epsilon){
	vector<double> new_v;
	for (int i = 0; i < v.size(); i++){
		new_v.push_back((double)(((double)1 + epsilon)*v[i]));
	}
	return new_v;
}

//codebook = stores all the M centroids :: vector<vector<double>>
//training_set = contains all L cepstral vectors :: vector<vector<double>>
//dimension = dimension of cepstral vector inputs, generally 12 :: int
//M = Upper limit of codebook size :: int
//epsilon = splitting parameter :: double
//threshold = Stop when average distance falls below preset threshold :: double
//max_loop = max iteration count after which KMeans will stop :: int
//LBG_loop_count = LBG loop number :: int
vector<vector<double>> LBG(vector<vector<double>> codebook, vector<vector<double>> training_set, int dimension, int M, double epsilon, double threshold, int max_loop, int LBG_loop_count){
	int m = codebook.size();
	int L = training_set.size();
	double min_distance = 9999, total_distance = 0;
	double curr_distortion = 0;
	vector<vector<double>> new_codebook, split_codebook;

	for (int i = 0; i < training_set.size(); i++){
		min_distance = 9999;
		//find min distance
		for (int j = 0; j < codebook.size(); j++){
			double curr_distance = distance(training_set[i], codebook[j], dimension);
			if (curr_distance < min_distance){
				min_distance = curr_distance;
			}
		}
		//assign new centroid
		total_distance += min_distance;
	}

	curr_distortion = total_distance / (double)training_set.size();

	ofstream graphfile;
	graphfile.precision(18);
	graphfile.open("graph.csv", ios::app);
	graphfile << codebook.size() << "," << (float)curr_distortion << "\n";
	graphfile.close();

	cout << endl << "Starting LBG :: Loop " << LBG_loop_count << " :: Codebook size = " << m << endl << endl;
	new_codebook = KMeans(codebook, training_set, dimension, MAX_DISTORTION, threshold, max_loop, 1);

	// Generate codebook_[i].csv for current codebook size
	ofstream outfile;
	outfile.open("codebook_" + to_string(codebook.size()) + ".csv");
	outfile.precision(18);
	for (int i = 0; i < codebook.size(); i++){
		for (int j = 0; j < codebook[0].size(); j++){
			outfile << codebook[i][j] << ",";
		}
		outfile << "\n";
	}
	outfile.close();

	if (m >= M){
		return new_codebook;
	}

	// Split the codebook in twice size
	cout << endl << "Splitting codebook :: Current codebook size = " << m << endl;
	for (int i = 0; i < new_codebook.size(); i++){
		split_codebook.push_back(split_vector(new_codebook[i], epsilon));
		split_codebook.push_back(split_vector(new_codebook[i], (double)(-1 * epsilon)));
	}
	cout << "Codebook split :: New codebook size = " << split_codebook.size() << endl;

	// Recurse LBG for split_codebook
	new_codebook = LBG(split_codebook, training_set, dimension, M, epsilon, threshold, max_loop, LBG_loop_count + 1);

	return new_codebook;
}

vector<vector<double>> codebook_generator(vector<vector<double>> modified_cepstral){
	cout << endl << "INITIATING CODEBOOK GENERATION :: LBG" << endl;
	int M = 32; //maximum number of codebook vectors
	double threshold = (double)0.01;
	double epsilon = (double)0.03;
	int max_loop = 10;
	vector<vector<double>> training_data, init_centroid, codebook;
	training_data = modified_cepstral;


	//equally likely random centroid generation
	cout << endl << "SELECTING AN INITIAL RANDOM CENTROID . . ." << endl;
	std::random_device rd; // obtain a random number from hardware
	std::mt19937 eng(rd()); // seed the generator
	std::uniform_int_distribution<> distr(0, training_data.size() - 1); // define the range
	init_centroid.push_back(training_data[distr(eng)]);

	//k-means & LBG
	ofstream outfile, graphfile;
	outfile.open("codebook.csv");
	outfile.precision(18);
	graphfile.open("graph.csv");
	graphfile << "";
	graphfile.close();
	cout << endl << "STARTING LBG ALGORITHM . . ." << endl;
	codebook = LBG(init_centroid, training_data, training_data[0].size(), M, epsilon, threshold, max_loop, 1);
	cout << endl << "PRINTING FINAL CODEBOOK . . ." << endl << endl;
	for (int i = 0; i < codebook.size(); i++){
		for (int j = 0; j < codebook[0].size(); j++){
			//cout << codebook[i][j] << " ";
			outfile << codebook[i][j] << ",";
		}
		//cout << endl;
		outfile << "\n";
	}
	outfile.close();


	return codebook;
}

// Generate observation sequence with states 1, 2, 3, .... , 32
vector<int> observation_generator(vector<vector<double>> codebook, vector<vector<double>> cepstral){
	vector<int> obs;
	int min_i = 0;
	double min = 9999, temp_dist;
	for (int i = 0; i < cepstral.size(); i++){
		min_i = 0;
		min = 9999;
		for (int j = 0; j < codebook.size(); j++){
			temp_dist = distance(cepstral[i], codebook[j], cepstral[i].size());
			if (temp_dist < min){
				min = temp_dist;
				min_i = (j + 1);
			}
		}
		obs.push_back(min_i);
	}
	return obs;
}

vector<vector<double>> forward_procedure(vector<int> O1, vector<vector<double>> A, vector<vector<double>> B, vector<double> PI){
	int T1 = O1.size();
	int N = HMM_N, M = HMM_M;
	vector<vector<double>> alpha(T1, vector<double>(N, 1));
	// Initialize
	//cout << "Initializing..." << endl;
	for (int i = 0; i < N; i++){
		alpha[0][i] = (double)((double)PI[i] * (double)B[i][O1[0] - 1]);
		//cout << "B[" << i << "][" << O1[0] - 1 << "]  ";
		//cout << PI[i] << " * " << B[i][O1[0] - 1] << " = " << alpha[0][i] << endl;
	}

	// Induction
	//cout << "Induction..." << endl;
	for (int t = 0; t < T1 - 1; t++){
		for (int j = 0; j < N; j++){
			alpha[t + 1][j] = 0.0;
			for (int i = 0; i < N; i++){
				alpha[t + 1][j] += (double)(alpha[t][i] * A[i][j]);
				//cout << alpha[t+1][i] << "*" << A[i][j] << endl;
			}
			//cout << j << " " << t + 1 << " " << O1[t + 1] << " " << B[j][O1[t + 1] - 1] << endl;
			//if (alpha[t + 1][j] > 0) cout << alpha[t + 1][j] <<" "<< B[j][O1[t + 1] - 1] << endl;
			alpha[t + 1][j] = alpha[t + 1][j] * B[j][O1[t + 1] - 1];
			//cout << alpha[t + 1][j] << endl;
		}
	}

	// Termination
	//cout << "Termination..." << endl;
	double P_O_given_lambda = 0;
	for (int i = 0; i < N; i++){
		P_O_given_lambda += alpha[T1 - 1][i];
	}

	//cout << "P(O|lambda) = " << P_O_given_lambda << endl;
	logstream << "P(O|lambda) = " << P_O_given_lambda << endl;
	return alpha;
}

vector<vector<double>> backward_procedure(vector<int> O, vector<vector<double>> A, vector<vector<double>> B, vector<double> PI){
	int T = O.size();
	int N = HMM_N, M = HMM_M;
	vector<vector<double>> beta(T, vector<double>(N, 1));

	// Initialization
	//cout << "Initialization..." << endl;
	for (int i = 0; i < N; i++){
		beta[T - 1][i] = 1;
	}

	// Induction
	//cout << "Induction..." << endl;
	for (int t = T - 2; t >= 0; t--){
		for (int i = 0; i < N; i++){
			beta[t][i] = 0;
			for (int j = 0; j < N; j++){
				beta[t][i] += (A[i][j] * B[j][O[t + 1] - 1] * beta[t + 1][j]);
			}
		}
	}

	// Termination
	//cout << "Termination..." << endl;
	double P_O_given_lambda = 0;
	for (int i = 0; i < N; i++){
		P_O_given_lambda += beta[0][i] * PI[i] * B[i][O[0] - 1];
	}

	//cout << "P(O|lambda) = " << P_O_given_lambda << endl;

	return beta;
}

double viterbi_algorithm(vector<int> O, vector<vector<double>> A, vector<vector<double>> B, vector<double> PI){
	int N = HMM_N, M = HMM_M;
	int T = O.size();
	vector<vector<double>> delta(T, vector<double>(N, 1));
	vector<vector<double>> shi(T, vector<double>(N, 1));

	// Initialization
	//cout << "Initialization..." << endl;
	for (int i = 0; i < N; i++){
		delta[0][i] = PI[i] * B[i][O[0] - 1];
		shi[0][i] = 0;
	}

	// Recursion
	//cout << "Recursion..." << endl;
	double max = -1;
	int max_i = -1;
	for (int t = 1; t < T; t++){
		for (int j = 0; j < N; j++){
			max = delta[t - 1][0] * A[0][j];
			max_i = 0;
			for (int i = 0; i < N; i++){
				if (delta[t - 1][i] * A[i][j] > max){
					max = delta[t - 1][i] * A[i][j];
					max_i = i;
				}
			}
			delta[t][j] = max * B[j][O[t] - 1];
			shi[t][j] = max_i;
		}
	}

	// Termination
	//cout << "Termination..." << endl;
	double P_star = -1;
	vector<int> q_star(T, -1);
	for (int i = 0; i < N; i++){
		if (delta[T - 1][i] > P_star){
			P_star = delta[T - 1][i];
			q_star[T - 1] = i;
		}
	}

	// Path Backtracking
	//cout << "Path Backtracking..." << endl;
	for (int t = T - 2; t >= 0; t--){
		q_star[t] = (int)shi[t + 1][q_star[t + 1]];
	}

	// Output
	//cout << "Output..." << endl;
	//cout << "P_star : " << P_star << endl;
	/*cout << "Q_star :" << endl;
	for (int i = 0; i < T; i++){
		cout << q_star[i] << " ";
	}
	cout << endl;*/
	

	// Return delta and phi
	//vector<vector<vector<double>>> delta_and_shi;
	//delta_and_shi.push_back(delta);
	//delta_and_shi.push_back(shi);
	//return delta_and_shi;
	return P_star;
}

vector<vector<vector<double>>> baum_welch(vector<int> O, vector<vector<double>> A, vector<vector<double>> B, vector<double> PI, vector<vector<double>> alpha, vector<vector<double>> beta){
	int N = HMM_N, M = HMM_M;
	int T = O.size();
	//cout << "Calculating Xi..." << endl;
	vector<vector<vector<double>>> Xi(T, vector<vector<double>>(N, vector<double>(N, 0)));
	for (int t = 0; t < T - 1; t++){
		double sum = 0;
		for (int i = 0; i < N; i++){
			for (int j = 0; j < N; j++){
				sum += (alpha[t][i] * A[i][j] * B[j][O[t + 1] - 1] * beta[t + 1][j]);
			}
		}
		for (int i = 0; i < N; i++){
			for (int j = 0; j < N; j++){
				Xi[t][i][j] = (alpha[t][i] * A[i][j] * B[j][O[t + 1] - 1] * beta[t + 1][j]) / sum;
			}
		}
	}

	return Xi;
}

int HMM_optimization(vector<int> O, vector<vector<double>> A, vector<vector<double>> B, vector<double> PI, string extension){
	int N = HMM_N, M = HMM_M;
	int T = O.size();
	vector<vector<double>> A_curr = A, B_curr = B;
	logstream.open("HMM_P_star" + extension + ".txt");
	vector<double> PI_curr = PI;
	double P_star_diff = 1, prev_P_star = -1, P_star = 0;
	int count = 1;

	while (P_star_diff > 0){
		// The Forward Procedure
		//cout << endl << "Initiate Hidden Markov Modelling . . . " << endl;
		//cout << endl << "Probability evaluation :: Forward procedure . . . " << endl << endl;
		vector<vector<double>> alpha;
		alpha = forward_procedure(O, A_curr, B_curr, PI_curr);

		// The Backward Procedure
		//cout << endl << "Probability evaluation :: Backward procedure . . . " << endl << endl;
		vector<vector<double>> beta;
		beta = backward_procedure(O, A_curr, B_curr, PI_curr);

		// Viterbi algo
		//cout << endl << "Viterbi algorithm :: Looking for optimal state sequence . . . " << endl << endl;
		/*vector<vector<vector<double>>> delta_n_shi;
		vector<vector<double>> delta, shi;
		delta_n_shi = viterbi_algorithm(O, A, B, PI);
		delta = delta_n_shi[0];
		shi = delta_n_shi[1];*/
		P_star = viterbi_algorithm(O, A_curr, B_curr, PI_curr);

		// Parameter estimation
		//cout << endl << "Baum-Welch Method/Expectation-Maximization Method :: Estimating parameters . . . " << endl << endl;
		vector<vector<vector<double>>> Xi;
		Xi = baum_welch(O, A_curr, B_curr, PI_curr, alpha, beta);

		//cout << "Calculating Gamma..." << endl;
		vector<vector<double>> gamma(T, vector<double>(N, 0));
		for (int t = 0; t < T; t++){
			for (int i = 0; i < N; i++){
				gamma[t][i] = 0;
				for (int j = 0; j < N; j++){
					gamma[t][i] += Xi[t][i][j];
				}
			}
		}

		//cout << "Calculating A_bar..." << endl;
		vector<vector<double>> A_bar(N, vector<double>(N, 0));
		for (int i = 0; i < N; i++){
			for (int j = 0; j < N; j++){
				double sum1 = 0, sum2 = 0;
				for (int t = 0; t < T - 1; t++){
					sum1 += Xi[t][i][j];
					sum2 += gamma[t][i];
				}
				A_bar[i][j] = sum1 / sum2;
			}
		}


		//cout << "Calculating B_bar..." << endl;
		vector<vector<double>> B_bar(N, vector<double>(M, 0));
		for (int j = 0; j < N; j++){
			for (int k = 0; k < M; k++){
				double sum1 = 0, sum2 = 0;
				for (int t = 0; t < T - 1; t++){
					if (O[t] - 1 == k){
						sum1 += gamma[t][j];
					}
					sum2 += gamma[t][j];
				}
				B_bar[j][k] = sum1 / sum2;
			}
		}

		//cout << "Calculating PI_bar..." << endl;
		vector<double> PI_bar;
		for (int i = 0; i < N; i++){
			PI_bar.push_back(gamma[0][i]);
		}

		// Find new P_O_given_lambda
		//cout << endl << "Finding optimized probability..." << endl;
		forward_procedure(O, A_bar, B_bar, PI_bar);

		A_curr = A_bar;
		B_curr = B_bar;
		PI_curr = PI_bar;
		P_star_diff = P_star - prev_P_star;
		cout << endl << "Loop " << count << " : P* = " << P_star << " ; P*Diff = " << P_star_diff << endl;
		logstream << "Loop " << count << " : P* = " << P_star << " ; P*Diff = " << P_star_diff << endl << endl;
		prev_P_star = P_star;
		count++;
	}

	/* Writing in file
	cout << endl << "Creating new a" + extension + ".txt b" + extension + ".txt pi" + extension + ".txt files..." << endl;
	string filename_A_bar = "a" + extension + ".txt";
	string filename_B_bar = "b" + extension + ".txt";
	string filename_PI_bar = "pi" + extension + ".txt";

	ofstream a_bar_stream, b_bar_stream, pi_bar_stream;

	a_bar_stream.open(filename_A_bar);
	a_bar_stream.precision(50);
	for (int i = 0; i < N; i++){
	for (int j = 0; j < N; j++){
	a_bar_stream << A_bar[i][j] << " ";
	}
	a_bar_stream << "\n";
	}
	a_bar_stream.close();

	b_bar_stream.open(filename_B_bar);
	b_bar_stream.precision(50);
	for (int j = 0; j < N; j++){
	for (int k = 0; k < M; k++){
	b_bar_stream << B_bar[j][k] << " ";
	}
	b_bar_stream << "\n";
	}
	b_bar_stream.close();

	pi_bar_stream.open(filename_PI_bar);
	pi_bar_stream.precision(50);
	for (int i = 0; i < N; i++){
	pi_bar_stream << PI_bar[i] << " ";
	}
	pi_bar_stream.close();

	cout << endl << "Write complete !" << endl << endl;
	*/

	logstream.close();
	return 1;
}

int HMM(vector<int> O, string extension){
	int N = HMM_N, M = HMM_M, T = O.size();
	string filename_A = "defaultHMM_A.txt";
	string filename_B = "defaultHMM_B.txt";
	string filename_PI = "defaultHMM_PI.txt";

	vector<vector<double>> A, B;
	vector<double> PI;

	ifstream a_stream;
	a_stream.open(filename_A);
	a_stream.precision(500);
	double a_val;
	for (int i = 0; i < N; i++){
		vector<double> temp;
		for (int j = 0; j < N; j++){
			a_stream >> a_val;
			temp.push_back(a_val);
		}
		A.push_back(temp);
	}
	a_stream.close();

	ifstream b_stream;
	b_stream.open(filename_B);
	b_stream.precision(500);
	double b_val;
	for (int i = 0; i < N; i++){
		vector<double> temp;
		for (int j = 0; j < M; j++){
			b_stream >> b_val;
			temp.push_back(b_val);
		}
		B.push_back(temp);
	}
	b_stream.close();

	ifstream pi_stream;
	pi_stream.open(filename_PI);
	pi_stream.precision(500);
	double pi_val;
	for (int i = 0; i < N; i++){
		pi_stream >> pi_val;
		PI.push_back(pi_val);
	}
	pi_stream.close();

	// HMM starts here
	HMM_optimization(O, A, B, PI, extension);
	// HMM ends here

	return 0;
}

int recognize(string filename){
	vector<vector<double>> frames, cepstral, codebook;
	vector<int> observation_sequence;
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
	cout << endl << "OUTPUT SAVED IN CSV" << endl << endl;

	vector<vector<double>> modified_cepstral;
	// cepstral contains c[0] to c[18]; remove c[0], c[13], ... , c[18]
	vector<double> modified_c;
	for (int i = 0; i < cepstral.size(); i++){
		//refine ci vector
		for (int j = 1; j <= 12; j++){
			modified_c.push_back(cepstral[i][j]);
		}
		modified_cepstral.push_back(modified_c);
		modified_c.clear();
	}
	cout << endl << "Generated training data for LBG." << endl;
	write_2d_csv_log(modified_cepstral, "recog_training_data.csv");

	// Generate codebook
	codebook = codebook_generator(modified_cepstral);

	// Generate observation sequence
	cout << endl << "Generating Observation Sequence" << endl;
	observation_sequence = observation_generator(codebook, modified_cepstral);
	cout << endl << "Generated Observation Sequence" << endl;
	write_1d_csv_log(observation_sequence, "recog_observation_seq.csv");
	cout << endl << "Observation Sequence Saved in .csv" << endl;

	// Run HMM
	cout << endl << "Initiate Hidden Markov Modelling . . . " << endl;
	HMM(observation_sequence, "_recog");

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


	cout << endl << endl;
	system("pause");
	return 0;
}

