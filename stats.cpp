/**
 *  Project:
 *
 *  File: stats.cpp
 *  Created: Oct 08, 2013
 *  Modified: Tue 08 Oct 2013 12:47:02 PM PDT
 *
 *  Author: Abhinav Sarje <asarje@lbl.gov>
 */

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <cstring>
#include <cstdlib>
#include <limits>

const int HIST_SIZE = 100;

bool read_param_data(std::ifstream& f, const char* param, std::vector<double>& data) {
	f.seekg(0);			// go to the beginning
	char buffer[101];	// buffer to read in a line
	while(f.good()) {
		f >> buffer;
		//std::cout << "'" << buffer << "'" << std::endl;
		if(strcmp(buffer, param) == 0) break;
	} // while
	if(!f.good()) return false;
	f >> buffer;		// ignore the '='
	while(f.good()) {
		double temp;
		f >> temp;		// read the number
		data.push_back(temp);
		f >> buffer;	// read the comma or semi-colon
		if(strcmp(buffer, ";") == 0) break;
	} // while
	return true;
} // real_param_data()

bool min_max(const std::vector<double>& data, double& min_val, double& max_val) {
	min_val = std::numeric_limits<double>::max();
	max_val = std::numeric_limits<double>::min();
	for(std::vector<double>::const_iterator i = data.begin(); i != data.end(); ++ i) {
		min_val = (min_val > *i) ? *i : min_val;
		max_val = (max_val < *i) ? *i : max_val;
	} // for
	//std::cout << "MIN: " << min_val << ", MAX: " << max_val << std::endl;
	return true;
} // min_max()

bool build_histogram(const std::vector<double>& data, double min_val, double max_val,
						std::map<long int, long int>& hist) {
	double range = max_val - min_val;
	double step = (range / HIST_SIZE);
	long int min_bucket = (min_val);
	long int max_bucket = min_bucket + step * HIST_SIZE;
	//std::cout << "MIN: " << min_bucket << ", MAX: " << max_bucket << std::endl;
	for(long int i = min_bucket; i <= max_bucket; i += step) hist[i] = 0;
	for(std::vector<double>::const_iterator i = data.begin(); i != data.end(); ++ i) {
		long int bucket = min_bucket;
		while(bucket < max_bucket) {
			if(*i > bucket && *i < bucket + step) {
				++ hist[bucket];
				break;
			} // if
			bucket += step;
		} // while
		if(bucket == max_bucket) ++ hist[max_bucket];
	} // for
	//std::cout << hist.size() << std::endl;
	return true;
} // build_histogram()

bool print_histogram(const std::map<long int, long int>& hist) {
	for(std::map<long int, long int>::const_iterator i = hist.begin(); i != hist.end(); ++ i) {
		std::cout << (*i).first << "\t" << (*i).second << std::endl;
	} // for
	return true;
} // print_histogram()

bool read_and_build(std::ifstream& f, const char* param, std::map<long int, long int>& hist) {
	f.seekg(0);			// go to the beginning
	char buffer[101];	// buffer to read in a line
	while(f.good()) {
		f >> buffer;
		//std::cout << "'" << buffer << "'" << std::endl;
		if(strcmp(buffer, param) == 0) break;
	} // while
	if(!f.good()) return false;
	f >> buffer;				// ignore the '='
	while(f.good()) {
		double temp;
		f >> temp;				// read the number
		long int key = (long int) temp;	// take the floor
		if(hist.count(key) == 0) hist[key] = 1;
		else ++ hist[key];
		f >> buffer;	// read the comma or semi-colon
		if(strcmp(buffer, ";") == 0) break;
	} // while
	return true;
} // read_and_build()

int main(int narg, char** args) {
	if(narg != 3) {
		std::cout << "usage: stats <ncdump_filename> <variable_name>" << std::endl;
		return 0;
	} // if
	std::ifstream f(args[1]);

	//std::vector<double> data;
	//read_param_data(f, args[2], data);

	//std::cout << "## " << args[2] << ": " << data.size() << std::endl;
	//double min_val, max_val;
	//min_max(data, min_val, max_val);

	std::map<long int, long int> histogram;
	//build_histogram(data, min_val, max_val, histogram);
	read_and_build(f, args[2], histogram);
	print_histogram(histogram);

	f.close();

	return 0;
} // main()
