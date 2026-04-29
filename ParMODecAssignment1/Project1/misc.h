#pragma once
#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED
#include <iostream>
#include <string>
#include <vector>
#include "ilcplex/ilocplex.h"
#include <omp.h> 
#include <tuple>
#include <ctime>
#include "dominatingNode.h"
#include "node.h"
using namespace std;


std::vector<dominatingNode> pareto_filter(vector<dominatingNode>&  solutionList, int num_obj)
{
	//Maximization type of objectives are assummed
	std::vector<dominatingNode>  filtered_solutionList;
	int num_dominated_obj = 0;
	for (int k = 0; k < solutionList.size();k++) {
		if (solutionList[k].dominated) { continue; }
		for (int kk = 0; kk < solutionList.size() && (kk!=k); kk++) {
			if (solutionList[kk].dominated) { continue; }
			num_dominated_obj = 0;
			for(int i =0;i<num_obj;i++){
				if (solutionList[kk].objValue[i]+0.00000001>= solutionList[k].objValue[i]) {
					num_dominated_obj += 1;
				}
				else { break; }
			}
			if (num_dominated_obj==num_obj) {
				solutionList[k].dominated = true;
				break;
			}
		}
	}
	for (int k = 0; k < solutionList.size(); k++) {
		if (!solutionList[k].dominated) {
			filtered_solutionList.push_back(solutionList[k]);
		}
	}
	return filtered_solutionList;

};
	
std::vector<dominatingNode> merge_parallels(vector<dominatingNode>&  solution_piece,vector<dominatingNode>&  solutionList, int num_obj) {
	for (int i = 0; i < solution_piece.size(); i++) { 
		solutionList.push_back(solution_piece[i]);
		//printf_s(solution_piece[i].objValue[0]);
	}
	return solutionList;
};
void log_it_date(ofstream& filename)
{
	time_t now = time(0);
	char* dt = ctime(&now);
	filename << "\n ---------------------------------- \n DATE: \t" << dt  << endl;
};
void log_it_float(const string &title,ofstream& filename, float output)
{
	filename<< title << ": \t"<< output << "\t" << endl;
};
void log_it_vector(ofstream& filename, vector<dominatingNode>&  solutionList, int num_obj)
{
	for (int i = 0; i < solutionList.size(); i++) {
		for (int j = 0; j < num_obj; j++) {
			filename << solutionList[i].objValue[j];
			if(j<num_obj-1){ filename <<  ","; }
		}
		filename << endl;
	}
};
#endif
