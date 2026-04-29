#pragma once
#include <iostream>
#include <string>
#include <vector>
#include "ilcplex/ilocplex.h"
#include <omp.h> 
#include <tuple>
using namespace std;

class dominatingNode //class definition
{
public:

	dominatingNode(IloEnv env);
	virtual ~dominatingNode();
	float objValue[3];
	int level;
	bool dominated;
	float Probability;
};

dominatingNode::dominatingNode(IloEnv env) {
	level = 0;
	//address=NULL;
	objValue[0] = 0;
	objValue[1] = 0;
	objValue[2] = 0;
	dominated = false;
	Probability = 0;
}

dominatingNode::~dominatingNode() {
	//delete[] objValue;
	//address
}
