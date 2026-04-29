#pragma once
#include <iostream>
#include <string>
#include <vector>
#include "ilcplex/ilocplex.h"
#include <omp.h> 
#include <tuple>
using namespace std;

#include "dominatingNode.h"

class node //class definition
{
public:

	node(IloEnv envMBC);
	virtual ~node();
	float *lb;
	float *ub;
	bool desofDC;
	bool DC;
	bool dominated;
	int level;
	int code;
	node* next;
	node* child[3];
	node* previous;
	node* parent;
	node* SONL;
	node* DCsibling;
	float objValue[3];
	float probability;
	int solRank;
	std::vector<node*>  dominantsListUpperLevel;
	std::vector<node*>  weakListUpperLevel;
	std::vector<int>  dominantsList;//list of solutions that can dominate the current solution (prev name dominationList)
	std::vector<int>  weakList;//list of solutions that can be dominated (prev name dominateList)

	//member function
	node copy(IloEnv envMBC);
};

node::node(IloEnv envMBC) {
	objValue[0] = IloInfinity;
	objValue[1] = IloInfinity;
	objValue[2] = IloInfinity;
	lb = new float[3];
	ub = new float[3];
	desofDC = false;
	DC = false;
	dominated = false;
	level = 0;
	code = 1;
	probability = 0;
	next = NULL;
	child[0] = NULL;
	child[1] = NULL;
	child[2] = NULL;
	previous = NULL;
	parent = NULL;
	SONL = NULL;
	DCsibling = NULL;
	solRank = IloInfinity;
};

node node::copy(IloEnv envMBC) {
	node *node2 = NULL;
	node2=new node(envMBC);
	node2->lb = new float[3];
	node2->ub = new float[3];
	//there wont be any child info copying b/c it is assumed that list returned by root_node has no children info
	for (int i = 0; i < 3; i++) {
		//pchild[i] = new node(envMBC);
		node2->lb[i] = this->lb[i];
		node2->ub[i] = this->ub[i];
		node2->objValue[i] = this->objValue[i];
	}
	node2->desofDC = this->desofDC;
	node2->DC = this->DC;
	node2->dominated = this->dominated;
	node2->level = this->level;
	node2->code = this->code;
	node2->previous = this->previous;
	node2->parent = this->parent;
	node2->SONL = this->SONL;
	node2->DCsibling = this->DCsibling;
	node2->probability = this->probability;
	node2->solRank = this->solRank;

	//std::vector<node*>  dominantsListUpperLevel;
	//std::vector<node*>  weakListUpperLevel;
	for (int i=0; i < this->dominantsList.size(); i++) {
		node2->dominantsList[i] =this->dominantsList[i];
	}
	for (int i=0; i < this->weakList.size(); i++) {
		node2->weakList[i] = this->weakList[i];
	}

	return *node2;
};
node::~node() {
	dominantsList.empty();
	weakList.empty();
	/*for(int i=0;i<dominationListUpperLevel.size();i++){
	if(dominationListUpperLevel[i]){
	delete dominationListUpperLevel[i];
	}
	}*/
	dominantsListUpperLevel.empty();
	weakListUpperLevel.empty();
};

