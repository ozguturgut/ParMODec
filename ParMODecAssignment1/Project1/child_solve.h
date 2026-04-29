#pragma once
#include <iostream>
#include <string>
#include <vector>
#include "ilcplex/ilocplex.h"
#include <omp.h> 
#include <tuple>
#include <time.h>
#include <fstream>
#include "dominatingNode.h"
#include "node.h"
#include "misc.h"
using namespace std;
#define NMax (4)

class child_solve
{
public:
	class child_solve
	(int id)
		: Id(id), model(env), model2(env), cplex(env), var(env), rng(env), lin_expr(env), lin_expr2(env), linCons2(env), linCons(env),obj(env)
	{
		cout << "Root for one parallel thread %d is initialized", id;

	}
	std::tuple<bool, std::vector<dominatingNode>, int> solve_main(int id, const string &body, const string &filename,
		node *root, node *root_sibling, int num_obj, vector<dominatingNode>&  solutionList,
		ofstream& BandBDataTime, ofstream& BandBDataProbability,
		int numVar, int instance)
	{
		loc_solutionList = solutionList;
		main_dc_solutions = solutionList;
		debug_mode = true;
		if (debug_mode){
			std::string outputfile1 =  "Solutions"+ std::to_string(num_obj) +"_"+ std::to_string(numVar) + "_" + std::to_string(instance) +"_thread_" + std::to_string(id) + ".xls";
			std::string outputfile2 = "LocTime" + std::to_string(num_obj) + "_" + std::to_string(numVar) + "_" + std::to_string(instance) + "_thread_" + std::to_string(id) + ".xls";
			BandBData.open(outputfile1, ios::out | ios::app);
			LocalTimeData.open(outputfile2, ios::out | ios::app);
		
			BandBData << "Available solutions at the root" <<id<< endl;
			log_it_vector(BandBData, loc_solutionList, num_obj);
			BandBData << "CODE" << "\t";
			for (int ii = 0; ii < num_obj; ii++) {
				BandBData << "lower bound for obj" <<ii<< "\t";
				BandBData << "obj value for obj" << ii <<  "\t";
				BandBData << "upper bound for obj" << ii << "\t";

			}
			BandBData << endl;
		}
		//printf_s("Lower bound of root %d's third objective %f\n", id, root->lb[2]);
		root_loc = root;
		root_sibling_loc = root_sibling;
		//printf_s("Lower bound of root of sibling's third objective %f\n", root_loc->lb[2]);
		bool tree_traversed = 0;
		int root_ready = 0;
		bool root_solved = 0;
		start_child = clock();
		if (debug_mode) {
			if (id == 3) {
				printf_s("Entered thread 2 \n");
			}
			LocalTimeData << "Started thread#" << std::to_string(id) << " is " << start_child;
		}
		numModels = 0;
		prepare_model(id, num_obj, body, filename);
		root_ready = prep_parallel(num_obj, id, body, filename, BandBDataTime, BandBDataProbability);
		node *current = NULL;
		current = new node(env);
		current = root_loc;
		if (root_ready) {
			tree_traversed=tree(model2, linCons2,current, num_obj, id, body, filename, numActNodes, BandBDataTime, BandBDataProbability);
			//node *current, int numObj, int id, const string &body, const string &filename
		}//end of if root solved
		else { tree_traversed = true; }
		end_tree_traverse = clock();
		
		BandBDataTime << std::to_string(id) << "," << "Total thread time" << "," << end_tree_traverse - start_child << " , Fathommed:" << fathommed << endl;
		if (debug_mode) {
			LocalTimeData << "Finished thread #" << std::to_string(id) << " in " << end_tree_traverse - start_child << endl;
		}
		return std::make_tuple(tree_traversed, loc_solutionList, numModels);
	}
	//################################################################################# MODEL READ FROM LP WILL BE RE-ORGANIZED ACCORDING TO THE P-1 FULL SPLIT MODEL OBJECT ############################
	void prepare_model(int th, int numobj, const string &body, const string &filename)
	{
		IloExtractable extr;
		int num_constraints = 0;
		int j = 0;
		IloNum epsilon = 0.0001;
		std::string filelocation(body + filename + ".lp");
		cplex.importModel(model, filelocation.c_str(), obj, var, rng);
		cplex.extract(model);
		nrows = cplex.getNrows();
		for (IloModel::Iterator it(model); it.ok(); ++it) {
			extr = it.operator *();
			//feasible region
			if (extr.isConstraint() == true) {
				if (num_constraints >= (nrows - numobj)) {//set of objective functions
					IloRangeI* impl = dynamic_cast<IloRangeI *>((*it).asConstraint().getImpl());
					if (impl) {
						IloRange rng(impl);
						IloExpr expr = rng.getExpr();
						//lin_expr.clear();
						for (IloExpr::LinearIterator it2 = expr.getLinearIterator(); it2.ok(); ++it2) {
							//lin_expr += it2.getCoef() * it2.getVar();
							//varVector.push_back(it2.getVar());
							obj_coeff_vector.push_back(it2.getCoef());
						}
					}
					//objList[j] = lin_expr;
					j += 1;
				}
				else {//constraints
					linCons2.add(extr.asConstraint());
					num_constraints += 1;

				}
			}
		}
		IloExpr xpr_main(env);
		for (int w = (numobj - 1); w >= 0; w--)
		{
			
			IloExpr xpr(env);
			
			if (w == 0) {
				for (int var_no = 0; var_no < var.getSize(); var_no++) {
					lin_expr += obj_coeff_vector[(w*var.getSize() + var_no)] * var[var_no];
					xpr_main += obj_coeff_vector[(w*var.getSize() + var_no)] * var[var_no];
				}
			}
			else {
				
				//printf_s("Coefficients read for child\n ");
				for (int var_no = 0; var_no < var.getSize(); var_no++) {
					//printf_s("Coefficient of w %d : %f\n", w, obj_coeff_vector[(w*var.getSize() + var_no)]);
					xpr += obj_coeff_vector[(w*var.getSize() + var_no)] * var[var_no];
					lin_expr += epsilon*(obj_coeff_vector[(w*var.getSize() + var_no)] * var[var_no]);
				}

				list_objs.push_back(xpr);
			}
		}
		list_objs.push_back(lin_expr); //main objective is at the end of this local vector
		list_objs.push_back(xpr_main); //pure form of main objective is at the end of this local vector
		
	}
	//################################################################################# SOLVE CURRENT NODE ############################################################################
	int Solve(int th, int num_obj, node *current, const string &body, const string &filename)
	{
		try
		{
			
			bool result = 1;
			//return objective values as well
			IloExtractable extr;
			vector<IloIntVar> varVector;
			
			
			int num_constraints = 0, status;
			int j = 0;
			IloNum epsilon = 0.0001;
			float optimal_value;
			std::string filelocation(body + filename + "_InsideTree.lp");
			
			model2.remove(linCons);
			linCons.endElements();
			for (int j = 0; j < num_obj; j++) {

				if (j<num_obj - 1)
				{
					linCons.add(list_objs[j] >= current->lb[num_obj - (j + 1)] + 0.01); //upper value of node
					linCons.add(list_objs[j] <= current->ub[num_obj - (j + 1)]); //upper value of node
				}
				else {
					linCons.add(list_objs[num_obj] <= current->ub[num_obj - (j + 1)]); //upper value of node
				}
			}

			obj.setExpr(list_objs[num_obj-1]);  //default  main objective is the FIRST objective in LP
			obj.setSense(IloObjective::Maximize);  //default type of objective is MAXIMIZATION type
			model2.add(obj);
			model2.add(linCons2);
			model2.add(linCons);
			cplex.extract(model2);

			cplex.exportModel(filelocation.c_str());

			cplex.setOut(envMOBB.getNullStream());
			cplex.setWarning(envMOBB.getNullStream());
			//SET THE # of threads to 1 HERE
			cplex.setParam(IloCplex::ParallelMode, 1);
			cplex.setParam(IloCplex::Param::Threads, 1);
			result = cplex.solve();
			status = cplex.getCplexStatus();
			numModels += 1;
			return status;
		}
		catch (IloException &ex)
		{
			cout << "Got exception: ";
			ex.print(cout);
			cout << endl;
			return false;
		}
		catch (...)
		{
			return false;
		}
		return true;
	}
	//############################################################# CALCULATE OBJECTIVE VALUES AND WRITE AS OUTPUT #####################################################################
	int Solve_Get_Results(int th, int num_obj, node *current, const string &body, const string &filename, ofstream& BandBData, ofstream& BandBDataTime, ofstream& BandBDataProbability)
	{
		IloNumArray	XValue(env, var.getSize());
		float dummySum = 0;
		int root_ok;
		root_ok = Solve(th, num_obj, current, body, filename);
		if (root_ok == 1 || root_ok == 101 || root_ok == 102) {
			dominatingNode tempDominatingNode(env);
			cplex.getValues(XValue, var);
			for (int w = (num_obj - 1); w >= 0; w--)
			{
				dummySum = 0;
				for (int var_no = 0; var_no < var.getSize(); var_no++) {
					dummySum += obj_coeff_vector[(w*var.getSize() + var_no)] * XValue[var_no];
					//printf_s("\n At TREE Coefficient of w %d : %f\n", w, obj_coeff_vector[(w*var.getSize() + var_no)]);
					//printf_s("At TREE value of binary %d : %f\n", w, XValue[var_no]);
				}
				if (debug_mode) {
					printf_s("Objective value of obj %d DC %f\n", w, dummySum);
				}
				tempDominatingNode.objValue[w] = dummySum;
				current->objValue[w] = dummySum;
			}
			if (debug_mode) {
				BandBData << current->code << "\t";
				for (int ii = 0; ii < num_obj; ii++) {
					BandBData << current->lb[ii] << "\t";
					BandBData << current->objValue[ii] << "\t";
					BandBData << current->ub[ii] << "\t";

				}
				BandBData << endl;
			}
			loc_solutionList.push_back(tempDominatingNode);
		}
		return root_ok;
	}
	//########################################################################DETAILS FOR SINGLE NODE###################################################################################
	int Create_Single_Node(int th, int num_obj, const string &body, const string &filename, ofstream& BandBData, ofstream& BandBDataTime, ofstream& BandBDataProbability)
	{	
		node *temp1 = NULL;
		float dummySum = 0;
		int root_ok;
		int num_children = 2 ^ (num_obj - 1) - 1;
		if (th==4){
			temp1 = root_sibling_loc;
		}else{ 
			temp1 = root_loc; }

		temp1->level = 1; //for the root node
		root_ok = Solve_Get_Results(th, num_obj, temp1, body, filename, BandBData, BandBDataTime, BandBDataProbability);
		if (root_ok==1 || root_ok == 101 || root_ok == 102) {

			temp1->solRank = loc_solutionList.size() - 1;
			//specifics for root node
			temp1->previous = NULL;
			temp1->DC = false;
			temp1->probability = 0;
			temp1->code = numActNodes;

			//create the children and replace the info inside 
			node *pchild[3] = { NULL };
			for (int j = 0; j<num_obj; j++) {
				pchild[j] = new node(env);
				temp1->child[j] = pchild[j];
			}
			temp1->level = 1; //for the root node
			pchild[0]->next = pchild[1];
			pchild[1]->next = pchild[2];
			for (int j = 0; j<num_obj; j++) {
				pchild[j]->parent = temp1;
				pchild[j]->level = (temp1->level) + 1;
			}
			//specifics for root node
			temp1->SONL = pchild[0];
			temp1->next = pchild[0];
			pchild[0]->previous = temp1;
			temp1->probability = 0;


			for (int ii = 0; ii < num_children; ii++) {
				for (int j = 0; j<num_obj; j++) {
					pchild[ii]->lb[j] = temp1->lb[j];             // store lower/upper bounds of objectives. It will be inherited to the child nodes
					pchild[ii]->ub[j] = temp1->ub[j];
				}
			}

			IloInt childCount, k1, k2;
			childCount = 0;
			for (k1 = 1; k1<3; k1++) {
				for (k2 = 1; k2<3; k2++) {
					if (k1 == 1 && k2 == 1) {
						pchild[childCount]->lb[1] = temp1->objValue[1];
						pchild[childCount]->lb[2] = temp1->objValue[2];
						numActNodes += 1;
						pchild[childCount]->code = numActNodes;
						pchild[childCount]->DC = true;
						childCount += 1;
					}
					else if (k1 == 1 && k2 == 2) {
						pchild[childCount]->lb[1] = temp1->objValue[1];
						pchild[childCount]->ub[2] = temp1->objValue[2];
						numActNodes += 1;
						pchild[childCount]->code = numActNodes;
						pchild[childCount]->DCsibling = pchild[0];
						childCount += 1;
					}
					else if (k1 == 2 && k2 == 1) {
						pchild[childCount]->ub[1] = temp1->objValue[1];
						pchild[childCount]->lb[2] = temp1->objValue[2];
						numActNodes += 1;
						pchild[childCount]->code = numActNodes;
						pchild[childCount]->DCsibling = pchild[0];
						childCount += 1;
					}
					else { continue; }
				}
			}
			loc_solutionList[loc_solutionList.size() - 1].level = temp1->level;
			loc_solutionList[loc_solutionList.size() - 1].dominated = temp1->dominated;
			numPareto += 1;
			//root_loc = temp1;
		}
		return root_ok;
	}
	//########################################################### PREPARE PARALLEL ROOTS BASED ONE OF TWO TYPES ########################################################
	int prep_parallel(int num_obj, int type, const string &body, const string &filename, ofstream& BandBDataTime, ofstream& BandBDataProbability) {
		//printf_s("Lower bound of DC grand child's third objective %f\n", root_loc->lb[2]);
		int num_children = 2 ^ (num_obj - 1) - 1;
		float dummySum = 0;
		int root_initial_ok = 0, sibling_protocol=0, root_initial_ok2 = 0;
		node *temp1 = NULL;
		//Type 1  single 
		//if (num_obj == 3 && type != 1) {
			bool root_dominated = 0;
			int dominated_axis = 0;
			for (int k = 0; k < loc_solutionList.size(); k++) {
				dominated_axis = 0;
				for (int j = 0; j < num_obj; j++) {
					if (loc_solutionList[k].objValue[j] >= root_loc->ub[j]) { dominated_axis += 1; }
				}
				if (dominated_axis == num_obj) { 
					root_dominated = 1;
					break;
				}
			}
			if(!root_dominated){
				numActNodes = 1;
				root_initial_ok = Create_Single_Node(type, num_obj, body, filename, BandBData, BandBDataTime, BandBDataProbability);
				if(root_initial_ok == 1 || root_initial_ok == 101 || root_initial_ok == 102){ root_loc = root_loc->next; }
				else { return 0; }
				
				}
		//}
		//Type 2-non DC siblings is removed for three objectives

		return 1;
	};
	//################################################################################# MAIN TREE TRAVERSAL ############################################################################
	int tree(IloModel model3, IloConstraintArray linCons,node *current, int numObj, int id, const string &body, const string &filename, int loc_numActNodes, ofstream& BandBDataTime, ofstream& BandBDataProbability) {
		int ActNodes = 1;
		if (id == 2) { ActNodes = 2; }
		node *best = NULL, *temp1 = NULL, *temp2 = NULL, *Dominating = NULL;
		int numModels = 0, nonDominatedCount, ii, k1, k2, kk, dominatesign=0, dummyCount, indexCount, i,j, opt_status,k;
		bool dominated, sign1;
		float dummy, minRHS, maxRHS, dummy2;
		IloNum  dummySum;
		IloInt childCount;
		IloNumArray	dummyObjectives(env, numObj), dummyProbSumNumerator(env, numObj), dummyProbSumDenominator(env, numObj);
		IloNumArray	XValue(env, var.getSize()), lower(env, numObj),upper(env, numObj);

		best = current;
		bool result = 0;
		do {  //main loop after root node
			if (current) {
				dominated = false;

				try {
					printf_s("\n*********************\n Code %d\n",current->code);
					for (ii = 0; ii < numObj;ii++) {
						if (debug_mode) {
							printf_s("Lower bound of current's %dth objective %f\n", ii, current->lb[ii]);
							printf_s("Upper bound of current's %dth objective %f\n", ii, current->ub[ii]);
						}
						lower[ii] = current->lb[ii];
						upper[ii] = current->ub[ii];
					}
					if (debug_mode) {
						if (current->code == 13) {
							cout << " ";
						}
					}
					start_child = clock();
					temp1 = current;
					opt_status = Solve(id, numObj, temp1, body, filename);//end of inner solve
					numModels += 1;
					endProblem = clock();
					LocalTimeData << std::to_string(id) << "," << numModels << "," << endProblem - start_child << endl;
					node *temp3 = NULL;
					if (opt_status != 1 && opt_status != 101 && opt_status != 102) {  //optimal is denoted by 1 here
						for (j = 0; j<numObj; j++) {
							current->objValue[j] = 0;
						}
						temp2 = current->next;
						sign1 = true;
						do {
							if (temp2) { //if next exists
								temp2->previous = current->previous;
								if (current->previous->SONL && current->previous->level == current->level) {
									temp2->SONL = current->previous->SONL;
								}
								if (temp2->dominated == false) {
									best = temp2;
									sign1 = false;
								}
								else {
									temp2 = temp2->next;
								}
							}
							else { //if next does not exist, at the end of level

								   //perform dominance test using dominateListUpperLevel
								//startFilter1 = clock();
								temp1 = current;
								while (temp1) {
									if (temp1->level == current->level) {
										if (temp1->dominated == false) {
											dummyCount = temp1->weakListUpperLevel.size();
											for (kk = 0; kk<dummyCount; kk++) {
												if (temp1->weakListUpperLevel[kk]->level == current->level) {
													if (temp1->weakListUpperLevel[kk]->objValue[0] == 0 && temp1->weakListUpperLevel[kk]->objValue[1] == 0) {
													}
													else {
														ii = temp1->weakListUpperLevel[kk]->solRank;
														if (temp1->weakListUpperLevel[kk]->dominated == false && loc_solutionList[ii].dominated == false && kk != -100) {
															temp1->weakList.push_back(temp1->weakListUpperLevel[kk]->solRank);
														}
													}
												}
											}
											dummyCount = temp1->weakList.size();
											for (j = 0; j<dummyCount>0; j++) {
												i = temp1->weakList[j];
												if (loc_solutionList[i].dominated == false) {
													dominatesign = 0;
													for (ii = 0; ii<numObj; ii++) {
														if (temp1->objValue[ii] >= loc_solutionList[i].objValue[ii]) {
															dominatesign += 1;
														}
													}
													if (dominatesign == numObj) {
														loc_solutionList[i].dominated = true;
														//Dominating=tempDominatingNode->address;
													}

												}
											}
										}
										temp1 = temp1->previous;
									}
									else { break; }
								}
								//endFilter1 = clock();
								//totalFilterTime += endFilter1 - startFilter1;

								temp1 = current->previous;
								if (temp1->level == current->level) { // if previous and current are on the same level
									if (temp1->SONL) {  // if SONL is defined on the current level
										if (temp1->SONL->dominated == false) {//if SONL of current level exists and not dominated, to delete two previous level
											current->next = temp1->SONL;
											best = current->next;
											best->previous = current;

											//delete two levels above
											temp1 = current->parent->parent;
											if (temp1) {
												temp2 = temp1->previous;
												do {
													if (temp2 && temp1->level == temp2->level) {
														temp1 = temp2;
														temp2 = temp1->previous;
													}
													else {
														temp2 = temp1;
														break;
													}
												} while (1);
											}
											do {
												if (temp2 && temp2->level == ((current->level) - 3)) {
													temp3 = temp2;
													temp1 = temp2->next;
													delete temp3;
													temp2 = temp1;
												}
												else { break; }
											} while (1);
											sign1 = false;
										}
										else { //if SONL of current level exists and but dominated
											best = temp1->SONL->next;
											do {
												if (best) {
													if (best->dominated) {
														best = best->next;
													}
													else {
														sign1 = false;
														best->previous = current;//Found the nondominated SNL in the next level, delete two previous level nodes

																				 //delete two levels above
														temp1 = current->parent->parent;
														if (temp1) {
															temp2 = temp1->previous;
															do {
																if (temp2 && temp1->level == temp2->level) {
																	temp1 = temp2;
																	temp2 = temp1->previous;
																}
																else {
																	temp2 = temp1;
																	break;
																}
															} while (1);
														}
														do {
															if (temp2 && temp2->level == ((current->level) - 3)) {
																temp3 = temp2;
																temp1 = temp2->next;
																delete temp3;
																temp2 = temp1;
															}
															else { break; }
														} while (1);


														break;
													}
												}
												else {
													sign1 = false;
													//best->previous = current;
													cout << " yes, this case was also possible!" << endl;
													BandBData << "exit4" << "\t" << endl;
													break;
												}
											} while (1);
										}
									}
									else { // if SONL is defined on the current level
										cout << "TREE TRAVERSAL IS COMPLETED! " << endl;
										BandBData << "exit1" << "\t" << endl;
										break;  //stopping of the algorithm
									}
								}
								else {  // if a node's next node does not exist an it is on the next level
									cout << " yes, this case was possible!" << endl;
									BandBData << "exit3" << "\t" << endl;
									sign1 = false;
								}
							}
						} while (sign1);
						
					}
					else {	// current node generates a feasible solution	
						dominatingNode tempDominatingNode(env);
						cplex.getValues(XValue, var);
						for (int w = (numObj - 1); w >= 0; w--)
						{
							dummySum = 0;
							for (int var_no = 0; var_no < var.getSize(); var_no++) {
								dummySum += obj_coeff_vector[(w*var.getSize() + var_no)] * XValue[var_no];
								//printf_s("\n At TREE Coefficient of w %d : %f\n", w, obj_coeff_vector[(w*var.getSize() + var_no)]);
								//printf_s("At TREE value of binary %d : %f\n", w, XValue[var_no]);
							}
							//printf_s("Objective value of obj %d DC %f\n",w, dummySum);
							tempDominatingNode.objValue[w] = dummySum;
							current->objValue[w] = dummySum;
						}
						

						BandBData << current->code << "\t";
						for (ii = 0; ii < numObj; ii++) {
							BandBData << current->lb[ii] << "\t";
							BandBData << current->objValue[ii] << "\t";
							BandBData << current->ub[ii] << "\t";

						}
						BandBData << endl;
						loc_solutionList.push_back(tempDominatingNode);
						current->solRank = loc_solutionList.size() - 1;

						//cout << loc_solutionList[loc_solutionList.size() - 1].objValue[0] << " " << loc_solutionList[loc_solutionList.size() - 1].objValue[1] << " " << loc_solutionList[loc_solutionList.size() - 1].objValue[2] << endl;
						if (loc_solutionList[loc_solutionList.size() - 1].objValue[0] == 1223 && loc_solutionList[loc_solutionList.size() - 1].objValue[1] == 1235) {
							cout << " ";
						}
						if (loc_solutionList[loc_solutionList.size() - 1].objValue[0] == 1246 && loc_solutionList[loc_solutionList.size() - 1].objValue[1] == 1252) {
							cout << " ";
						}
						if (loc_solutionList[loc_solutionList.size() - 1].objValue[0] == 1264 && loc_solutionList[loc_solutionList.size() - 1].objValue[1] == 1240) {
							cout << " ";
						}
						//Update domination list with children of previous list***************************************************************************************************************************************************************
						dummyCount = 0;
						if (current->parent){dummyCount = current->parent->dominantsListUpperLevel.size(); } //find size of list
						if (dummyCount>0) {
							current->dominantsList = current->parent->dominantsList;
							for (ii = 0; ii<dummyCount; ii++) {

								temp1 = current->parent->dominantsListUpperLevel[ii];
								if (temp1) {
									if (temp1->level == ((current->level) - 1)) {
										for (j = 0; j<numObj ; j++) {

											temp2 = temp1->child[j];
											if (temp2) {// && (temp2->objValue[0]!=0
												current->dominantsListUpperLevel.push_back(new node(env));
												current->dominantsListUpperLevel.back() = temp2;
												current->dominantsList.push_back(temp2->solRank); // add the content of tempDominatingNode2 to the dominationList of curentNode
											}
										}
									}
									else if (temp1->level<((current->level) - 1)) {
										current->dominantsListUpperLevel.erase(current->dominantsListUpperLevel.begin() + ii);
									}
								}
							}
						}
						//if node itself is not DC add the DCsibling to the dominantsListUpperLevel
						if (!current->DC) {
							dummyCount = current->dominantsListUpperLevel.size();
							if (current->DCsibling) {
								//if(current->DCsibling->dominated==false && (current->DCsibling->objValue[0]!=0)){  // if this child node is not dominated
								current->dominantsListUpperLevel.push_back(new node(env));
								current->dominantsListUpperLevel.back() = current->DCsibling;
								current->dominantsList.push_back(current->DCsibling->solRank); // add the content of tempDominatingNode2 to the dominationList of curentNode											
																								//}
							}
						}
						//Update weakListUpperLevel with children of previous list***************************************************************************************************************************************************************
						dummyCount = 0;
						if (current->parent) { dummyCount = current->parent->weakListUpperLevel.size(); } //find size of list						
						
						if (dummyCount>0) {
							current->weakList = current->parent->weakList;
							for (ii = 0; ii<dummyCount; ii++) {

								temp1 = current->parent->weakListUpperLevel[ii];
								if (temp1) {
									if (temp1->level == ((current->level) - 1)) {
										for (j = 0; j<numObj ; j++) {

											temp2 = temp1->child[j];
											if (temp2 && (temp2->objValue[0] != 0)) {
												current->weakListUpperLevel.push_back(new node(env));
												current->weakListUpperLevel.back() = temp2;
											}
										}
									}
									else if (temp1->level<((current->level) - 1)) {
										current->weakListUpperLevel.erase(current->weakListUpperLevel.begin() + ii);
									}
								}
							}
						}
						//if node itself is DC add the non DCsibling to the weakListUpperLevel
						if (current->DC) {
							dummyCount = current->weakListUpperLevel.size();
							temp2 = current;
							for (j = 0; j<numObj - 1; j++) {
								temp1 = temp2->next;
								if (temp1) {
									if (temp1->DCsibling) {
										if (temp1->DCsibling->code == current->code) {
											if (temp1->dominated == false && (temp1->objValue[0] != 0)) {  // if this child node is not dominated
												current->weakListUpperLevel.push_back(new node(env));
												current->weakListUpperLevel.back() = temp1;
											}
										}
									}
									temp2 = temp1;
								}
							}
						}


						//domination check with dominatingList
						dummyCount = 0;
						dummyCount = current->dominantsList.size();
						for (j = 0; j<dummyCount; j++) {
							//tempDominatingNode=new dominatingNode(envMBC);
							i = current->dominantsList[j];
							if (dummyCount>i && i>0) {
								if (loc_solutionList[i].dominated == false) {
									dominatesign = 0;
									for (ii = 0; ii<numObj; ii++) {
										if (loc_solutionList[i].objValue[ii] >= current->objValue[ii]) {
											dominatesign += 1;
										}
									}
									if (dominatesign == numObj) {
										dominated = true;
										current->dominated = true;
										loc_solutionList[current->solRank].dominated = true;
										break;
									}
								}
							}
							//delete tempDominatingNode2;									
						}
						//PROBABILITY calculation*******************************************************************************************************************************
						if (numModels % 200 == 0) {
							nonDominatedCount = 0;
							dummyCount = loc_solutionList.size();
							dummySum = 0;
							for (k = 0; k<dummyCount; k++) {
								if (loc_solutionList[k].dominated == false) {
									nonDominatedCount += 1;
									for (ii = 0; ii<numObj; ii++) {
										dummyObjectives[ii] = loc_solutionList[k].objValue[ii];
									}

									//Probability based on the area of previous child nodes on the same level
									temp1 = current->previous;
									for (ii = 0; ii<numObj; ii++) {
										dummyProbSumNumerator[ii] = 0;
										dummyProbSumDenominator[ii] = 0;
									}
									while (temp1) {
										if (temp1->level == current->level) {
											for (ii = 0; ii<numObj; ii++) {
												temp2 = temp1->child[ii];
												if (temp2) {
													if (temp2->ub[0] >= dummyObjectives[0] && temp2->dominated == false) {
														indexCount = 0;
														for (j = 1; j<numObj; j++) {
															if (temp2->ub[j] >= dummyObjectives[j]) {
																indexCount += 1;
															}
														}
														if (indexCount == numObj - 1) {
															for (j = 0; j<numObj; j++) {
																//cout<<dummyObjectives[j]<<" "<<temp2->ub[j]<<endl;
																if (temp2->lb[j]<dummyObjectives[j]) {
																	dummyProbSumNumerator[j] += (temp2->ub[j] - dummyObjectives[j]);
																}
																else {
																	dummyProbSumNumerator[j] += (temp2->ub[j] - temp2->lb[j]);
																}
																dummyProbSumDenominator[j] += (temp2->ub[j] - temp2->lb[j]);
																//cout<<dummyProbSumNumerator[j]<<" "<<dummyProbSumDenominator[j]<<endl;
															}
														}

													}
												}
											}
											temp1 = temp1->previous;
											if (!temp1) {
												break;
											}
										}
										else { break; }
									};
									//Probability based on the area of next nodes on the same level***********************************************************************
									temp1 = current->next;
									while (temp1) {
										if (temp1->level == current->level) {
											if (temp1->ub[0] >= dummyObjectives[0] && temp1->dominated == false) {
												indexCount = 0;
												for (j = 1; j<numObj; j++) {
													if (temp1->ub[j] >= dummyObjectives[j]) {
														indexCount += 1;
													}
												}
												if (indexCount == numObj - 1) {
													for (j = 0; j<numObj; j++) {
														//cout<<dummyObjectives[j]<<" "<<temp1->ub[j]<<endl;
														if (temp1->lb[j]<dummyObjectives[j]) {
															dummyProbSumNumerator[j] += (temp1->ub[j] - dummyObjectives[j]);
														}
														else {
															dummyProbSumNumerator[j] += (temp1->ub[j] - temp1->lb[j]);
														}
														dummyProbSumDenominator[j] += (temp1->ub[j] - temp1->lb[j]);
														//cout<<dummyProbSumNumerator[j]<<" "<<dummyProbSumDenominator[j]<<endl;
													}
												}
											}
											temp1 = temp1->next;
											if (!temp1) {
												break;
											}
										}
										else { break; }
									};
									loc_solutionList[k].Probability = 1;
									for (ii = 0; ii<numObj; ii++) {
										//cout<<dummyProbSumNumerator[ii]<<" "<<dummyProbSumDenominator[ii]<<endl;
										dummy = dummyProbSumNumerator[ii] / dummyProbSumDenominator[ii];
										if (dummy>0) {
											loc_solutionList[k].Probability = loc_solutionList[k].Probability*dummy;
										}
										else {
											loc_solutionList[k].Probability = 0;
											break;
										}
									}
									if (k<200) {
										BandBDataProbability << k << "\t" << loc_solutionList[k].Probability << endl;
									}
									dummySum += loc_solutionList[k].Probability;
								}
								else {
									loc_solutionList[k].Probability = 1;
								}
							}
							BandBDataProbability << "\t" << numModels << "\t" << dummySum / nonDominatedCount << "\t" << nonDominatedCount << endl;
						}
						//end of probability calculation****************************************************
						//create the children and replace the info inside (NORMAL)
						//create the children and replace the info inside (NORMAL)
						node *pchild[3] = { NULL };
						for (j = 0; j<numObj; j++) {
							pchild[j] = new node(env);
							current->child[j] = pchild[j];
							for (ii = 0; ii<numObj; ii++) {
								//printf_s("Lower bound of child %d third objective %f\n", j, lower[ii]);
								//printf_s("Upper bound of child %d third objective %f\n", j, upper[ii]);
								pchild[j]->lb[ii] = lower[ii];             // store lower/upper bounds of objectives. It will be inherited to the child nodes
								if (ii == 0) { pchild[j]->ub[ii] = current->objValue[ii]; }
								else {
									pchild[j]->ub[ii] = upper[ii];
								}
							}
						}

						childCount = 0;
						for (k1 = 1; k1<3; k1++) {
							for (k2 = 1; k2<3; k2++) {
								if (k1 == 1 && k2 == 1) {
									pchild[childCount]->lb[1] = current->objValue[1];
									pchild[childCount]->lb[2] = current->objValue[2];
									loc_numActNodes += 1;
									pchild[childCount]->code = loc_numActNodes;
									pchild[childCount]->DC = true;
									childCount += 1;
								}
								else if (k1 == 1 && k2 == 2) {
									pchild[childCount]->lb[1] = current->objValue[1];
									pchild[childCount]->ub[2] = current->objValue[2];
									loc_numActNodes += 1;
									pchild[childCount]->code = loc_numActNodes;
									pchild[childCount]->DCsibling = pchild[0];
									childCount += 1;
								}
								else if (k1 == 2 && k2 == 1) {
									pchild[childCount]->ub[1] = current->objValue[1];
									pchild[childCount]->lb[2] = current->objValue[2];
									loc_numActNodes += 1;
									pchild[childCount]->code = loc_numActNodes;
									pchild[childCount]->DCsibling = pchild[0];
									childCount += 1;
								}
								else { continue; }
							}
						}

						for (k1 = 0; k1<numObj; k1++) {
							temp1 = current->child[k1];
							for (j = 0; j<numObj; j++) {
								/*
								printf_s("Lower bound of child %d third objective %f\n", k1, temp1->lb[2]);
								printf_s("Upper bound of child %d third objective %f\n", k1, temp1->ub[2]);
								printf_s("Lower bound of child %d second objective %f\n", k1, temp1->lb[1]);
								printf_s("Upper bound of child %d second objective %f\n", k1, temp1->ub[1]);
								printf_s("Lower bound of child %d first objective %f\n", k1, temp1->lb[0]);
								printf_s("Upper bound of child %d first objective %f\n", k1, temp1->ub[0]);*/

								if (j == 0) {
									minRHS = current->objValue[j];
									maxRHS = 0;
								}
								else {
									minRHS = IloMin(temp1->ub[j], current->ub[j]);
									maxRHS = IloMax(temp1->lb[j], current->lb[j]);

								}
								temp1->lb[j] = maxRHS;
								temp1->ub[j] = minRHS;
								if (temp1->ub[j] == temp1->lb[j]) {
									temp1->dominated = true;
									continue;
								}
							}
						}
						
						//FATHOMING due to main_dc_solutions ******************************************
						//if (current->dominated == true) {// if current node is dominated
						dummyCount = main_dc_solutions.size();
						for (i = 0; i<childCount; i++) {
							temp2 = current->child[i];
							if (temp2 && !(temp2->dominated)) {
								for (j = 0; j<dummyCount; j++) {
									dominatesign = 0;
									for (ii = 0; ii<numObj; ii++) {
										if (main_dc_solutions[j].objValue[ii] >= temp2->ub[ii]) {
											dominatesign += 1;
										}
									}
									if (dominatesign == numObj) {
										//printf_s("A node fathomed!!!");
										fathommed += 1;
										dominated = true;
										temp2->dominated = true;
										break;
									}
								}
							}

						}
						//FATHOMING******************************************
						dummyCount = current->dominantsList.size();
						for (i = 0; i<childCount; i++) {
							temp2 = current->child[i];
							if (temp2 && !(temp2->dominated)) {
								for (j = 0; j<dummyCount; j++) {
									kk = current->dominantsList[j];
									if (kk>0 && kk<dummyCount) {
										if (loc_solutionList[kk].dominated == false) {
											dominatesign = 0;
											for (ii = 0; ii<numObj; ii++) {
												if (loc_solutionList[kk].objValue[ii] >= temp2->ub[ii]) {
													dominatesign += 1;
												}
											}
											if (dominatesign == numObj) {
												//printf_s("A node fathomed!!!");
												fathommed += 1;
												dominated = true;
												temp2->dominated = true;
												break;
											}
										}
									}
								}
							}
								//}
							}// end of specialties for childrenof current, if node itself is dominated


						for (j = 0; j<numObj; j++) {
							current->child[j]->parent = current;
							current->child[j]->level = (current->level) + 1;
						}
						//connect siblings
						current->child[0]->next = current->child[1];
						current->child[1]->next = current->child[2];
						numPareto += 1;

						//enter SONL information to the solved node; connect last child of previous node to the current node's first child if the levels are same

						temp1 = current->previous;
						if (temp1->level == current->level) {
							temp2 = temp1->child[numObj - 1];
							temp2->next = current->child[0];
							current->SONL = temp1->SONL;
						}
						
						else if (current->child[0]) {
							current->SONL = current->child[0];
						}
						
						//else{
						//	current->SONL=NULL;
						//}
						loc_solutionList[loc_solutionList.size() - 1].level = current->level;
						loc_solutionList[loc_solutionList.size() - 1].dominated = current->dominated;

						do {
							temp1 = current->next;
							if ((!temp1) || temp1->dominated == false) {// current next does not exist or temp1 exists but nondominated
								break;
							}
							else if (temp1->dominated == true) {
								current->next = temp1->next;
							}
						} while (1);

						if (!current->next) {
							current->next = current->SONL;

							//perform dominance test using dominateListUpperLevel
							//startFilter2 = clock();
							temp1 = current;
							while (temp1) {
								if (temp1->code == 7) {
									cout << " ";
								}
								if (temp1->level == current->level) {
									if (temp1->dominated == false) {
										dummyCount = temp1->weakListUpperLevel.size();
										for (kk = 0; kk < dummyCount; kk++) {
											if (temp1->weakListUpperLevel[kk]->level == current->level) {
												if (temp1->weakListUpperLevel[kk]->objValue[0] == 0 && temp1->weakListUpperLevel[kk]->objValue[1] == 0) {
												}
												else {
													ii = temp1->weakListUpperLevel[kk]->solRank;
													if (temp1->weakListUpperLevel[kk]->dominated == false && loc_solutionList[ii].dominated == false && kk != -100) {
														temp1->weakList.push_back(temp1->weakListUpperLevel[kk]->solRank);
													}
												}
											}
										}
										dummyCount = temp1->weakList.size();
										for (j = 0; j < dummyCount>0; j++) {
											i = temp1->weakList[j];
											if (loc_solutionList[i].dominated == false) {
												dominatesign = 0;
												for (ii = 0; ii < numObj; ii++) {
													if (temp1->objValue[ii] >= loc_solutionList[i].objValue[ii]) {
														dominatesign += 1;
													}
												}
												if (dominatesign == numObj) {
													loc_solutionList[i].dominated = true;
													//Dominating=tempDominatingNode->address;
												}

											}
										}
									}
									temp1 = temp1->previous;
								}
								else { break; }
							}
							//endFilter2 = clock();
							//totalFilterTime += endFilter2 - startFilter2;


							//delete two levels above
						if(current->level >= 4){
							temp1 = current->parent->parent;
							if (temp1) {
								temp2 = temp1->previous;
								do {
									if (temp2 && temp1->level == temp2->level) {
										temp1 = temp2;
										temp2 = temp1->previous;
									}
									else {
										temp2 = temp1;
										break;
									}
								} while (1);
							}
							do {
								if (temp2 && temp2->level == ((current->level) - 2)) {
									temp3 = temp2;
									temp1 = temp2->next;
									delete temp3;
									temp2 = temp1;
								}
								else { break; }
							} while (1);
							}
						}
						best = current->next;
						best->previous = current;

					}//end of else which claims current node generated a solution
					
					//BandBDataTime << current->code << "\t" << endProblem - start_child << "\t" << end_tree_traverse - endProblem << "\t" << endl;
					current = best;
					best = NULL;

				}
				catch (IloException& e) {
					cerr << "Concert exception caught: " << e << endl;
				}
				catch (...) {
					cerr << "Unknown exception caught" << endl;
				}
			}//end of if current node exists
			else {
				cout << "no current node exists!" << endl;
				BandBData << "exit2" << "\t" << endl;
				break;
			}
		} while (1);
		result = 1;
		//BandBData.close();
		//BandBDataTime.close();
		cplex.end();
		env.end();
		envMOBB.end();
		return result;
	};

public:
	int Id;
private:
	IloEnv env;
	IloEnv envMOBB;
	IloModel model;
	IloModel model2;
	IloCplex cplex;
	IloObjective obj;
	IloNumVarArray var;
	IloRangeArray rng;
	IloNumArray2 objCoef;
	IloExpr	lin_expr, lin_expr2;
	std::vector<dominatingNode>  loc_solutionList;
	std::vector<dominatingNode>  main_dc_solutions;
	node *root_loc;
	node *root_sibling_loc;
	vector<IloNum> obj_coeff_vector;
	int nrows,numModels;
	IloConstraintArray	linCons2;
	IloConstraintArray	linCons;
	std::vector<IloExpr>  list_objs;
	clock_t endProblem;
	clock_t startFilter;
	clock_t start_child;
	clock_t end_solve;
	clock_t end_tree_traverse;
	ofstream BandBData;
	ofstream LocalTimeData;
	//ofstream BandBDataProbability;
	int numActNodes;
	IloInt numPareto;
	bool debug_mode;
	int fathommed;
};
