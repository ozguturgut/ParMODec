#pragma once
#include <iostream>
#include <string>
#include <vector>
#include "ilcplex/ilocplex.h"
#include <omp.h> 
#include <tuple>

#include "dominatingNode.h"
#include "node.h"
#include "misc.h"
#include "child_solve.h"
using namespace std;
#define NMax (4)
//typedef IloArray<IloExpr> IloExprArray;
class root_solve
{
public:
	class root_solve
	(int id, const string &fileLocation, int numobj)
		: Id(id), model(env), cplex(env), var(env), rng(env), lin_expr(env),
		obj(env),objArray(env,numobj)
	{
		cplex.importModel(model, fileLocation.c_str(), obj, var, rng);
		debug_mode = true;
		if(debug_mode){cout << model;}
	}

	//################################################################################# FIND ANCHORS ############################################################################
	std::tuple<float**, vector<node*>, std::vector<dominatingNode>, int > find_anchors(int th, int num_obj, const string &fileName,
		const string &body, int numVar,int instance,
		ofstream& BandBDataTime, ofstream& BandBDataProbability) {
		//initialize the anchors array
		float** anchor = new float*[2];
		for (int h = 0; h < 2; h++)
		{
			anchor[h] = new float[num_obj];
			for (int w = 0; w < num_obj; w++)
			{
				anchor[h][w] = 0;
			}
		}
		
		//read LP file and convert i to the required model
		try
		{
			cplex.setOut(env.getNullStream());
			IloBool start_tree;
			long long myth = 0;
			envMOBB = model.getEnv();
			IloObjective objMOBB(envMOBB);
			IloConstraintArray	linCons2(envMOBB);
			IloConstraintArray	linCons(envMOBB);
			vector<IloIntVar> varVector;
			vector<IloNum> obj_coeff_vector;
			IloCplex cplexMOBB(envMOBB);
			IloModel model2(envMOBB);
			int nrows;
			myth = th;
			std::string lpfilelocation;
			lpfilelocation = body + fileName + "_Anchor_" + std::to_string(myth) + ".lp";
			cplex.extract(model);
			nrows = cplex.getNrows();
			if (debug_mode) {
				cout << "Number of rows:" << nrows;
			}
			IloExtractable extr;
			IloExpr expr;
			IloExprArray	expr_obj(env, num_obj);
			IloExprArray  objList(envMOBB, num_obj);
			bool collected_variables=false;
			bool bresult;
			float optimal_value;
			//read LP file and convert i to the required model
			//objList = new IloExpr[num_obj];
			int num_constraints, num_variables, j;
			j = 0;
			num_constraints = 0;
			//num_variables = 0;
			for (IloModel::Iterator it(model); it.ok(); ++it) {
				extr = it.operator *();
				//feasible region
				if (extr.isConstraint() == true) {
					if (num_constraints >= (nrows - num_obj)) {//set of objective functions
						IloRangeI* impl = dynamic_cast<IloRangeI *>((*it).asConstraint().getImpl());
						if (impl) {
							IloRange rng(impl);
							IloExpr expr = rng.getExpr();
							//lin_expr.clear();
							for (IloExpr::LinearIterator it2 = expr.getLinearIterator(); it2.ok(); ++it2) {
								//lin_expr += it2.getCoef() * it2.getVar();
								varVector.push_back(it2.getVar()); 
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
			model2.add(linCons2);
			cplexMOBB.setOut(envMOBB.getNullStream());
			cplexMOBB.setWarning(envMOBB.getNullStream());
			IloIntVarArray varArray(envMOBB, var.getSize());
			for (int var_no = 0; var_no < var.getSize(); var_no++) {
				varArray[var_no] = varVector[var_no];
			}
			if (debug_mode) {
				printf_s("Size:%d", var.getSize());
			}
			IloNumArray	XValue(envMOBB, var.getSize());

			for (int w = (num_obj - 1); w >= 0; w--)
			{	
				lin_expr.clear();
				for (int var_no = 0; var_no < var.getSize(); var_no++) {
					lin_expr += obj_coeff_vector[(var.getSize()*w + var_no)] * var[var_no];
					if (debug_mode) {
						printf_s("Coefficient %f\n", obj_coeff_vector[(var.getSize()*w + var_no)]);
					}
				}
				objList[w] = lin_expr;
				for (int h = 0; h < 2; h++)
				{
					IloObjective	obj(envMOBB);
					obj.setExpr(objList[w]);
					if (h == 0) {
						obj.setSense(IloObjective::Maximize);
					}
					else {
						obj.setSense(IloObjective::Minimize);
					}
					model2.add(obj);
					cplexMOBB.extract(model2);
					//cplexMOBB.exportModel(lpfilelocation.c_str());
					optimal_value = 0;
					bresult = cplexMOBB.solve();
					if (debug_mode) {
						printf_s("Variable Vector Size %d\n", varVector.size());
						printf_s("original Variable Vector Size %d\n", var.getSize());
					}
					if (bresult) {
						cplexMOBB.getValues(XValue, varArray);
						for (int var_no = 0; var_no < var.getSize(); var_no++) {
							//printf_s("Coefficient in ANCHOR calculation %f\n", obj_coeff_vector[(var.getSize()*w + var_no)]);
							optimal_value += obj_coeff_vector[(var.getSize()*w + var_no)] * XValue[var_no];
						}
						if (debug_mode) {
							printf_s("Anchor %d of %d value %f\n", w, h, optimal_value);
						}
						anchor[h][w] = optimal_value;
					}
					obj.end();
				}
			}
			int result=root_node_solve(model2, linCons2, anchor, num_obj, fileName, body, numVar, obj_coeff_vector, var.getSize(), BandBDataTime,BandBDataProbability,instance);
			if (debug_mode) {
				printf_s("Variable Vector Size %d\n", solutionList.size());
				printf_s("Variable Vector Size %d\n", node_list.size());
				printf_s("Variable Vector Size %f\n", node_list[node_list.size() - 1]->lb[1]);
			}

			return std::make_tuple(anchor, node_list, solutionList, numModels);
			
		}
		catch (IloException &ex)
		{
			cout << "Got exception: ";
			ex.print(cout);
			cout << endl;
			return std::make_tuple(anchor, node_list, solutionList, numModels);
		}
		catch (...)
		{
			return std::make_tuple(anchor, node_list, solutionList, numModels);
		}
		return std::make_tuple(anchor, node_list, solutionList, numModels);
		
	};
	//################################################################################# SOLVE ROOT ############################################################################
	int root_node_solve(IloModel model, IloConstraintArray linCons, float **anchor,int num_obj, const string &fileName,
		const string &body, int numVar, vector<IloNum>& obj_coeff_vector, int act_varSize,
		ofstream& P1FS_Time, ofstream& P1FS_Probability,int instance) {

		cplex.setOut(env.getNullStream());
		IloBool start_tree;
		long long myth = 11;
		envMOBB = model.getEnv();
		IloObjective objMOBB(envMOBB);
		IloConstraintArray	linCons2(envMOBB);
		IloExpr	lin_expr2(envMOBB);
		IloCplex cplexMOBB(envMOBB);
		IloModel model2(envMOBB);
		IloNum epsilon = 0.0001;
		int nrows,child_count,kk1,kk2;
		cplex.extract(model);
		nrows = cplex.getNrows();
		IloExtractable extr;
		IloExpr expr;
		IloExprArray	expr_obj(env, num_obj);
		IloNumArray	XValue(envMOBB, act_varSize);
		IloNumArray		lower(envMOBB, num_obj);
		IloNumArray		upper(envMOBB, num_obj);
		model2.add(linCons);
		lin_expr.clear();
		for (int w = (num_obj - 1); w >= 0; w--)
		{	
			lin_expr2.clear();	
				//lin_expr2 = objList[w];
			if (debug_mode) {
				printf_s("\n Coefficients read for root\n ");
			}
			for (int var_no = 0; var_no < act_varSize; var_no++) {
				//printf_s("Coefficient of w %d : %f\n", w, obj_coeff_vector[(w*act_varSize + var_no)]);
				lin_expr2 += obj_coeff_vector[(w*act_varSize + var_no)] * var[var_no];
				if (w == 0) {
					lin_expr += obj_coeff_vector[(w*act_varSize + var_no)] * var[var_no];
				}
				else {
					lin_expr += epsilon*(obj_coeff_vector[(w*act_varSize + var_no)] * var[var_no]);
				}
			}
			lower[w] = anchor[1][w];
			upper[w] = anchor[0][w];
			linCons2.add(lin_expr2 <= anchor[0][w]); //objectives were maximized for second row of anchor
			linCons2.add(lin_expr2 >= anchor[1][w]); //objectives were minimized for first row of anchor
			//printf_s("Max of anchor %d %f\n", w, anchor[0][w]);
			//printf_s("Min of anchor %d %f\n", w, anchor[1][w]);
			
		}
		objMOBB.setExpr(lin_expr);  //default  main objective is the FIRST objective in LP
		objMOBB.setSense(IloObjective::Maximize);  //default type of objective is MAXIMIZATION type
		model2.add(objMOBB);
		model2.add(linCons2);
		std::string lpfilelocation;
		lpfilelocation = body + fileName + "_Anchor_" + std::to_string(myth) + ".lp";
		cplexMOBB.setOut(envMOBB.getNullStream());
		cplexMOBB.setWarning(envMOBB.getNullStream());
		cplexMOBB.extract(model2);
		if (debug_mode) {
			cplexMOBB.exportModel(lpfilelocation.c_str());}
		bool result,presolved;
		float optimal_value = 0;
		try {
			result = cplexMOBB.solve();
			//child_solve* pre_dc_solve;
			//pre_dc_solve=new child_solve(1000);
			//std::tie(if_solved, solutionList, numModels)=pre_dc_solve->solve_main(i, body, filename, &root, &root_sibling, num_obj, solutionList, P1FS_Time, P1FS_Probability,numVar,instance);
			if (result) {
				node *current = NULL, *best = NULL, *temp1 = NULL, *temp2 = NULL, *Dominating = NULL;
				dominatingNode tempDominatingNode(envMOBB);
				current = new node(envMOBB);
				cplexMOBB.getValues(XValue, var);
				for (int w = (num_obj - 1); w >= 0; w--)
				{
					optimal_value = 0;
					for (int var_no = 0; var_no < act_varSize; var_no++) {
						//printf_s("Coefficient %f\n", obj_coeff_vector[(w*act_varSize + var_no)]);
						//printf_s("Variable Value %f\n", XValue[var_no]);
						optimal_value += obj_coeff_vector[(w*act_varSize + var_no)] * XValue[var_no];
					}
					if (debug_mode) {
						printf_s("Objective value of DC %f\n", optimal_value);
					}
					current->objValue[w] = optimal_value;
					tempDominatingNode.objValue[w] = optimal_value;
				}

				current->lb = new float[num_obj];
				current->ub = new float[num_obj];
				for (int w = (num_obj - 1); w >= 0; w--) {
					current->lb[w] = lower[w]-1;             // store lower/upper bounds of objectives. It will be inherited to the child nodes
					current->ub[w] = upper[w]+1;
				}
				solutionList.push_back(tempDominatingNode);
				current->solRank = solutionList.size() - 1;
				numModels += 1;
				child_count = 1;
				int presolved_numModels = 0;
				int num_children = 2 ^ (num_obj - 1) - 1;
				node *pchild[3] = { NULL }; //REQUIRES STATIC
				for (int j = 0; j < num_obj; j++) {
					pchild[j] = new node(envMOBB);
					current->child[j] = pchild[j];
				}
				//specifics for root node
				current->level = 1;
				pchild[0]->next = pchild[1];
				pchild[1]->next = pchild[2];
				for (int j = 0; j < num_obj; j++) {
					pchild[j]->parent = current;
					pchild[j]->level = (current->level) + 1;
				}

				current->previous = NULL;
				current->SONL = pchild[0];
				current->next = pchild[0];  // first child is Dominating Child
				pchild[0]->previous = current;
				current->probability = 0;

				lpfilelocation = body + fileName + "_Anchor_" + std::to_string(child_count) + ".lp";

				for (int ii = 0; ii < num_children; ii++) {
					pchild[ii]->lb = new float[num_obj];
					pchild[ii]->ub = new float[num_obj];
					for (int w = (num_obj - 1); w >= 0; w--) {
						//printf_s("Inherited lower %f\n", pchild[ii]->lb[w]);
						pchild[ii]->lb[w] = lower[w]; // store lower/upper bounds of objectives(Anchors for root). It will be inherited to the child nodes
						pchild[ii]->ub[w] = upper[w];
					}
				}

				IloInt childCount, k1, k2, numActNodes;
				numActNodes = 1;
				childCount = 0;
				//This loop is specific for 3 Objectives
				for (k1 = 1; k1 < num_children; k1++) {
					for (k2 = 1; k2 < num_children; k2++) {
						if (k1 == 1 && k2 == 1) {
							pchild[childCount]->lb[1] = current->objValue[1];
							pchild[childCount]->lb[2] = current->objValue[2];
							numActNodes += 1;
							pchild[childCount]->code = numActNodes;
							pchild[childCount]->DC = true;
							childCount += 1;
							//for 3 obj: DC of root node *****************************************************************************************
							model2.remove(linCons2);
							linCons2.endElements();
							for (int w = (num_obj - 1); w >= 0; w--)
							{
								lin_expr2.clear();
								if (w != 0) {
									//lin_expr2 = objList[w];
									for (int var_no = 0; var_no < act_varSize; var_no++) {
										//printf_s("Coefficient %f\n", obj_coeff_vector[(w*act_varSize + var_no)]);
										lin_expr2 += obj_coeff_vector[(w*act_varSize + var_no)] * var[var_no];
									}
									linCons2.add(lin_expr2 <= pchild[childCount-1]->ub[w]); //objectives are maximized for second row of anchor
									linCons2.add(lin_expr2 >= pchild[childCount-1]->lb[w]+1); //objectives are minimized for first row of anchor
								}
							}
							model2.add(linCons2);
							std::string lpfilelocation;
							myth = 111; //DC of root node
							lpfilelocation = body + fileName + "_Anchor_" + std::to_string(myth) + ".lp";
							cplexMOBB.extract(model2);
							if (debug_mode) {
								cplexMOBB.exportModel(lpfilelocation.c_str());
							}
							result = cplexMOBB.solve();
							if (result) {
								if (debug_mode) {
									printf_s("One level proceeded for the DC child \n");
								}
								cplexMOBB.getValues(XValue, var);
								for (int w = (num_obj - 1); w >= 0; w--)
								{
									optimal_value = 0;
									for (int var_no = 0; var_no < act_varSize; var_no++) {
										//printf_s("Coefficient %f\n", obj_coeff_vector[(w*act_varSize + var_no)]);
										//printf_s("Variable Value %f\n", XValue[var_no]);
										optimal_value += obj_coeff_vector[(w*act_varSize + var_no)] * XValue[var_no];
									}
									if (debug_mode) {
										printf_s("Objective value %f\n", optimal_value);
									}
									pchild[childCount - 1]->objValue[w] = optimal_value;
									tempDominatingNode.objValue[w] = optimal_value;
								}

								solutionList.push_back(tempDominatingNode);
								pchild[childCount - 1]->solRank = solutionList.size() - 1;
								solutionList[solutionList.size() - 1].level = pchild[childCount - 1]->level;
								solutionList[solutionList.size() - 1].dominated = pchild[childCount - 1]->dominated;
								numModels += 1;
								child_count = 1;
								node *ppchild[3] = { NULL }; //REQUIRES STATIC
								node temp_dc(envMOBB);
								
								for (int j = 0; j < num_obj; j++) {
									ppchild[j] = new node(envMOBB);
									pchild[childCount - 1]->child[j] = ppchild[j];
								}
								//specifics for root node
								ppchild[0]->next = ppchild[1];
								ppchild[1]->next = ppchild[2];
								for (int j = 0; j < num_obj; j++) {
									pchild[j]->parent = pchild[childCount - 1];
									pchild[j]->level = (pchild[childCount - 1]->level) + 1;
								}

								pchild[childCount - 1]->previous = current;
								pchild[childCount - 1]->SONL = ppchild[0];
								ppchild[0]->previous = current;
								pchild[childCount - 1]->probability = 0;

								lpfilelocation = body + fileName + "_Anchor_" + std::to_string(child_count) + ".lp";

								for (int ii = 0; ii < num_children; ii++) {
									ppchild[ii]->lb = new float[num_obj];
									ppchild[ii]->ub = new float[num_obj];
									for (int w = (num_obj - 1); w >= 0; w--) {
										//printf_s("Inherited lower %f\n", pchild[ii]->lb[w]);
										ppchild[ii]->lb[w] = pchild[childCount - 1]->lb[w]; // store lower/upper bounds of objectives(parent). It will be inherited to the child nodes
										ppchild[ii]->ub[w] = pchild[childCount - 1]->ub[w];
									}
								}
								int childCount2 = 0;
								for (kk1 = 1; kk1<3; kk1++) {
									for (kk2 = 1; kk2<3; kk2++) {
										if (kk1 == 1 && kk2 == 1) {
											ppchild[childCount2]->lb[1] = pchild[childCount - 1]->objValue[1];
											ppchild[childCount2]->lb[2] = pchild[childCount - 1]->objValue[2];
											numActNodes += 1;
											ppchild[childCount2]->code = 1;
											ppchild[childCount2]->DC = true;
											
											//node_list.push_back(ppchild[childCount2 - 1]);
											temp_dc = ppchild[childCount2]->copy(env);
											child_solve* pre_dc_solve;
											pre_dc_solve=new child_solve(1000);
											std::tie(presolved, solutionList, presolved_numModels)=pre_dc_solve->solve_main(1000, body, fileName, &temp_dc, &temp_dc, num_obj, solutionList, P1FS_Time, P1FS_Probability,numVar,instance);
											numModels = presolved_numModels + numModels;
											childCount2 += 1;
										}
										else if (kk1 == 1 && kk2 == 2) {
											ppchild[childCount2]->lb[1] = pchild[childCount - 1]->objValue[1];
											ppchild[childCount2]->ub[2] = pchild[childCount - 1]->objValue[2];
											numActNodes += 1;
											ppchild[childCount2]->code = 1;
											ppchild[childCount2]->DCsibling = ppchild[0];
											childCount2 += 1;
											node_list.push_back(ppchild[childCount2 - 1]);
										}
										else if (kk1 == 2 && kk2 == 1) {
											ppchild[childCount2]->ub[1] = pchild[childCount - 1]->objValue[1];
											ppchild[childCount2]->lb[2] = pchild[childCount - 1]->objValue[2];
											numActNodes += 1;
											ppchild[childCount2]->code = 2;
											ppchild[childCount2]->DCsibling = ppchild[0];
											childCount2 += 1;
											node_list.push_back(ppchild[childCount2 - 1]);
										}
										else { continue; }
									}
								}
								if (debug_mode) {
									printf_s("Lower bound of root 1's sibling first objective %f\n", ppchild[1]->next->lb[1]);
									printf_s("Lower bound of root 1's sibling second objective %f\n", ppchild[1]->next->lb[2]);
									printf_s("Upper bound of root 1's sibling first objective %f\n", ppchild[1]->next->ub[1]);
									printf_s("Upper bound of root 1's sibling second objective %f\n", ppchild[1]->next->ub[2]);
								}


							}// end of lower level for first DC child *****************************************************************************************
							else {//if first level DC solve is infeasible
								node_list.push_back({ NULL });
								node_list.push_back({ NULL });
							}
						}
						else if (k1 == 1 && k2 == 2) {
							pchild[childCount]->lb[1] = current->objValue[1];
							pchild[childCount]->ub[2] = current->objValue[2];
							numActNodes += 1;
							pchild[childCount]->code = numActNodes;
							pchild[childCount]->DCsibling = pchild[0];
							childCount += 1;
							node_list.push_back(pchild[childCount - 1]);
						}
						else if (k1 == 2 && k2 == 1) {
							pchild[childCount]->ub[1] = current->objValue[1];
							pchild[childCount]->lb[2] = current->objValue[2];
							numActNodes += 1;
							pchild[childCount]->code = numActNodes;
							pchild[childCount]->DCsibling = pchild[0];
							childCount += 1;
							node_list.push_back(pchild[childCount - 1]);
						}
						else { continue; }
					}
				}
				solutionList[solutionList.size() - 1].level = current->level;
				solutionList[solutionList.size() - 1].dominated = current->dominated;
				if (debug_mode) {
				cplexMOBB.extract(model2);
				cplexMOBB.exportModel(lpfilelocation.c_str());}
			}
			else {
				printf_s("root node infeasible \n");
				exit(EXIT_FAILURE);
			}
		}
		catch (IloException& ex) {
			cerr << "Error IloException:" << ex << endl;
		}
		catch (...) {
			cerr << "Error" << endl;
		}
		//env.end();
		envMOBB.end();
		//return solutionList & current node as tuple
		return 1; 
	};
	
public:
	int Id;
private:
	IloEnv env;
	IloEnv envMOBB;
	IloModel model;
	IloCplex cplex;
	IloObjective obj;
	IloNumVarArray var;
	IloRangeArray rng;
	IloNumArray2 objCoef;
	IloExpr	lin_expr;
	std::vector<dominatingNode>  solutionList;
	std::vector<node*>  node_list;
	int numModels;
	IloExprArray objArray;
	bool debug_mode;
};