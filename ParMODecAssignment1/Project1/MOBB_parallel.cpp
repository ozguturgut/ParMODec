#include <iostream>
#include <string>
#include <vector>
#include "ilcplex/ilocplex.h"
#include <omp.h> 
#include <tuple>
#include "dominatingNode.h"
#include "node.h"
#include "root_solve.h"
#include "child_solve.h"
#include "misc.h"

using namespace std;
#define NMax (4)


typedef IloArray<IloNumVarArray> IloNumVarArray2;
typedef IloArray<IloNumExprArray> IloExprArray2;
//typedef IloArray<IloExpr> IloExprArray;

vector<root_solve*> solvers;
vector<child_solve*> child_solvers;
vector<node*> parallel_roots;

vector<dominatingNode>  solutionList_pieces[NMax];
std::vector<dominatingNode>  solutionList;
std::vector<dominatingNode>  filtered_solutionList;

ofstream P1FS_Pareto;
ofstream P1FS_Time;
ofstream P1FS_Probability;

int main(void)
{
	bool parallel = 0;
	int numVar;
	long long numObj_l = 3;
	long long numVar_l = 5;
	long long instance_l = 1;
	bool debug_mode=true;

	P1FS_Time.open("Timelog.txt", ios::out | ios::app);
	log_it_date(P1FS_Time);

	std::string body("C:\\Users\\Ozgu\\Documents\\MOBB Collection\\MOBB_paper\\Inputs Sept 2016\\Ozgu Assignment");
	//std::string body("C:\\Users\\ozgutr\\Documents\\MOBB\\MOBB_paper\\Inputs Sept 2016");
	std::string prob_type("Assignment");
	for (int numObj = 3; numObj < 4; numObj++) {
		numObj_l = numObj;
		for (int myj = 1; myj < 4; myj++) {
			if (myj == 1) { numVar = 5; }
			if (myj == 2) { numVar = 7; }
			if (myj == 3) { numVar = 10; }
			numVar_l = numVar;
			numVar = numVar*numVar;
			for (int instance = 1; instance < 6; instance++) {
				instance_l = instance;
				//numObj_l = 3;
				//numVar_l = 5;
				//instance_l = 1;
				int num_children;
				char stop_flag;
				bool solvedAll = true;
				//int num_threads = 1;
				int num_obj = numObj_l;
				float** anchors;
				int num_models_solved[4];
				int total_models_solved = 0;
				clock_t startproblem, endfilter, endProblem;
				startproblem = clock();
				std::string filename("\\AP_p-" + std::to_string(numObj_l) + "_n-" + std::to_string(numVar_l) + "_ins-" + std::to_string(instance_l));
				//std::string filename("\\LokmanKoksalan");
				//std::string filename("\\Sayin_Knapsack_sample20");
				std::string PARETO_filename("ParMoDec_AP-" + std::to_string(numObj_l) + "_n-" + std::to_string(numVar_l) + "_ins-" + std::to_string(instance_l) + ".txt");

				P1FS_Pareto.open(PARETO_filename, ios::out | ios::app);
				P1FS_Probability.open("Stat_log.txt", ios::out | ios::app);

				std::string lpfilelocation(body + filename + ".lp");
				printf_s("File %s\n", lpfilelocation.c_str());
				root_solve* solve_at_a_time;
				//###################################################################
				//anchor
				
				solve_at_a_time = new root_solve(10, lpfilelocation, num_obj); // id 10 is used for anchor, i.e. single thread solve
				std::tie(anchors, parallel_roots, solutionList, total_models_solved) = solve_at_a_time->find_anchors(10, num_obj, filename, body, numVar,instance, P1FS_Time, P1FS_Probability);
				P1FS_Time << prob_type << "," << numObj << "," << numVar << "," << instance << endl;
				P1FS_Time << "Root solve time:" << clock()- startproblem<< endl;
				//root_node_solve is called in find_anchors 
				printf_s("File %f\n", anchors[0][0]);
				if (debug_mode) {
					printf_s("# of Parallel Roots: %d\n", parallel_roots.size());
					for (int i = 0; i < parallel_roots.size(); i++) {
						if (parallel_roots[i]){
							printf_s("\n");
							printf_s("Lower bound of root %d first objective %f\n", i, parallel_roots[i]->lb[1]);
							printf_s("Lower bound of root %d second objective %f\n", i, parallel_roots[i]->lb[2]);
							printf_s("Lower bound of root %d main objective %f\n", i, parallel_roots[i]->lb[0]);
							printf_s("Upper bound of root %d first objective %f\n", i, parallel_roots[i]->ub[1]);
							printf_s("Upper bound of root %d second objective %f\n", i, parallel_roots[i]->ub[2]);
							printf_s("Upper bound of root %d main objective %f\n", i, parallel_roots[i]->ub[0]);
							printf_s("  ****************  ");
							/*if (parallel_roots[1] && (i == 1)) {
								printf_s("\n");
								printf_s("Lower bound of root %d 's sibling first objective %f\n", i, parallel_roots[i]->next->lb[1]);
								printf_s("Lower bound of root %d  's sibling second objective %f\n", i, parallel_roots[i]->next->lb[2]);
								printf_s("Lower bound of root %d  's sibling main objective %f\n", i, parallel_roots[i]->next->lb[0]);
								printf_s("Upper bound of root %d  's sibling first objective %f\n", i, parallel_roots[i]->next->ub[1]);
								printf_s("Upper bound of root %d  's sibling second objective %f\n", i, parallel_roots[i]->next->ub[2]);
								printf_s("Upper bound of root %d  's sibling main objective %f\n", i, parallel_roots[i]->next->ub[0]);
							}*/
						}
					}
				}
				//###################################################################
				//child_node_solve-PARALLEL
				for (int i = 0; i < NMax; i++) {
					child_solvers.push_back(new child_solve(i));
				}

#pragma omp parallel num_threads(NMax)
				{
					bool solvedHere = true;
					bool if_solved = true;
					int th;
#pragma omp for
					for (int i = 0; i < NMax; i++) {
						th = omp_get_thread_num();
						node *temp = NULL, *temp2 = NULL;
						IloEnv env;
						node root(env);
						node root_sibling(env);
						if (parallel_roots[i]){
							temp = parallel_roots[i];
							root = temp->copy(env);
							//printf_s("Lower bound of root %d's second objective %f\n", i, root_sibling.lb[1]);
							//printf_s("Lower bound of root %d's third objective %f\n", i, root_sibling.lb[2]);
							root_sibling = temp->copy(env); // we send the same node as if sibling to the nodes other than 2 
						
							std::tie(if_solved, solutionList_pieces[i], num_models_solved[i]) = child_solvers[i]->solve_main(i, body, filename, &root, &root_sibling, num_obj, solutionList, P1FS_Time, P1FS_Probability,numVar,instance);
							total_models_solved = total_models_solved + num_models_solved[i];
							solvedHere = solvedHere && if_solved;
						}

					}
//#pragma omp critical
					//{
					//	solvedAll = solvedAll && solvedHere;

					//}
#pragma omp barrier 
				}//end of pragma parallel
				endProblem = clock();
				cout << (solvedAll ? "Success" : "Fail") << endl;
				if (solvedAll) {
					for (int i = 0; i < NMax; i++) {
						//solutionList.insert(solutionList.end(), solutionList_pieces[i].begin(), solutionList_pieces[i].end());
						solutionList = merge_parallels(solutionList_pieces[i], solutionList, num_obj);
						
					}
					filtered_solutionList = pareto_filter(solutionList, num_obj);
				}
				//cin.get();
				endfilter = clock();
				cout << "Total solve time: " << endProblem - startproblem << endl;
				cout << "Total number of models: " << total_models_solved << endl;
				P1FS_Time  << "# of Models Solved" << "," << "# of Pareto" << "," << "Total Run Time" << "," << "Filter Time"<< endl;
				P1FS_Time <<total_models_solved << "," << filtered_solutionList.size() << "," << endProblem - startproblem << "," << endfilter- endProblem <<endl;
				log_it_vector(P1FS_Pareto, filtered_solutionList, num_obj);
				//system("pause");

				for (auto it = child_solvers.begin(); it != child_solvers.end(); it++) {
					delete * it;
					it = child_solvers.erase(it);
				}
				child_solvers.clear();
				for (auto it = solvers.begin(); it != solvers.end(); it++) {
					delete * it;
					it = solvers.erase(it);
				}
				solvers.clear();
				for (auto it = parallel_roots.begin(); it != parallel_roots.end(); it++) {
					delete * it;
					it = parallel_roots.erase(it);
				}
				parallel_roots.clear();
				for (int i = 0; i < NMax; i++) {
					solutionList_pieces[i].clear();
				}
				solutionList.clear();
				filtered_solutionList.clear();
				
			}
		}
	}
	return 0;
}