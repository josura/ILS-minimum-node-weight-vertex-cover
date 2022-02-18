#include <algorithm>
#include <cstddef>
#include <iostream>
#include <string>
#include <fstream>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>
#include <chrono>

#include "IteratedLocalSearch.h"
#include "LocalSearch.h"
#include "WeightedVertexGraph.h"
#include "utilities.h"
#include "heuristics.h"


using namespace std;

int main(int argc, char** argv)
{
    if (argc <2) {
        printUsage("experiments");
        return 0;
    }
    double * nodesWeight;
    WeightedVertexGraph* graph;
    ifstream input (argv[1]);
    if (input.fail()) {
        cerr <<"[MAIN]: "<< "input failed";
        printUsage("experiments");
        return 0;
    }
    if (input.is_open())
    {
        string numNodesStr;
        input >> numNodesStr;
        stringstream numstream(numNodesStr);
        int numNodes;
        numstream >> numNodes;
        nodesWeight = new double[numNodes];
        for (int i =0; i < numNodes; i++) {
            string weightStr;
            input >> weightStr;
            stringstream weightstream(weightStr);
            weightstream >> nodesWeight[i];
        }


        graph = new WeightedVertexGraph(numNodes,nodesWeight);

        for (int i =0; i < numNodes; i++) {
            for (int j =0; j < numNodes; j++) {    
                string edgeStr;
                input >> edgeStr;
                if (edgeStr == "1") {                
                    graph->addEdge(i, j);
                }
            }    
        }
        input.close();



        ofstream results;
        results.open ("resultsILS.csv",ios_base::app);  //appending to file

        //LOCAL SEARCH HILL CLIMBING EXPERIMENTS
        // for (double samplingFactor = 1 ; samplingFactor <  20; samplingFactor+=1) {
        //     LocalSearch* locSear = new LocalSearch(graph,samplingFactor);
        //     double cost;

        //     auto started = std::chrono::high_resolution_clock::now();
        //     NodeBitArray solution = locSear->startResolveFinal(cost);
        //     auto done = std::chrono::high_resolution_clock::now();
            
        //     //NodeList* solutionList = nodeBitArrayToList(solution,graph->getNumNodes());
        //     std::cout << "obj eval after: "<<locSear->getCurrentNumberOfObjectiveFunctionEval();

        //     //if(vertexCoverValidityEdgescheckBitArray(graph, solution))
        //     //    cout << "final solution: " << *solutionList << endl;
        //     results << "localSearch," << argv[1] << "," << std::chrono::duration_cast<std::chrono::microseconds>(done-started).count() << "," <<
        //                             locSear->getSolutionWeight() << "," <<
        //                             locSear->getCurrentNumberOfObjectiveFunctionEval() << "," <<
        //                             nodeBitArrayToString(solution,graph->getNumNodes()) <<
        //                             endl;

        //     delete locSear;
        
        // }

        //ILS EXPERIMENTS
        for (double propagationProportion = 5 ; propagationProportion <  20; propagationProportion+=1) {
            IteratedLocalSearch* ILS = new IteratedLocalSearch(graph,propagationProportion);
            double cost;

            auto started = std::chrono::high_resolution_clock::now();
            NodeBitArray solution = ILS->startResolve();
            auto done = std::chrono::high_resolution_clock::now();
            
            //NodeList* solutionList = nodeBitArrayToList(solution,graph->getNumNodes());
            std::cout << "obj eval after: "<<ILS->getObjectiveEvalForOptimum();

            //if(vertexCoverValidityEdgescheckBitArray(graph, solution))
            //    cout << "final solution: " << *solutionList << endl;
            results << "localSearch," << argv[1] << "," << std::chrono::duration_cast<std::chrono::microseconds>(done-started).count() << "," <<
                                    ILS->getSolutionWeight() << "," <<
                                    ILS->getObjectiveEvalForOptimum() << "," <<
                                    nodeBitArrayToString(solution,graph->getNumNodes()) <<
                                    endl;

            delete ILS;
        
        }
        results.close();

        //ILS testing



        // IteratedLocalSearch* ILS = new IteratedLocalSearch(graph);

        // auto started = std::chrono::high_resolution_clock::now();
        // NodeBitArray solution = ILS->startResolve();
        // auto done = std::chrono::high_resolution_clock::now();
        
        // //NodeList* solutionList = nodeBitArrayToList(solution,graph->getNumNodes());

        // //if(vertexCoverValidityEdgescheckBitArray(graph, solution))
        // //    cout << "final solution: " << *solutionList << endl;

        // cout <<"[MAIN]: "<< "solution weight:" << ILS->getSolutionWeight() << endl;

        // cout <<"[MAIN]: "<< "execution time: " << std::chrono::duration_cast<std::chrono::microseconds>(done-started).count() << " microseconds" << endl;

    }
    cout << endl;

    delete graph;
    delete [] nodesWeight;
}