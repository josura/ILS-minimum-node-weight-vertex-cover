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
    double * nodesWeight;
    WeightedVertexGraph* graph;
    ifstream input (argv[1]);
    if (input.fail()) {
        cerr <<"[MAIN]: "<< "input failed"<<endl;
        printUsage("vertex_cover");
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
        
        
        /*
        //cout << *graph << endl;
        size_t numberNodes = graph->getNumNodes();
        NodeList* solution= new NodeList{1,2,4,5,6,7,8,9,12,14,15,17,18,19};
        NodeBitList* bitSolution = new NodeBitList{false,true,true,false,true,true,true,true,true,true,false,false,true,false,true,true,false,true,true,true};
        NodeList* notsolution= new NodeList{1,2,4};
        NodeBitList* notbitsolution = new NodeBitList{false,true,true,false,true,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false};
        //bool edgescovered = vertexCoverValidityEdgescheck(graph, solution);
        bool edgescoveredBit = vertexCoverValidityEdgescheckBitList(graph, bitSolution);
        double cost = costFunctionBitList(graph, bitSolution);
        cost = costFunction(graph, solution);
        bool notcoveredBit = vertexCoverValidityEdgescheckBitList(graph, notbitsolution);
        */

        //GENERATORS OF SOLUTIONS
        /*graph->makeEdgesArray();  // REALLY IMPORTANT
        //NodeBitArray solution13 = randomBooleanArray(graph->getNumNodes());
        NodeBitArray solution13 = new bool[graph->getNumNodes()];  
        fill_n(solution13, graph->getNumNodes(), false);
        NodeList* solutionList13 = nodeBitArrayToList(solution13,graph->getNumNodes());
        cout << "pre array: " << *solutionList13 << endl;
        delete solutionList13;
        solution13 = randomSmartSolutionBitArrayFromPartial(graph,solution13 );
        solutionList13 = nodeBitArrayToList(solution13,graph->getNumNodes());
        if (vertexCoverValidityEdgescheckBitArray(graph, solution13)) cout << "post array: " << *solutionList13 << endl;
        */


        //LOCAL SEARCH HILL CLIMBING TESTING
        // LocalSearch* localSearch = new LocalSearch(graph);
        // cout << "greedy initial solution weight:" << localSearch->getSolutionWeight() << endl;

        // double cost;

        // auto started = std::chrono::high_resolution_clock::now();
        // NodeBitArray solution = localSearch->startResolveWithLimitAndSamplingOptimized(cost);
        // auto done = std::chrono::high_resolution_clock::now();
        
        // NodeList* solutionList = nodeBitArrayToList(solution,graph->getNumNodes());

        // if(vertexCoverValidityEdgescheckBitArray(graph, solution))
        //    cout << "final solution: " << *solutionList << endl;

        // cout << "solution weight:" << localSearch->getSolutionWeight() << " algorithm weight"<< cost << endl;

        // cout << "execution time: " << std::chrono::duration_cast<std::chrono::microseconds>(done-started).count() << " microseconds" << endl;
        //ILS testing
        IteratedLocalSearch* ILS = new IteratedLocalSearch(graph);

            
        auto started = std::chrono::high_resolution_clock::now();
        NodeBitArray solution = ILS->startResolve();
        auto done = std::chrono::high_resolution_clock::now();
        

        if(vertexCoverValidityEdgescheckBitArray(graph, solution))
            printNodeBitArray(solution, graph->getNumNodes());

        cout << "nodes: "<<graph->getNumNodes() << ", edges: " << graph->getNumEdges()<< ", objevalopit:"<<ILS->getObjectiveEvalForOptimum()<< endl;

        cout <<"[MAIN]: "<< "solution weight:" << ILS->getSolutionWeight() << endl;

        cout <<"[MAIN]: "<< "execution time: " << std::chrono::duration_cast<std::chrono::microseconds>(done-started).count() << " microseconds" << endl;

    }
    cout << endl;

    delete graph;
    delete [] nodesWeight;
}