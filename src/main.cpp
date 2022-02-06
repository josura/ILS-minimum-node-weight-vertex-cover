#include <cstddef>
#include <iostream>
#include <string>
#include <fstream>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>
#include <chrono>

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
        cerr << "input failed";
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
        graph->makeEdgesArray();  // REALLY IMPORTANT
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
        LocalSearch* localSearch = new LocalSearch(graph);
        cout << "greedy initial solution weight:" << localSearch->getSolutionWeight() << endl;

        auto started = std::chrono::high_resolution_clock::now();
        NodeBitArray solution = localSearch->startResolveOptimized();
        auto done = std::chrono::high_resolution_clock::now();
        
        NodeList* solutionList = nodeBitArrayToList(solution,graph->getNumNodes());

        //cout << "final solution: " << *solutionList << endl;

        cout << "solution weight:" << localSearch->getSolutionWeight() << endl;

        cout << "execution time: " << std::chrono::duration_cast<std::chrono::microseconds>(done-started).count() << " microseconds" << endl;

        input.close();
    }
    cout << endl;

    delete graph;
    delete [] nodesWeight;
}