#include <cstddef>
#include <iostream>
#include <string>
#include <fstream>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

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
        graph->makeEdgesArray();

        cout << *graph << endl;
        size_t numberNodes = graph->getNumNodes();
        NodeList* solution= new NodeList{1,2,4,5,6,7,8,9,12,14,15,17,18,19};
        NodeBitList* bitSolution = new NodeBitList{false,true,true,false,true,true,true,true,true,true,false,false,true,false,true,true,false,true,true,true};
        NodeList* notsolution= new NodeList{1,2,4};
        NodeBitList* notbitsolution = new NodeBitList{false,true,true,false,true,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false};
        //bool edgescovered = vertexCoverValidityEdgescheck(graph, solution);
        bool edgescoveredBit = vertexCoverValidityEdgescheckBit(graph, bitSolution);
        double cost = costFunctionBit(graph, bitSolution);
        cost = costFunction(graph, solution);
        bool notcoveredBit = vertexCoverValidityEdgescheckBit(graph, notbitsolution);

        NodeBitList* firstSolution = greedySolution(graph);
        cout << "first solution :" <<  *firstSolution << endl;

        input.close();
    }
    cout << endl;

    delete graph;
    delete [] nodesWeight;
}