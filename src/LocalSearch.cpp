#include "LocalSearch.h"
#include "utilities.h"

LocalSearch::LocalSearch(WeightedVertexGraph* _graph, uint numberOfIterations){
    this->graph = _graph;
    this->solution = greedySolutionBitArray(graph);
    /*std::cout << "greedy solution: ";
    for (int i = 0; i< graph->getNumNodes(); i++) {
            if (solution[i]) {
                std::cout << " true ";
            }else{
                std::cout << " false ";
            }
        }
        std::cout << std::endl;*/
    this->numberOfIterations = numberOfIterations;
}

LocalSearch::~LocalSearch(){
    delete [] this->solution;
}

NodeBitArray LocalSearch::getSolution()const{
    return solution;
}

double LocalSearch::getSolutionWeight()const{
    return costFunctionBitArray(graph, solution);
}

NodeBitArray LocalSearch::startResolveOptimized(){
    uint numNodesGraph = graph->getNumNodes();
    NodeBitArray candidateSolution = new bool[numNodesGraph];
    double currentMinimumWeight = getSolutionWeight();
    for (int numIter = 0; numIter < numberOfIterations; numIter++) {
        double previousMinimumWeight = currentMinimumWeight;
        std::copy(solution, solution + numNodesGraph, candidateSolution);
        for(uint i = 0; i < numNodesGraph ; i++){
            //removing a node from the solution and seeing if it is valid and better than the current solution
            if(candidateSolution[i]) {
                candidateSolution[i] = false;
                if (vertexCoverValidityEdgescheckBitArray(graph,candidateSolution) && (currentMinimumWeight > costFunctionBitArray(graph, candidateSolution))) {
                    //better solution
                    std::copy(candidateSolution, candidateSolution + numNodesGraph, solution);
                    currentMinimumWeight = costFunctionBitArray(graph, candidateSolution);
                }
                //reestablishing the original solution
                candidateSolution[i] = true;


                // substitution of 1 element in the solution with another element, only one considered for candidate in the neighborhood
                for (uint j = i+1; j < numNodesGraph ; j++) {
                    if(!candidateSolution[j]){
                        candidateSolution[j] = true;
                        candidateSolution[i] = false;

                        if (vertexCoverValidityEdgescheckBitArray(graph,candidateSolution) && (currentMinimumWeight > costFunctionBitArray(graph, candidateSolution))) {
                            //better solution with sostitution of node i with node j
                            std::copy(candidateSolution, candidateSolution + numNodesGraph, solution);
                            currentMinimumWeight = costFunctionBitArray(graph, candidateSolution);
                        }

                        candidateSolution[j] = false;
                        candidateSolution[i] = true;
                    }
                }
            }
        }
        if(previousMinimumWeight <= currentMinimumWeight){
            //local minimum
            delete [] candidateSolution;
            return solution;
        }
    }
    delete [] candidateSolution;
    return solution;
}   