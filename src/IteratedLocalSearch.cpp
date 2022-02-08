
#include "IteratedLocalSearch.h"
#include "LocalSearch.h"
#include "utilities.h"
#include <array>
IteratedLocalSearch::IteratedLocalSearch(WeightedVertexGraph* _graph,uint numberOfIterations, uint numberofIterationsWithouthImprovement, double startingPower, double startingRarity, double startingProbabilityBadSolution){
    this->graph = _graph;
    this->power = startingPower;
    this->rarity = startingRarity;
    this->probabilityBadSolution = startingProbabilityBadSolution;
    this->solution = nullptr;
    this->graph->makeEdgesArray();
    this->maxNumberOfIterations = numberOfIterations;
    this->localSearchMethod = new LocalSearch(_graph);
}

IteratedLocalSearch::~IteratedLocalSearch(){
    delete [] solution;
}  

NodeBitArray IteratedLocalSearch::simplePerturbationSolution(NodeBitArray solution){
    NodeBitArray ret = new bool[graph->getNumNodes()];

    return ret;
}

double IteratedLocalSearch::adjustPerturbationPowerAndRarity(){

}


double IteratedLocalSearch::adjustProbabilityOfBadSolutions(){

}

NodeBitArray IteratedLocalSearch::startResolve(){
    solution = localSearchMethod->startResolveOptimized(solutionCost);
    for (uint i = 0; i < maxNumberOfIterations; i++) {
        NodeBitArray perturbedSolutionBefore = simplePerturbationSolution(solution);
        double perturbedSolutionCost;
        NodeBitArray perturbedSolutionAfter = localSearchMethod->startResolveOptimized(perturbedSolutionCost); 
        
        if(perturbedSolutionCost<solutionCost || ){
            solutionCost = perturbedSolutionCost;
            std::swap(solution,perturbedSolutionAfter);
            delete [] perturbedSolutionAfter;
        } else {
            delete [] perturbedSolutionAfter;
        }
        
    }

}

NodeBitArray IteratedLocalSearch::getSolution()const{

}

double getSolutionWeight()const{

}