
#include "IteratedLocalSearch.h"
#include "LocalSearch.h"
#include "utilities.h"
#include <array>
#include <iostream>
#include <ostream>
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

void IteratedLocalSearch::adjustPerturbationPowerAndRarity(uint iterationWithoutImprovement, uint iterationWithImprovement){

}


void IteratedLocalSearch::adjustProbabilityOfBadSolutions(uint iterationWithoutImprovement, uint iterationWithImprovement,uint badSolutionsAccepted){

}

NodeBitArray IteratedLocalSearch::startResolve(){
    solution = localSearchMethod->startResolveWithLimit(solutionCost);
    for (uint i = 0; i < maxNumberOfIterations; i++) {
        NodeBitArray perturbedSolutionBefore = simplePerturbationSolution(solution);
        perturbedSolutionBefore = randomSmartSolutionBitArrayFromPartial(graph, perturbedSolutionBefore);
        double perturbedSolutionCost;
        NodeBitArray perturbedSolutionAfter = localSearchMethod->startResolveWithLimit(perturbedSolutionCost,perturbedSolutionBefore); 
        if (!perturbedSolutionAfter) {
            std::cout << "maximum number of objective function evaluation reached"<<std::endl;
            if(perturbedSolutionCost<solutionCost ){
                solutionCost = perturbedSolutionCost;
                std::swap(solution,perturbedSolutionAfter);
                delete [] perturbedSolutionAfter;
            }
            std::cout << "cost of the solution found: "<< solutionCost <<std::endl;
            return solution;
        }
        double testBadSolution = randomRealNumber(0, 100) ;
        if(testBadSolution < probabilityBadSolution  || perturbedSolutionCost<solutionCost ){
            solutionCost = perturbedSolutionCost;
            std::swap(solution,perturbedSolutionAfter);
            delete [] perturbedSolutionAfter;
        } else {
            delete [] perturbedSolutionAfter;
        }
        
    }
    return solution;

}

NodeBitArray IteratedLocalSearch::getSolution()const{
    return solution;
}

double IteratedLocalSearch::getSolutionWeight()const{
    return graph->costFunction(solution);
}