
#include "IteratedLocalSearch.h"
#include "LocalSearch.h"
#include "utilities.h"
#include <algorithm>
#include <array>
#include <cmath>
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
    this->maxNumberOfIterationsWithouthImprovement = numberofIterationsWithouthImprovement;
    this->localSearchMethod = new LocalSearch(_graph);
}

IteratedLocalSearch::~IteratedLocalSearch(){
    delete [] solution;
}  

NodeBitArray IteratedLocalSearch::simplePerturbationSolution(NodeBitArray solution){
    NodeBitArray ret = new bool[graph->getNumNodes()];
    std::copy(solution, solution+graph->getNumNodes(), ret);
    uint nodesChanged =0;
    uint maximumNodesChanged = std::floor((graph->getNumNodes() / 100.0) * power);

    for (uint i = 0; nodesChanged<maximumNodesChanged && i< graph->getNumNodes(); i++) {
        if(randomRealNumber(0, 100) <= rarity){
            ret[i]=!ret[i];
            nodesChanged++;
        }
    }

    return randomSmartSolutionBitArrayFromPartial(graph, ret);
}

void IteratedLocalSearch::adjustPerturbationPowerAndRarity(uint iterationWithoutImprovement, uint iterationWithImprovement){
    if (power>100) {
        power = 100;
    }
    if (rarity>100) {
        rarity = 100;
    }
}


void IteratedLocalSearch::adjustProbabilityOfBadSolutions(uint iterationWithoutImprovement, uint iterationWithImprovement,uint badSolutionsAccepted){

}

NodeBitArray IteratedLocalSearch::startResolve(){
    solution = localSearchMethod->startResolveWithLimit(solutionCost);
    uint iterationWithImprovement = 0 , iterationWithoutImprovement = 0 , badSolutionsAccepted = 0;
    for (uint i = 0; i < maxNumberOfIterations && (iterationWithoutImprovement < maxNumberOfIterationsWithouthImprovement); i++) {
        NodeBitArray perturbedSolutionBefore = simplePerturbationSolution(solution);
        std::cout <<" [ILS]: "<< "ILS iteration number " << i<< std::endl;
        std::cout <<" [ILS]: "<< "\t weight of perturbed solution before localSearch : " << graph->costFunction(perturbedSolutionBefore)<< std::endl;

        double perturbedSolutionCost;
        NodeBitArray perturbedSolutionAfter = localSearchMethod->startResolveWithLimit(perturbedSolutionCost,perturbedSolutionBefore); 
        std::cout <<" [ILS]: "<< "\t weight of perturbed solution after  localSearch : " << graph->costFunction(perturbedSolutionAfter)<< std::endl<<std::endl;
        if (!perturbedSolutionAfter) {
            std::cout <<" [ILS]: "<< "maximum number of objective function evaluation reached"<<std::endl;
            if(perturbedSolutionCost<solutionCost ){
                solutionCost = perturbedSolutionCost;
                std::swap(solution,perturbedSolutionAfter);
                delete [] perturbedSolutionAfter;
            }
            std::cout <<" [ILS]: "<< "cost of the solution found: "<< solutionCost <<std::endl;
            return solution;
        }
        double testBadSolution = randomRealNumber(0, 100) ;
        if(testBadSolution < probabilityBadSolution  || perturbedSolutionCost<solutionCost ){
            if (testBadSolution < probabilityBadSolution && perturbedSolutionCost>=solutionCost) {
                iterationWithoutImprovement++;
                badSolutionsAccepted++;
            } else {
                iterationWithImprovement++;
            }
            solutionCost = perturbedSolutionCost;
            std::swap(solution,perturbedSolutionAfter);
            delete [] perturbedSolutionAfter;
        } else {
            iterationWithoutImprovement++;
            delete [] perturbedSolutionAfter;
        }
        adjustPerturbationPowerAndRarity(iterationWithoutImprovement, iterationWithImprovement);
        adjustProbabilityOfBadSolutions(iterationWithoutImprovement, iterationWithImprovement, badSolutionsAccepted);
    }
    if (iterationWithoutImprovement >= maxNumberOfIterationsWithouthImprovement) {
        std::cout <<" [ILS]: "<< "maximum number of iteration without improvement reached for ILS, TERMINATED."<<std::endl;
    }
    return solution;

}

NodeBitArray IteratedLocalSearch::getSolution()const{
    return solution;
}

double IteratedLocalSearch::getSolutionWeight()const{
    return graph->costFunction(solution);
}