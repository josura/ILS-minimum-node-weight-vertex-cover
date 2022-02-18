
#include "IteratedLocalSearch.h"
#include "LocalSearch.h"
#include "utilities.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <ostream>
IteratedLocalSearch::IteratedLocalSearch(WeightedVertexGraph* _graph, double startingPropagationProportion,uint numberOfIterations, uint numberofIterationsWithouthImprovement, double startingPower, double startingRarity, double startingProbabilityBadSolution){
    this->graph = _graph;
    this->power = startingPower;
    this->rarity = startingRarity;
    this->probabilityBadSolution = startingProbabilityBadSolution;
    this->propagationProportion = startingPropagationProportion;
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

    return randomGreedySolutionBitArrayFromPartial(graph, ret);
}

NodeBitArray IteratedLocalSearch::randomPerturbationSolution(NodeBitArray solution){
    NodeBitArray ret = new bool[graph->getNumNodes()];
    std::copy(solution, solution+graph->getNumNodes(), ret);
    uint nodesChanged =0;
    uint maximumNodesChanged = std::floor((graph->getNumNodes() / 100.0) * power);
    NodeList randomNodes = randomVector(0, graph->getNumNodes(),graph->getNumNodes());

    for (auto it = randomNodes.begin(); it != randomNodes.end() && nodesChanged<maximumNodesChanged; it++) {
        if(randomRealNumber(0, 100) <= rarity){
            ret[*it]=!ret[*it];
            nodesChanged++;
        }

    }

    return randomSmartSolutionBitArrayFromPartial(graph, ret);
}

NodeBitArray IteratedLocalSearch::randomDistributedPerturbationSolution(NodeBitArray solution){
    std::vector<uint> nodeInSolution,nodeNotInSolution;
    for (uint i = 0; i < graph->getNumNodes(); i++) {
        if(solution[i]) {
            nodeInSolution.push_back(i);
        } else {
            nodeNotInSolution.push_back(i);
        }
    }
    
    NodeBitArray ret = new bool[graph->getNumNodes()];
    std::copy(solution, solution+graph->getNumNodes(), ret);
    uint nodesChangedRemoved =0, nodesChangedInserted = 0;
    uint maximumNodesChanged = std::ceil((graph->getNumNodes() / 100.0) * power);
    uint maxNodesRemoved = std::floor((maximumNodesChanged/100.0) * propagationProportion);
    uint maxNodesInserted = maximumNodesChanged - maxNodesRemoved;
    NodeList randomNodesRemoved = randomVector(0, nodeInSolution.size(),nodeInSolution.size());
    NodeList randomNodesInserted = randomVector(0, nodeNotInSolution.size(),nodeNotInSolution.size());

    for (auto it = randomNodesRemoved.begin(); it != randomNodesRemoved.end() && nodesChangedRemoved<maxNodesRemoved; it++) {
        if(randomRealNumber(0, 100) <= rarity){
            ret[nodeInSolution[*it]]= false;   //!ret[nodeInSolution[*it]];
            nodesChangedRemoved++;
        }

    }

    for (auto it = randomNodesInserted.begin(); it != randomNodesInserted.end() && nodesChangedInserted<maxNodesInserted; it++) {
        if(randomRealNumber(0, 100) <= rarity){
            ret[nodeNotInSolution[*it]]=true;    //!ret[nodeNotInSolution[*it]];
            nodesChangedInserted++;
        }

    }

    return randomGreedySolutionBitArrayFromPartial(graph, ret);
}


NodeBitArray IteratedLocalSearch::randomDistributedPerturbationAndRemoves(NodeBitArray solution){
    return removeRedundantNodesBest(graph, randomDistributedPerturbationSolution(solution));
}

void IteratedLocalSearch::adjustPerturbationPowerAndRarity(int iterationWithoutImprovement, int iterationWithImprovement){
    power += power*(((iterationWithoutImprovement-iterationWithImprovement)*graph->getMaxDegree()*1.0)/(maxNumberOfIterationsWithouthImprovement *  graph->getNumNodes()));
    if (power>20) {
        power = 20;
    }
    if (power<0.1) {
        power = 0.1;
    }

    rarity += rarity*(((iterationWithoutImprovement-iterationWithImprovement)*graph->getAverageDegree())/(maxNumberOfIterationsWithouthImprovement *  graph->getNumNodes()));
    if (rarity>100) {
        rarity = 100;
    }
    if (rarity<0.1) {
        rarity = 0;
    }
}


void IteratedLocalSearch::adjustProbabilityOfBadSolutions(int iterationWithoutImprovement, int iterationWithImprovement,int badSolutionsAccepted){
    probabilityBadSolution += probabilityBadSolution*(((iterationWithoutImprovement-iterationWithImprovement-badSolutionsAccepted)*graph->getAverageDegree())/(maxNumberOfIterationsWithouthImprovement *  graph->getNumNodes()));
    if (probabilityBadSolution > 10) {
        probabilityBadSolution = 10;
    }
    if (probabilityBadSolution < 0) {
        probabilityBadSolution = 0;
    }
}

NodeBitArray IteratedLocalSearch::startResolve(){
    solution = localSearchMethod->startResolveFinal(solutionCost);
    objectiveEvalForOptimum = localSearchMethod->getCurrentNumberOfObjectiveFunctionEval();
    solutionCost = graph->costFunction(solution);
    uint iterationWithImprovement = 0 , iterationWithoutImprovement = 0 , badSolutionsAccepted = 0;
    for (uint i = 0; i < maxNumberOfIterations && (iterationWithoutImprovement < maxNumberOfIterationsWithouthImprovement); i++) {
        NodeBitArray perturbedSolutionBefore = randomDistributedPerturbationSolution(solution);
        //NodeBitArray perturbedSolutionBefore = randomDistributedPerturbationAndRemoves(solution);
        
        std::cout <<" [ILS]: "<< "ILS iteration number " << i<< std::endl;
        std::cout <<" [ILS]: "<< "\t weight of perturbed solution before localSearch : " << graph->costFunction(perturbedSolutionBefore)<< std::endl;

        double perturbedSolutionCost;
        NodeBitArray perturbedSolutionAfter = localSearchMethod->startResolveFinal(perturbedSolutionCost,perturbedSolutionBefore); 
        if (!perturbedSolutionAfter) {
            perturbedSolutionAfter = localSearchMethod->getSolution();
            std::cout <<" [ILS]: "<< "maximum number of objective function evaluation reached"<<std::endl;
            if(perturbedSolutionCost<solutionCost ){
                solutionCost = perturbedSolutionCost;
                std::swap(solution,perturbedSolutionAfter);
                delete [] perturbedSolutionAfter;
            }
            std::cout <<" [ILS]: "<< "cost of the solution found: "<< solutionCost <<std::endl;
            //TESTING
            std::cout << "remove Globali : " << localSearchMethod->globalRemove << ", swap out :"<< localSearchMethod->globalSwap << ", swap in" << localSearchMethod->globalHybrid<<std::endl;
            //TESTING
            return solution;
        }
        std::cout <<" [ILS]: "<< "\t weight of perturbed solution after  localSearch : " << graph->costFunction(perturbedSolutionAfter)<< std::endl<<std::endl;
        std::cout <<" [ILS]: "<< "cost function evaluation done : " << localSearchMethod->getCurrentNumberOfObjectiveFunctionEval()<<std::endl;
        double testBadSolution = randomRealNumber(0, 100) ;
        if(testBadSolution < probabilityBadSolution  || perturbedSolutionCost<solutionCost ){
            if (testBadSolution < probabilityBadSolution && perturbedSolutionCost>=solutionCost) {
                iterationWithoutImprovement++;
                badSolutionsAccepted++;
            } else {
                objectiveEvalForOptimum = localSearchMethod->getCurrentNumberOfObjectiveFunctionEval();
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
        std::cout <<" [ILS]: "<< "power of perturbation adjusted: "<<power<<std::endl;
        std::cout <<" [ILS]: "<< "rarity of perturbation adjusted: "<<rarity<<std::endl;
        std::cout <<" [ILS]: "<< "probability of bad solution adjusted: "<<probabilityBadSolution<<std::endl<<std::endl;
    }
    if (iterationWithoutImprovement >= maxNumberOfIterationsWithouthImprovement) {
        std::cout <<" [ILS]: "<< "maximum number of iteration without improvement reached for ILS, TERMINATED."<<std::endl;
    }

    //TESTING
    std::cout << "remove Globali : " << localSearchMethod->globalRemove << ", swap :"<< localSearchMethod->globalSwap << ", hybrid" << localSearchMethod->globalHybrid<<std::endl;
    //TESTING
    return solution;

}

NodeBitArray IteratedLocalSearch::getSolution()const{
    return solution;
}

double IteratedLocalSearch::getSolutionWeight()const{
    return graph->costFunction(solution);
}


double IteratedLocalSearch::getPower()const{
    return power;
}

double IteratedLocalSearch::getRarity()const{
    return rarity;
}

double IteratedLocalSearch::getProbabilityBadSolution()const{
    return probabilityBadSolution;
}

int IteratedLocalSearch::getObjectiveEvalForOptimum()const{
    return objectiveEvalForOptimum;
}