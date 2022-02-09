#include "LocalSearch.h"
#include "utilities.h"
#include <iostream>
#include <new>
#include <ostream>
#include <vector>

LocalSearch::LocalSearch(WeightedVertexGraph* _graph, uint numberOfIterations, double samplingFactor){
    this->graph = _graph;
    this->samplingFactor = samplingFactor;
    this->solution = new bool[graph->getNumNodes()];
    this->graph->makeEdgesArray();
    this->solutionSet = *greedySolutionBitArray(_graph,this->solution);
    this->numberOfIterations = numberOfIterations;
    this->currentObjectiveFunEvalution = 0;
}

LocalSearch::~LocalSearch(){
    delete [] this->solution;
}

void LocalSearch::setSolution(NodeBitArray solution){
    this->solution = solution;
}

NodeBitArray LocalSearch::getSolution()const{
    return solution;
}


NodeSet LocalSearch::getSolutionSet()const{
    return solutionSet;
}

double LocalSearch::getSolutionWeight(){
    currentObjectiveFunEvalution++;
    return graph->costFunction( solution);
}

NodeBitArray LocalSearch::startResolve(){
    uint numNodesGraph = graph->getNumNodes();
    NodeBitArray candidateSolution = new bool[numNodesGraph];
    NodeSet* candidateSolutionSet = new NodeSet;
    double currentMinimumWeight = getSolutionWeight();
    for (int numIter = 0; numIter < numberOfIterations; numIter++) {
        double previousMinimumWeight = currentMinimumWeight;
        std::copy(solution, solution + numNodesGraph, candidateSolution);
        for(uint i = 0; i < numNodesGraph ; i++){
            //removing a node from the solution and seeing if it is valid and better than the current solution
            if(candidateSolution[i]) {
                candidateSolution[i] = false;
                candidateSolutionSet->erase(i);
                if (vertexCoverValidityEdgescheckBitArray(graph,candidateSolution) ) {
                    double candidateWeight = costFunction(graph, candidateSolution);
                    //double candidateWeight = costFunction(graph, candidateSolutionSet);
                    if(currentMinimumWeight > candidateWeight){
                        //better solution
                        std::copy(candidateSolution, candidateSolution + numNodesGraph, solution);
                        solutionSet = *candidateSolutionSet;
                        currentMinimumWeight = candidateWeight;
                    }
                }
                //reestablishing the original solution
                candidateSolution[i] = true;
                candidateSolutionSet->insert(i);


                // substitution of 1 element in the solution with another element, only one considered for candidate in the neighborhood
                for (uint j = i+1; j < numNodesGraph ; j++) {
                    if(!candidateSolution[j]){
                        candidateSolution[j] = true;
                        candidateSolutionSet->insert(j);
                        candidateSolution[i] = false;
                        candidateSolutionSet->erase(i);
                        if (vertexCoverValidityEdgescheckBitArray(graph,candidateSolution) ) {
                            double candidateWeight = costFunction(graph, candidateSolution);
                            //double candidateWeight = costFunction(graph, candidateSolutionSet);
                            if(currentMinimumWeight > candidateWeight){    
                                //better solution with sostitution of node i with node j
                                std::copy(candidateSolution, candidateSolution + numNodesGraph, solution);
                                solutionSet = *candidateSolutionSet;
                                currentMinimumWeight = candidateWeight;
                            }
                        }

                        candidateSolution[j] = false;
                        candidateSolutionSet->erase(j);
                        candidateSolution[i] = true;
                        candidateSolutionSet->insert(i);
                    }
                }
            }
        }
        if(previousMinimumWeight <= currentMinimumWeight){
            //local minimum
            std::cout << "iteration: " << numIter << std::endl;
            delete [] candidateSolution;
            delete candidateSolutionSet;
            return solution;
        }
    }
    delete [] candidateSolution;
    delete candidateSolutionSet;
    return solution;
}   


NodeBitArray LocalSearch::startResolveStrange(){
    uint numNodesGraph = graph->getNumNodes();
    double* nodeWeights = graph->getNodeWeights();
    NodeBitArray candidateSolution = new bool[numNodesGraph];
    double currentMinimumWeight = getSolutionWeight();
    for (int numIter = 0; numIter < numberOfIterations; numIter++) {
        double previousMinimumWeight = currentMinimumWeight;

        //LOGGING
        std::cout << "weight for iteration n " << numIter << " : " <<currentMinimumWeight<<std::endl;
        //LOGGING

        std::copy(solution, solution + numNodesGraph, candidateSolution);
        for(uint i = 0; i < numNodesGraph ; i++){
            //removing a node from the solution and seeing if it is valid and better than the current solution
            if(candidateSolution[i]) {
                candidateSolution[i] = false;
                if (vertexCoverValidityEdgescheckBitArray(graph,candidateSolution) ) {
                    double candidateWeight = currentMinimumWeight - nodeWeights[i];
                    if(currentMinimumWeight > candidateWeight){
                        //better solution
                        std::copy(candidateSolution, candidateSolution + numNodesGraph, solution);
                        currentMinimumWeight = candidateWeight;
                    }
                }
                //reestablishing the original solution
                candidateSolution[i] = true;


                // substitution of 1 element in the solution with another element, only one considered for candidate in the neighborhood
                for (uint j = i+1; j < numNodesGraph ; j++) {
                    if(!candidateSolution[j]){
                        candidateSolution[j] = true;
                        candidateSolution[i] = false;
                        if (vertexCoverValidityEdgescheckBitArray(graph,candidateSolution) ) {
                            double candidateWeight = currentMinimumWeight - nodeWeights[i] + nodeWeights[j] ;
                            if(currentMinimumWeight > candidateWeight){    
                                //better solution with sostitution of node i with node j
                                std::copy(candidateSolution, candidateSolution + numNodesGraph, solution);
                                currentMinimumWeight = candidateWeight;
                            }
                        }

                        candidateSolution[j] = false;
                        candidateSolution[i] = true;
                    }
                }
            }
        }
        if(previousMinimumWeight <= currentMinimumWeight){
            //local minimum
            std::cout << "local minimum weight: " << currentMinimumWeight << std::endl;
            std::cout << "iteration: " << numIter << std::endl;
            delete [] candidateSolution;
            return solution;
        }
    }
    delete [] candidateSolution;
    return solution;
}   

NodeBitArray LocalSearch::startResolveOptimized(){
    uint numNodesGraph = graph->getNumNodes();
    double* nodeWeights = graph->getNodeWeights();
    NodeBitArray candidateSolution = new bool[numNodesGraph];
    double currentMinimumWeight = getSolutionWeight();
    for (int numIter = 0; numIter < numberOfIterations; numIter++) {
        double currentIterationMinimumWeight = currentMinimumWeight;
        double previousMinimumWeight = currentMinimumWeight;

        //LOGGING
        std::cout << "weight for iteration n " << numIter << " : " <<currentMinimumWeight<<std::endl;
        //LOGGING

        std::copy(solution, solution + numNodesGraph, candidateSolution);
        for(uint i = 0; i < numNodesGraph ; i++){
            //removing a node from the solution and seeing if it is valid and better than the current solution
            if(candidateSolution[i]) {
                candidateSolution[i] = false;
                if (graph->vertexCoverValidityEdgescheckBitArray(candidateSolution) ) {
                    double candidateWeight = currentIterationMinimumWeight - nodeWeights[i];
                    if(currentMinimumWeight > candidateWeight){
                        //better solution
                        std::copy(candidateSolution, candidateSolution + numNodesGraph, solution);
                        currentMinimumWeight = candidateWeight;
                    }
                }
                //reestablishing the original solution
                candidateSolution[i] = true;


                // substitution of 1 element in the solution with another element, only one considered for candidate in the neighborhood
                for (uint j = i+1; j < numNodesGraph ; j++) {
                    if(!candidateSolution[j]){
                        candidateSolution[j] = true;
                        candidateSolution[i] = false;
                        if (graph->vertexCoverValidityEdgescheckBitArray(candidateSolution) ) {
                            double candidateWeight = currentIterationMinimumWeight - nodeWeights[i] + nodeWeights[j] ;
                            if(currentMinimumWeight > candidateWeight){    
                                //better solution with sostitution of node i with node j
                                std::copy(candidateSolution, candidateSolution + numNodesGraph, solution);
                                currentMinimumWeight = candidateWeight;
                            }
                        }

                        candidateSolution[j] = false;
                        candidateSolution[i] = true;
                    }
                }
            }
        }
        if(previousMinimumWeight <= currentMinimumWeight){
            //local minimum
            std::cout << "local minimum weight: " << currentMinimumWeight << std::endl;
            std::cout << "iteration: " << numIter << std::endl;
            delete [] candidateSolution;
            return solution;
        }
    }
    delete [] candidateSolution;
    return solution;
}   


NodeBitArray LocalSearch::startResolveOptimized(double &finalCost ,NodeBitArray startSolution){
    uint numNodesGraph = graph->getNumNodes();
    double* nodeWeights = graph->getNodeWeights();
    NodeBitArray candidateSolution = new bool[numNodesGraph];
    double currentMinimumWeight = getSolutionWeight();
    //changing solution pointer to startSolution, freeing memory from solution should be done outside of the class
    if(startSolution)this->solution = startSolution;
    else{this->solution = new bool[numNodesGraph];}

    for (int numIter = 0; numIter < numberOfIterations; numIter++) {
        double currentIterationMinimumWeight = currentMinimumWeight;
        double previousMinimumWeight = currentMinimumWeight;

        //LOGGING
        std::cout << "weight for iteration n " << numIter << " : " <<currentMinimumWeight<<std::endl;
        //LOGGING

        std::copy(solution, solution + numNodesGraph, candidateSolution);
        for(uint i = 0; i < numNodesGraph ; i++){
            //removing a node from the solution and seeing if it is valid and better than the current solution
            if(candidateSolution[i]) {
                candidateSolution[i] = false;
                if (graph->vertexCoverValidityEdgescheckBitArray(candidateSolution) ) {
                    double candidateWeight = currentIterationMinimumWeight - nodeWeights[i];
                    if(currentMinimumWeight > candidateWeight){
                        //better solution
                        std::copy(candidateSolution, candidateSolution + numNodesGraph, solution);
                        currentMinimumWeight = candidateWeight;
                    }
                }
                //reestablishing the original solution
                candidateSolution[i] = true;


                // substitution of 1 element in the solution with another element, only one considered for candidate in the neighborhood
                for (uint j = i+1; j < numNodesGraph ; j++) {
                    //if(!candidateSolution[j]){   //C++ has gone mad because this condition passes sometimes when It should not
                    if(candidateSolution[j] == false){
                        candidateSolution[j] = true;
                        candidateSolution[i] = false;
                        if (graph->vertexCoverValidityEdgescheckBitArray(candidateSolution) ) {
                            double candidateWeight = currentIterationMinimumWeight - nodeWeights[i] + nodeWeights[j] ;
                            if(currentMinimumWeight > candidateWeight){    
                                //better solution with sostitution of node i with node j
                                std::copy(candidateSolution, candidateSolution + numNodesGraph, solution);
                                currentMinimumWeight = candidateWeight;
                            }
                        }

                        candidateSolution[j] = false;
                        candidateSolution[i] = true;
                    }
                }
            }
        }
        if(previousMinimumWeight <= currentMinimumWeight){
            //local minimum
            std::cout << "local minimum weight: " << currentMinimumWeight << std::endl;
            std::cout << "iteration: " << numIter << std::endl;
            delete [] candidateSolution;
            finalCost = currentMinimumWeight;
            return solution;
        }
    }
    finalCost = currentMinimumWeight;
    delete [] candidateSolution;
    return solution;
}   


NodeBitArray LocalSearch::startResolveWithLimit(double &finalCost ,NodeBitArray startSolution){
    uint numNodesGraph = graph->getNumNodes();
    double* nodeWeights = graph->getNodeWeights();
    NodeBitArray candidateSolution = new bool[numNodesGraph];
    //changing solution pointer to startSolution, freeing memory from solution should be done outside of the class
    if(startSolution)this->solution = startSolution;
    else{ 
        solution = new bool[numNodesGraph];
        this->solution = randomSmartSolutionBitArrayFromPartial(graph, solution);
    }
    double currentMinimumWeight = getSolutionWeight();

    for (int numIter = 0; numIter < numberOfIterations; numIter++) {
        double currentIterationMinimumWeight = currentMinimumWeight;
        double previousMinimumWeight = currentMinimumWeight;

        //LOGGING
        std::cout <<"[LSEA]: "<< "weight for iteration n " << numIter << " : " <<currentMinimumWeight<<std::endl;
        //LOGGING

        std::copy(solution, solution + numNodesGraph, candidateSolution);
        for(uint i = 0; i < numNodesGraph ; i++){
            //removing a node from the solution and seeing if it is valid and better than the current solution
            if(candidateSolution[i]) {
                candidateSolution[i] = false;
                if (graph->vertexCoverValidityEdgescheckBitArray(candidateSolution) ) {
                    double candidateWeight = currentIterationMinimumWeight - nodeWeights[i];
                    //OBJECTIVE FUNCTION EVALUATION
                    currentObjectiveFunEvalution++;
                    //OBJECTIVE FUNCTION EVALUATION
                    
                    if(currentMinimumWeight > candidateWeight){
                        //better solution
                        std::copy(candidateSolution, candidateSolution + numNodesGraph, solution);
                        currentMinimumWeight = candidateWeight;
                    }
                    
                    if (currentObjectiveFunEvalution >= maximumObjectiveFunEvalution) {
                        std::cout <<"[LSEA]: "<< "MAXIMUM NUMBER OF COST FUNCTION EVALUATION REACHED (removal), LEAVING METHOD " << currentMinimumWeight << std::endl;
                        if(previousMinimumWeight <= currentMinimumWeight){
                            //local minimum
                            std::cout <<"[LSEA]: "<< "local minimum weight: " << currentMinimumWeight << std::endl;
                            std::cout <<"[LSEA]: "<< "iteration: " << numIter << std::endl;
                            finalCost = currentMinimumWeight;
                            return solution;
                        }
                        delete [] candidateSolution;
                        return nullptr;
                    }
                    
                }
                //reestablishing the original solution
                candidateSolution[i] = true;


                // substitution of 1 element in the solution with another element, only one considered for candidate in the neighborhood
                for (uint j = i+1; j < numNodesGraph ; j++) {
                    //if(!candidateSolution[j]){   //C++ has gone mad because this condition passes sometimes when It should not
                    if(candidateSolution[j] == false){
                        candidateSolution[j] = true;
                        candidateSolution[i] = false;
                        if (graph->vertexCoverValidityEdgescheckBitArray(candidateSolution) ) {
                            double candidateWeight = currentIterationMinimumWeight - nodeWeights[i] + nodeWeights[j] ;
                            //OBJECTIVE FUNCTION EVALUATION
                            currentObjectiveFunEvalution++;
                            //OBJECTIVE FUNCTION EVALUATION

                            if(currentMinimumWeight > candidateWeight){    
                                //better solution with sostitution of node i with node j
                                std::copy(candidateSolution, candidateSolution + numNodesGraph, solution);
                                currentMinimumWeight = candidateWeight;
                            }

                            if (currentObjectiveFunEvalution >= maximumObjectiveFunEvalution) {
                                std::cout <<"[LSEA]: "<< "MAXIMUM NUMBER OF COST FUNCTION EVALUATION REACHED (substitution), LEAVING METHOD " << currentMinimumWeight << std::endl;
                                if(previousMinimumWeight <= currentMinimumWeight){
                                    //local minimum
                                    std::cout <<"[LSEA]: "<< "local minimum weight: " << currentMinimumWeight << std::endl;
                                    std::cout <<"[LSEA]: "<< "iteration: " << numIter << std::endl;
                                    finalCost = currentMinimumWeight;
                                    return solution;
                                }
                                delete [] candidateSolution;
                                return nullptr;
                            }
                        }

                        candidateSolution[j] = false;
                        candidateSolution[i] = true;
                    }
                }
            }
        }
        if(previousMinimumWeight <= currentMinimumWeight){
            //local minimum
            std::cout <<"[LSEA]: "<< "local minimum weight: " << currentMinimumWeight << std::endl;
            std::cout <<"[LSEA]: "<< "iteration: " << numIter << std::endl;
            delete [] candidateSolution;
            finalCost = currentMinimumWeight;
            return solution;
        }
    }
    finalCost = currentMinimumWeight;
    delete [] candidateSolution;
    return solution;
}   


NodeBitArray LocalSearch::startResolveWithLimitAndSampling(double &finalCost ,NodeBitArray startSolution){
    uint numNodesGraph = graph->getNumNodes();
    double* nodeWeights = graph->getNodeWeights();
    NodeBitArray candidateSolution = new bool[numNodesGraph];
    if(startSolution)this->solution = startSolution;
    else{ 
        solution = new bool[numNodesGraph];
        this->solution = randomSmartSolutionBitArrayFromPartial(graph, solution);
    }
    double currentMinimumWeight = getSolutionWeight();

    for (int numIter = 0; numIter < numberOfIterations; numIter++) {
        double currentIterationMinimumWeight = currentMinimumWeight;
        double previousMinimumWeight = currentMinimumWeight;

        //LOGGING
        std::cout <<"[LSEA]: "<< "weight for iteration n " << numIter << " : " <<currentMinimumWeight<<std::endl;
        //LOGGING

        std::copy(solution, solution + numNodesGraph, candidateSolution);

        std::vector<uint> nodeInCandidateSolution;
        std::vector<uint> nodeNotInCandidateSolution;
        for (uint i = 0; i < numNodesGraph; i++) {
            if(candidateSolution[i]) nodeInCandidateSolution.push_back(i);
            else nodeNotInCandidateSolution.push_back(i);
        }

        for(uint i = 0; i < numNodesGraph ; i++){
            //removing a node from the solution and seeing if it is valid and better than the current solution
            if(candidateSolution[i]) {
                candidateSolution[i] = false;
                if (graph->vertexCoverValidityEdgescheckBitArray(candidateSolution) ) {
                    double candidateWeight = currentIterationMinimumWeight - nodeWeights[i];
                    //OBJECTIVE FUNCTION EVALUATION
                    currentObjectiveFunEvalution++;
                    //OBJECTIVE FUNCTION EVALUATION
                    
                    if(currentMinimumWeight > candidateWeight){
                        //better solution
                        std::copy(candidateSolution, candidateSolution + numNodesGraph, solution);
                        currentMinimumWeight = candidateWeight;
                    }
                    
                    if (currentObjectiveFunEvalution >= maximumObjectiveFunEvalution) {
                        std::cout <<"[LSEA]: "<< "MAXIMUM NUMBER OF COST FUNCTION EVALUATION REACHED (removal), LEAVING METHOD " << currentMinimumWeight << std::endl;
                        if(previousMinimumWeight <= currentMinimumWeight){
                            //local minimum
                            std::cout <<"[LSEA]: "<< "local minimum weight: " << currentMinimumWeight << std::endl;
                            std::cout <<"[LSEA]: "<< "iteration: " << numIter << std::endl;
                            finalCost = currentMinimumWeight;
                            return solution;
                        }
                        delete [] candidateSolution;
                        return nullptr;
                    }
                    
                }
                //reestablishing the original solution
                candidateSolution[i] = true;


                // substitution of 1 element in the solution with another element, only one considered for candidate in the neighborhood
                for (uint j = i+1; j < numNodesGraph ; j++) {
                    //if(!candidateSolution[j]){   //C++ has gone mad because this condition passes sometimes when It should not
                    if(candidateSolution[j] == false){
                        candidateSolution[j] = true;
                        candidateSolution[i] = false;
                        if (graph->vertexCoverValidityEdgescheckBitArray(candidateSolution) ) {
                            double candidateWeight = currentIterationMinimumWeight - nodeWeights[i] + nodeWeights[j] ;
                            //OBJECTIVE FUNCTION EVALUATION
                            currentObjectiveFunEvalution++;
                            //OBJECTIVE FUNCTION EVALUATION

                            if(currentMinimumWeight > candidateWeight){    
                                //better solution with sostitution of node i with node j
                                std::copy(candidateSolution, candidateSolution + numNodesGraph, solution);
                                currentMinimumWeight = candidateWeight;
                            }

                            if (currentObjectiveFunEvalution >= maximumObjectiveFunEvalution) {
                                std::cout <<"[LSEA]: "<< "MAXIMUM NUMBER OF COST FUNCTION EVALUATION REACHED (substitution), LEAVING METHOD " << currentMinimumWeight << std::endl;
                                if(previousMinimumWeight <= currentMinimumWeight){
                                    //local minimum
                                    std::cout <<"[LSEA]: "<< "local minimum weight: " << currentMinimumWeight << std::endl;
                                    std::cout <<"[LSEA]: "<< "iteration: " << numIter << std::endl;
                                    finalCost = currentMinimumWeight;
                                    return solution;
                                }
                                delete [] candidateSolution;
                                return nullptr;
                            }
                        }

                        candidateSolution[j] = false;
                        candidateSolution[i] = true;
                    }
                }
            }
        }
        if(previousMinimumWeight <= currentMinimumWeight){
            //local minimum
            std::cout <<"[LSEA]: "<< "local minimum weight: " << currentMinimumWeight << std::endl;
            std::cout <<"[LSEA]: "<< "iteration: " << numIter << std::endl;
            delete [] candidateSolution;
            finalCost = currentMinimumWeight;
            return solution;
        }
    }
    finalCost = currentMinimumWeight;
    delete [] candidateSolution;
    return solution;
}   

uint LocalSearch::getCurrentNumberOfObjectiveFunctionEval()const{
    return currentObjectiveFunEvalution;
}