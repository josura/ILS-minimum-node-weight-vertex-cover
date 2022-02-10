#include "LocalSearch.h"
#include "utilities.h"
#include <cmath>
#include <iostream>
#include <math.h>
#include <new>
#include <ostream>
#include <utility>
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

/*
* Returns: solution pointer or nullptr when the number of objective function evalutation is greater than the max
*/
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


    uint sampleSizeOpt = std::ceil( numNodesGraph / 7.0);


    for (int numIter = 0; numIter < numberOfIterations; numIter++) {
        double previousMinimumWeight = currentMinimumWeight;

        //LOGGING
        //std::cout <<"[LSEA]: "<< "weight for iteration n " << numIter << " : " <<currentMinimumWeight<<std::endl;
        //LOGGING

        std::copy(solution, solution + numNodesGraph, candidateSolution);

        //maybe I should use unorderedSets or something with low complexity for erase-remove operation because I erase a lot
        // for the nodes considered, I should stick with vector because I need to iterate over all elements at the end of every operation
        //REDACTED vector is ok because I need random access with indexes, I could also create a class that memorizes a vector and an unordered set but for now
        //I will stick with only vectors
        // std::unordered_set<uint> nodeInCandidateSolution;
        // std::unordered_set<std::pair<uint, uint>,boost::hash<std::pair<uint, uint>> > swappableInCandidateSolution;
        std::vector<uint> nodeInCandidateSolution;
        std::vector<std::pair<uint, uint>> swappableInCandidateSolution;
        std::vector<uint>  nodeRemovedFromCandidateSolution;
        std::vector<std::pair<uint, uint>> swappableRemovedFromCandidateSolution;
        for (uint i = 0; i < numNodesGraph; i++) {
            if(candidateSolution[i]) {
                nodeInCandidateSolution.push_back(i);
                for (uint j = i+1; j < numNodesGraph; j++) {
                    if(candidateSolution[j]==false){
                        swappableInCandidateSolution.push_back(std::pair<uint, uint>(i,j));
                    }
                }
            }
        }

        int currentNumberSamples = 0;

        //remove nodes from the solution to create a part of the neighbourhood
        while ((currentNumberSamples)<sampleSizeOpt && nodeInCandidateSolution.size()>0) {
            uint indexRemoved1 = randomNumber(0, nodeInCandidateSolution.size()-1);
            uint nodeRemoved = nodeInCandidateSolution[indexRemoved1];
            nodeInCandidateSolution.erase(nodeInCandidateSolution.begin()+indexRemoved1);
            nodeRemovedFromCandidateSolution.push_back(nodeRemoved);
            candidateSolution[nodeRemoved] = false;

            if(graph->vertexCoverValidityEdgescheckBitArray(candidateSolution)){
                double candidateWeight = previousMinimumWeight - nodeWeights[nodeRemoved];
                //TESTING
                // uint testWeight = graph->costFunction(candidateSolution);
                // std::cout << "candidate weight in remove : "<<candidateWeight<<std::endl;
                // std::cout << "test weight in remove : "<<testWeight<<std::endl<<std::endl;
                //TESTING

                if(currentMinimumWeight > candidateWeight){
                    //better solution
                    std::copy(candidateSolution, candidateSolution + numNodesGraph, solution);
                    currentMinimumWeight = candidateWeight;
                }

                currentNumberSamples++;
                currentObjectiveFunEvalution++;

                if (currentObjectiveFunEvalution >= maximumObjectiveFunEvalution) {
                    std::cout <<"[LSEA]: "<< "MAXIMUM NUMBER OF COST FUNCTION EVALUATION REACHED (remove), LEAVING METHOD " << currentMinimumWeight << std::endl;
                    if(previousMinimumWeight <= currentMinimumWeight){
                        //local minimum
                        std::cout <<"[LSEA]: "<< "local minimum weight: " << currentMinimumWeight << std::endl;
                        std::cout <<"[LSEA]: "<< "iteration: " << numIter << std::endl;
                        finalCost = currentMinimumWeight;
                        return nullptr;
                    }
                    delete [] candidateSolution;
                    return nullptr;
                }
            
            }
            candidateSolution[nodeRemoved] = true;

        }

        //reinserting node removed from candidate solution
        for (uint i = 0; i<nodeRemovedFromCandidateSolution.size(); i++) {
            nodeInCandidateSolution.push_back(nodeRemovedFromCandidateSolution[i]);
        }
        nodeRemovedFromCandidateSolution.clear();

        //swap two nodes from the solution
        while ((currentNumberSamples - sampleSizeOpt)<sampleSizeOpt && swappableInCandidateSolution.size() > 0) {
            uint indexSwappable = randomNumber(0, swappableInCandidateSolution.size()-1);
            std::pair<uint, uint> indexes = swappableInCandidateSolution[indexSwappable]; 
            swappableInCandidateSolution.erase(swappableInCandidateSolution.begin()+indexSwappable);
            swappableRemovedFromCandidateSolution.push_back(indexes);
            candidateSolution[indexes.first] = false;
            candidateSolution[indexes.second] = true;

            if(nodeWeights[indexes.first]>nodeWeights[indexes.second] && graph->vertexCoverValidityEdgescheckBitArray(candidateSolution)){
                double candidateWeight = previousMinimumWeight - nodeWeights[indexes.first] + nodeWeights[indexes.second];
                //TESTING
                // uint testWeight = graph->costFunction(candidateSolution);
                // std::cout << "candidate weight in substitute : "<<candidateWeight<<std::endl;
                // std::cout << "test weight in substitute : "<<testWeight<<std::endl<<std::endl;
                //TESTING

                if(currentMinimumWeight > candidateWeight){
                    //better solution
                    std::copy(candidateSolution, candidateSolution + numNodesGraph, solution);
                    currentMinimumWeight = candidateWeight;
                }

                currentNumberSamples++;
                currentObjectiveFunEvalution++;

                if (currentObjectiveFunEvalution >= maximumObjectiveFunEvalution) {
                    std::cout <<"[LSEA]: "<< "MAXIMUM NUMBER OF COST FUNCTION EVALUATION REACHED (substitution), LEAVING METHOD " << currentMinimumWeight << std::endl;
                    if(previousMinimumWeight <= currentMinimumWeight){
                        //local minimum
                        std::cout <<"[LSEA]: "<< "local minimum weight: " << currentMinimumWeight << std::endl;
                        std::cout <<"[LSEA]: "<< "iteration: " << numIter << std::endl;
                        finalCost = currentMinimumWeight;
                        return nullptr;
                    }
                    delete [] candidateSolution;
                    return nullptr;
                }
            
            }
            candidateSolution[indexes.first] = true;
            candidateSolution[indexes.second] = false;
        }

        //reinserting swapped removed from candidate solution
        for (uint i = 0; i<swappableRemovedFromCandidateSolution.size(); i++) {
            swappableInCandidateSolution.push_back(swappableRemovedFromCandidateSolution[i]);
        }
        swappableRemovedFromCandidateSolution.clear();

        //remove and swap nodes(1 remove + 1 swap) from the solution
        //not optimized because if the remove and the swapped out node are the same, the iteration is repeated
        //TODO possibility of an endless loop
        while ((currentNumberSamples - 2*sampleSizeOpt)<sampleSizeOpt && nodeInCandidateSolution.size()>0 && swappableInCandidateSolution.size()>0 ) {
            uint indexSwappable = randomNumber(0, swappableInCandidateSolution.size()-1);
            uint indexRemoved1 = randomNumber(0, nodeInCandidateSolution.size()-1);
            uint nodeRemoved = nodeInCandidateSolution[indexRemoved1];
            std::pair<uint, uint> indexes = swappableInCandidateSolution[indexSwappable];
            //TODO better controls over this iteration, because the possibility of an endless loop is not impossible 
            if(indexes.first != nodeRemoved){
                nodeInCandidateSolution.erase(nodeInCandidateSolution.begin() + indexRemoved1);
                swappableInCandidateSolution.erase(swappableInCandidateSolution.begin()+indexSwappable);
                candidateSolution[nodeRemoved] = false;
                candidateSolution[indexes.first] = false;
                candidateSolution[indexes.second] = true;

                if( (nodeWeights[nodeRemoved]+nodeWeights[indexes.first])>nodeWeights[indexes.second] && graph->vertexCoverValidityEdgescheckBitArray(candidateSolution)){
                    double candidateWeight = previousMinimumWeight - nodeWeights[nodeRemoved] - nodeWeights[indexes.first] + nodeWeights[indexes.second];
                    //TESTING
                    // uint testWeight = graph->costFunction(candidateSolution);
                    // std::cout << "candidate weight in hybrid : "<<candidateWeight<<std::endl;
                    // std::cout << "test weight in hybrid : "<<testWeight<<std::endl<<std::endl;
                    //TESTING
                    if(currentMinimumWeight > candidateWeight){
                        //better solution
                        std::copy(candidateSolution, candidateSolution + numNodesGraph, solution);
                        currentMinimumWeight = candidateWeight;
                    }

                    currentNumberSamples++;
                    currentObjectiveFunEvalution++;

                    if (currentObjectiveFunEvalution >= maximumObjectiveFunEvalution) {
                        std::cout <<"[LSEA]: "<< "MAXIMUM NUMBER OF COST FUNCTION EVALUATION REACHED (remove + substitution), LEAVING METHOD " << currentMinimumWeight << std::endl;
                        if(previousMinimumWeight <= currentMinimumWeight){
                            //local minimum
                            std::cout <<"[LSEA]: "<< "local minimum weight: " << currentMinimumWeight << std::endl;
                            std::cout <<"[LSEA]: "<< "iteration: " << numIter << std::endl;
                            finalCost = currentMinimumWeight;
                            return nullptr;
                        }
                        delete [] candidateSolution;
                        return nullptr;
                    }
                
                }
                candidateSolution[nodeRemoved] = true;
                candidateSolution[indexes.first] = true;
                candidateSolution[indexes.second] = false;
            }
        }
        //TODO generalization with p remove and 1 substitution

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