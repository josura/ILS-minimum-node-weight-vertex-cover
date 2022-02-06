#include "heuristics.h"
#include "utilities.h"
#include <algorithm>
#include <utility>

//controlling all edges
/*
 * IMPORTANT always call makeEdgesArray on the graph before this function
 */
bool vertexCoverValidityEdgescheckBitList(WeightedVertexGraph *graph, NodeBitList* nodeSubset){
    //TODO see if edges vector is quicker or not
    std::pair<uint, uint>** edges = graph->getEdgesArray();

    for(int i = 0; i < graph->getNumEdges(); ++i){
        if (!nodeSubset->at(edges[i]->first) && !nodeSubset->at(edges[i]->second) ) return false;
    }
    return true;
}


//controlling all edges
/*
 * IMPORTANT always call makeEdgesArray on the graph before this function
 */
bool vertexCoverValidityEdgescheckBitArray(WeightedVertexGraph *graph, NodeBitArray nodeSubset){
    //TODO see if edges vector is quicker or not
    std::pair<uint, uint>** edges = graph->getEdgesArray();

    for(int i = 0; i < graph->getNumEdges(); ++i){
        if (!nodeSubset[edges[i]->first] && !nodeSubset[edges[i]->second] ) return false;
    }
    return true;
}



double costFunction(WeightedVertexGraph* graph,NodeList* NodeSubset){
    double* nodeweights = graph->getNodeWeights();
    double sum =0;
    for (auto it =NodeSubset->begin(); it != NodeSubset->end(); it++ ) {
        sum += nodeweights[*it];
    }
    return sum;
}

double costFunctionBitList(WeightedVertexGraph* graph,NodeBitList* NodeSubset){
    //std::transform(NodeSubset->begin(), NodeSubset->end(),)  //TODO simil enumerate for indexes, filtering only true values, and passing It to costFunction

    double* nodeweights = graph->getNodeWeights();
    double sum =0;
    uint i =0;
    for (auto it =NodeSubset->begin(); it != NodeSubset->end();i++, it++ ) {
        if (*it) {
            sum += nodeweights[i];
        }
    }
    return sum;
}


double costFunctionBitArray(WeightedVertexGraph* graph,NodeBitArray NodeSubset){

    double* nodeweights = graph->getNodeWeights();
    double sum =0;
    uint numNodes = graph->getNumNodes();
    for (uint i = 0 ; i < numNodes ;i++ ) {
        if (NodeSubset[i]) {
            sum += nodeweights[i];
        }
    }
    return sum;
}

NodeBitList* greedySolutionBitList(WeightedVertexGraph *graph){
    NodeBitList* ret = new NodeBitList(graph->getNumNodes(),false);
    std::pair<uint, uint>** edges = graph->getEdgesArray();
    std::vector<uint> partSol;

    for (int i =0 ; i < graph->getNumEdges(); i++) {
        if (std::count(partSol.begin(), partSol.end(), edges[i]->first ) == 0 && std::count(partSol.begin(), partSol.end(), edges[i]->second ) == 0 ) {
            partSol.push_back(edges[i]->first);
        }
    }

    for (auto it = partSol.begin(); it != partSol.end(); it++) {
        ret->at(*it) = true;
    }

    return ret;
}


NodeBitArray greedySolutionBitArray(WeightedVertexGraph *graph){

    NodeBitArray ret = new bool[graph->getNumNodes()];
    std::pair<uint, uint>** edges = graph->getEdgesArray();
    std::vector<uint> partSol;

    for (int i =0 ; i < graph->getNumEdges(); i++) {
        if (std::count(partSol.begin(), partSol.end(), edges[i]->first ) == 0 && std::count(partSol.begin(), partSol.end(), edges[i]->second ) == 0 ) {
            partSol.push_back(edges[i]->first);
        }
    }

    for (auto it = partSol.begin(); it != partSol.end(); it++) {
        ret[*it] = true;
    }

    return ret;
}

