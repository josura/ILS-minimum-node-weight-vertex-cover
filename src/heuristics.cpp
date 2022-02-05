#include "heuristics.h"
#include "utilities.h"
#include <unordered_set>
#include <utility>

//controlling all edges
bool vertexCoverValidityEdgescheckBit(WeightedVertexGraph *graph, NodeBitList* nodeSubset){
    //TODO see if edges vector is quicker or not
    std::pair<uint, uint>** edges = graph->getEdgesArray();

    for(int i = 0; i < graph->getNumEdges(); ++i){
        if (!nodeSubset->at(edges[i]->first) && !nodeSubset->at(edges[i]->second) ) return false;
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

double costFunctionBit(WeightedVertexGraph* graph,NodeBitList* NodeSubset);