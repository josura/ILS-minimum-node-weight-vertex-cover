#include "heuristics.h"
#include "utilities.h"
#include <unordered_set>
#include <utility>


//controlling if adjacency list of all nodes contain at least one node in the NodeSubset that represents vertex cover, if nodes are more or less = than edges, this algorithm is slower
bool vertexCoverValidityNodescheckBit(WeightedVertexGraph* graph,NodeBitList* NodeSubset){
    for(uint i = 0; i< graph->getNumNodes();i++){
        //std::unordered_set<uint>* adjListNode = graph->getAdjList(i);
        bool nodeCovered = false;
        if (graph->getAdjList(i)->size() == 0) {
            nodeCovered = true;
        }
        int nodeInd = 0;
        for (auto it = NodeSubset->begin(); it != NodeSubset->end() && !nodeCovered;nodeInd++, it++) {
            if((*it) && graph->adjNodes(i, nodeInd)){
                nodeCovered = true;
            }
        }
        if(!nodeCovered){
            return false;
        }
    }

    return true;
}

//controlling all edges
bool vertexCoverValidityEdgescheckBit(WeightedVertexGraph *graph, NodeBitList* nodeSubset){
    //TODO see if edges vector is quicker or not
    std::pair<uint, uint>** edges = graph->getEdgesArray();

    for(int i = 0; i < graph->getNumEdges(); ++i){
        if (!nodeSubset->at(edges[i]->first) && !nodeSubset->at(edges[i]->second) ) return false;
    }
    return true;
}

bool vertexCoverValidityNodescheck(WeightedVertexGraph* graph,NodeList* NodeSubset){
    for(uint i = 0; i< graph->getNumNodes();i++){
        //std::unordered_set<uint>* adjListNode = graph->getAdjList(i);
        bool nodeCovered = false;
        if (graph->getAdjList(i)->size() == 0) {
            nodeCovered = true;
        }
        for (auto it = NodeSubset->begin(); it != NodeSubset->end() && !nodeCovered; it++) {
            if(graph->adjNodes(i, *it)){
                nodeCovered = true;
            }
        }
        if(!nodeCovered){
            return false;
        }
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