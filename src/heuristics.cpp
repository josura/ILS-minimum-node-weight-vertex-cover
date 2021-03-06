#include "heuristics.h"
#include "utilities.h"
#include <algorithm>
#include <cstddef>
#include <functional>
#include <iostream>
#include <iterator>
#include <queue>
#include <unordered_set>
#include <utility>
#include <vector>

//controlling all edges
/*
 * IMPORTANT always call makeEdgesArray on the graph before this function
 */
bool vertexCoverValidityEdgescheckBitList(WeightedVertexGraph *graph, NodeBitList* nodeSubset){
    //TODO see if edges vector is quicker or not
    std::pair<uint, uint>* edges = graph->getEdgesArray();

    for(int i = 0; i < graph->getNumEdges(); ++i){
        if (!nodeSubset->at(edges[i].first) && !nodeSubset->at(edges[i].second) ) return false;
    }
    return true;
}


bool contains(std::vector<uint> vec, const uint & elem)
{
    return any_of(vec.begin(), vec.end(), [&](const auto & x){
        return x == elem;
    });
}

bool vertexCoverValidityEdgescheckNodeList(WeightedVertexGraph *graph, NodeList& nodeSubset){
    //TODO see if edges vector is quicker or not
    std::pair<uint, uint>* edges = graph->getEdgesArray();

    for(int i = 0; i < graph->getNumEdges(); ++i){
        if (!contains(nodeSubset, edges[i].first ) && !contains(nodeSubset, edges[i].second )) {
            return false;
        }
    }
    return true;
}


//controlling all edges
/*
 * IMPORTANT always call makeEdgesArray on the graph before this function
 */
bool vertexCoverValidityEdgescheckBitArray(WeightedVertexGraph *graph, NodeBitArray nodeSubset){
    //TODO see if edges vector is quicker or not
    std::pair<uint, uint>* edges = graph->getEdgesArray();

    for(int i = 0; i < graph->getNumEdges(); ++i){
        if (!nodeSubset[edges[i].first] && !nodeSubset[edges[i].second] ) return false;
    }
    return true;
}

double costFunction(WeightedVertexGraph* graph,NodeList& NodeSubset){
    double* nodeweights = graph->getNodeWeights();
    double sum =0;
    for (auto it =NodeSubset.begin(); it != NodeSubset.end(); it++ ) {
        sum += nodeweights[*it];
    }
    return sum;
}

double costFunction(WeightedVertexGraph* graph,NodeSet* NodeSubset){
    double* nodeweights = graph->getNodeWeights();
    double sum =0;
    for (auto it =NodeSubset->begin(); it != NodeSubset->end(); it++ ) {
        sum += nodeweights[*it];
    }
    return sum;
}

double costFunction(WeightedVertexGraph* graph,NodeBitList* NodeSubset){
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


double costFunction(WeightedVertexGraph* graph,NodeBitArray NodeSubset){

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
    std::pair<uint, uint>* edges = graph->getEdgesArray();
    std::vector<uint> partSol;

    for (int i =0 ; i < graph->getNumEdges(); i++) {
        if (std::count(partSol.begin(), partSol.end(), edges[i].first ) == 0 && std::count(partSol.begin(), partSol.end(), edges[i].second ) == 0 ) {
            partSol.push_back(edges[i].first);
        }
    }

    for (auto it = partSol.begin(); it != partSol.end(); it++) {
        ret->at(*it) = true;
    }

    return ret;
}


NodeSet* greedySolutionBitArray(WeightedVertexGraph *graph, NodeBitArray solution){
    std::pair<uint, uint>* edges = graph->getEdgesArray();
    NodeSet* partSol = new NodeSet;
    for (int i =0 ; i < graph->getNumEdges(); i++) {
        if (std::count(partSol->begin(), partSol->end(), edges[i].first ) == 0 && std::count(partSol->begin(), partSol->end(), edges[i].second ) == 0 ) {
            partSol->insert(edges[i].first);
        }
    }

    for (auto it = partSol->begin(); it != partSol->end(); it++) {
        solution[*it] = true;
    }

    return partSol;
}


//OPTIMIZATION
NodeBitArray removeRedundantNodes(WeightedVertexGraph *graph, NodeBitArray solution){
    std::vector<uint> redundantNodes;
    for (uint i =0 ; i< graph->getNumNodes(); i++) {
        if(solution[i]){
            solution[i] = false;
            std::vector<uint> tmpNodesRemoved{i};
            if (graph->vertexCoverValidityNodesRemoved(solution, tmpNodesRemoved)) {
                redundantNodes.push_back(i);
            }

            solution[i] = true;
        }
    }

    while (redundantNodes.size() > 0) {
        uint index = randomNumber(0, redundantNodes.size());
        solution[redundantNodes[index]] = false;
        redundantNodes.clear();
        for (uint i =0 ; i< graph->getNumNodes(); i++) {
            if(solution[i]){
                solution[i] = false;
                std::vector<uint> tmpNodesRemoved{i};
                if (graph->vertexCoverValidityNodesRemoved(solution, tmpNodesRemoved)) {
                    redundantNodes.push_back(i);
                }

                solution[i] = true;
            }
        }
        
    }

    return solution;
}

NodeBitArray removeRedundantNodesBest(WeightedVertexGraph *graph, NodeBitArray solution){
    double * nodeWeights = graph->getNodeWeights();
    bool candidatesPresence = true;
    

    do {
        double max = 0;
        uint indexMax=-1;
        candidatesPresence = false;
        for (uint i =0 ; i< graph->getNumNodes(); i++) {
            if(solution[i]){
                solution[i] = false;
                std::vector<uint> tmpNodesRemoved{i};
                if (graph->vertexCoverValidityNodesRemoved(solution, tmpNodesRemoved) && nodeWeights[i] > max) {
                    candidatesPresence = true;
                    max = nodeWeights[i];
                    indexMax = i;
                }

                solution[i] = true;
            }
        }
        if (candidatesPresence) {
            solution[indexMax] = false;
        }
        
    } while (candidatesPresence);

    return solution;
}
//OPTIMIZATION


NodeBitArray randomGreedySolutionBitArrayFromPartial(WeightedVertexGraph *graph, NodeBitArray solution){
    std::pair<uint, uint>* edges = graph->getEdgesArray();
    std::unordered_set<std::pair<uint, uint>,boost::hash<std::pair<uint, uint>>> edgesNotCovered;

    for (int i =0 ; i < graph->getNumEdges(); i++) {
        if (!solution[edges[i].first] && !solution[edges[i].second] ) {
            edgesNotCovered.insert(edges[i]);
        }
    }

    do {
        int i = randomNumber(0,edgesNotCovered.size());

        auto it = edgesNotCovered.begin();
        for(uint j=0;j<i;j++){it++;}
        if (it!=edgesNotCovered.end()) {
            std::pair<uint, uint> candidatePair = *it;
            double test = randomRealNumber(0, 100);
            if(test < 33.333){
                solution[candidatePair.first] = true;
                for ( auto it = graph->getAdjList(candidatePair.first)->begin() ; it != graph->getAdjList(candidatePair.first)->end(); it++){
                    std::pair<uint, uint> pair (candidatePair.first,*it);
                    edgesNotCovered.erase(pair);
                }
            }else if(test >= 33.333 && test < 66.666){
                solution[candidatePair.second] = true;
                for ( auto it = graph->getAdjList(candidatePair.second)->begin() ; it != graph->getAdjList(candidatePair.second)->end(); it++){
                    std::pair<uint, uint> pair (candidatePair.second,*it);
                    edgesNotCovered.erase(pair);
                }
            }else {
                solution[candidatePair.first] = true;
                for ( auto it = graph->getAdjList(candidatePair.first)->begin() ; it != graph->getAdjList(candidatePair.first)->end(); it++){
                    std::pair<uint, uint> pair (candidatePair.first,*it);
                    edgesNotCovered.erase(pair);
                }
                solution[candidatePair.second] = true;
                for ( auto it = graph->getAdjList(candidatePair.second)->begin() ; it != graph->getAdjList(candidatePair.second)->end(); it++){
                    std::pair<uint, uint> pair (candidatePair.second,*it);
                    edgesNotCovered.erase(pair);
                }        
            }
        }

    }while (edgesNotCovered.size()>0);

    bool test = graph->vertexCoverValidityEdgescheckBitArray(solution);

    //return removeRedundantNodesBest(graph, solution);
    return solution;
}

NodeSet* randomGreedySolutionBitArray(WeightedVertexGraph *graph, NodeBitArray solution){
    for (uint i = 0; i < graph->getNumNodes(); i++) {
        solution[i] = false;
    }

    solution = randomGreedySolutionBitArrayFromPartial(graph,solution);

    return nodeBitArrayToSet(solution, graph->getNumNodes());
}
struct Comp{
    bool operator()(const std::pair<double, uint> &a, const std::pair<double, uint> &b) { 
        return a.first < b.first;
    }
};
bool pairSecondFunctorIncr(std::pair<double, uint> const &a, std::pair<double, uint> const &b) { 
       return a.second < b.second;
}

//taking nodes that have an high (degree of edges not covered / weight)
NodeBitArray randomSmartSolutionBitArrayFromPartial(WeightedVertexGraph *graph, NodeBitArray solution){
    std::pair<uint, uint>* edges = graph->getEdgesArray();
    std::unordered_set<std::pair<uint, uint>,boost::hash<std::pair<uint, uint>>> edgesNotCovered;
    std::vector<std::pair<double, uint> > nodeDegreesRatioOfNotCoveredEdges;

    for (uint i = 0; i < graph->getNumNodes(); i++) {
        nodeDegreesRatioOfNotCoveredEdges.push_back(std::pair<double, uint>(0,i));
    }

    for (uint i =0 ; i < graph->getNumEdges(); i++) {
        if (!solution[edges[i].first] && !solution[edges[i].second] ) {
            edgesNotCovered.insert(edges[i]);
            nodeDegreesRatioOfNotCoveredEdges[edges[i].first].first+=1;
            nodeDegreesRatioOfNotCoveredEdges[edges[i].second].first+=1;
        }
    }

    for (uint i = 0; i < nodeDegreesRatioOfNotCoveredEdges.size(); i++) {
        nodeDegreesRatioOfNotCoveredEdges[i].first /= graph->getNodeWeight(i); 
    }
    std::make_heap(nodeDegreesRatioOfNotCoveredEdges.begin(),nodeDegreesRatioOfNotCoveredEdges.end(),Comp());


    do {

        std::pair<double , uint> element = nodeDegreesRatioOfNotCoveredEdges.front(); 
        solution[element.second] = true;

        for ( auto it = graph->getAdjList(element.second)->begin() ; it != graph->getAdjList(element.second)->end(); it++){
            std::pair<uint, uint> pair (element.second,*it);
            std::pair<uint, uint> pair2 (*it,element.second);
            edgesNotCovered.erase(pair);
            edgesNotCovered.erase(pair2);
        }

        for (uint i = 0; i < graph->getNumNodes(); i++) {
            nodeDegreesRatioOfNotCoveredEdges[i] = std::pair<double, uint>(0,i);
        }

        for (auto it = edgesNotCovered.begin() ; it != edgesNotCovered.end(); it++) {
            nodeDegreesRatioOfNotCoveredEdges[it->first].first+=1;
            nodeDegreesRatioOfNotCoveredEdges[it->second].first+=1;
        }

        for (uint i = 0; i < nodeDegreesRatioOfNotCoveredEdges.size(); i++) {
            nodeDegreesRatioOfNotCoveredEdges[i].first /= graph->getNodeWeight(i); 
        }

        std::make_heap(nodeDegreesRatioOfNotCoveredEdges.begin(),nodeDegreesRatioOfNotCoveredEdges.end(),Comp());

    }while (edgesNotCovered.size()>0);

    bool test = graph->vertexCoverValidityEdgescheckBitArray(solution);


    //return removeRedundantNodesBest(graph, solution);
    return  solution;
}


NodeBitArray randomGreedySolutionAndRemoves(WeightedVertexGraph *graph, NodeBitArray solution){
    return removeRedundantNodesBest(graph, randomGreedySolutionBitArrayFromPartial(graph, solution));
}

NodeBitArray randomSmartSolutionAndRemoves(WeightedVertexGraph *graph, NodeBitArray solution){
    return removeRedundantNodesBest(graph, randomSmartSolutionBitArrayFromPartial(graph, solution));
}