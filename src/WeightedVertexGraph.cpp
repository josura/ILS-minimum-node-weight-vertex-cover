#include "WeightedVertexGraph.h"
#include "utilities.h"
#include <algorithm>
#include <iostream>
#include <iterator>
#include <ostream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>


bool WeightedVertexGraph::adjNodes(uint node1, uint node2){
    return ( (adjList[node1].find(node2) != adjList[node1].end()) || (adjList[node2].find(node1) != adjList[node2].end())) ;
}


WeightedVertexGraph::WeightedVertexGraph(){
    numberOfNodes = 0;
    this->nodeWeights = nullptr;
    adjList = nullptr;
    adjVector = nullptr;
    
}

WeightedVertexGraph::WeightedVertexGraph(uint numNodes, double* nodeWeights){
    numberOfNodes = numNodes;
    this->nodeWeights = nodeWeights;
    adjList = new std::unordered_set<uint>[numNodes];
    adjVector = new std::vector<uint>[numNodes];
    //edgesVector = new std::vector<std::pair<uint, uint>>;
    
}

WeightedVertexGraph::~WeightedVertexGraph(){
    //for (int i=0; i<numberOfNodes; i++) {
    //    delete this->adjList[i];
    //}
    delete [] this->adjList;
    
}

uint WeightedVertexGraph::degreeOfNode(uint node)const{
    return adjList[node].size();
}

WeightedVertexGraph* WeightedVertexGraph::addEdge(uint node1, uint node2){
    if(node1 >= numberOfNodes || node2 >= numberOfNodes){
        std::cerr << "add edge failed for edges " << node1 << " and " << node2 << std::endl;
    } else if (adjNodes(node1, node2)) {
        //edge already added
    } else {
        numberOfEdges++;
        edgesVector.push_back(std::pair<uint, uint>(node1,node2));
        adjList[node1].insert(node2);
        adjList[node2].insert(node1);
        adjVector[node1].push_back(node2);
        adjVector[node2].push_back(node1);
    }

    return this;
}


std::pair<uint, uint>* WeightedVertexGraph::makeEdgesArray(){
    edgesArray = new std::pair<uint, uint>[numberOfEdges];
    for (int i =0 ; i<numberOfEdges; i++) {
        edgesArray[i] = edgesVector.at(i);
    }
    return edgesArray;
}


// accessory functions

uint WeightedVertexGraph::getNumNodes()const {
    return numberOfNodes;
}

uint WeightedVertexGraph::getNumEdges()const {
    return numberOfEdges;
}

double* WeightedVertexGraph::getNodeWeights()const{
    return nodeWeights;
}

std::string WeightedVertexGraph::getNodeWeightsStr()const{
    std::string stringa = "";
    for (int i = 0; i<numberOfNodes; i++) {
        stringa += std::to_string(nodeWeights[i]) + std::string(" ");
    }
    return stringa;
}


std::unordered_set<uint>* WeightedVertexGraph::getAdjList(uint node)const{
    if(node>=numberOfNodes){
        std::cerr << "trying to get an adjacent list of an unknown node: " << node << ">=" << numberOfNodes << std::endl;
        return NULL;
    }
    return &adjList[node];
}

std::string WeightedVertexGraph::getAdjListStr(uint node)const{
    std::string stringa;
    for(auto it = adjList[node].cbegin(); it != adjList[node].cend();it++){
                stringa += std::to_string(*it) + " ";
            }
    return stringa;
}

std::vector<std::pair<uint, uint>> WeightedVertexGraph::getEdgesVector()const{
    return edgesVector;
}

std::pair<uint, uint>* WeightedVertexGraph::getEdgesArray()const{
    return edgesArray;
}

//optimization

bool WeightedVertexGraph::vertexCoverValidityEdgescheckBitArray(bool* nodeSubset){
    if (!nodeSubset) {
        throw std::invalid_argument("cost function on null pointer");
    }
    for(int i = 0; i < numberOfEdges; ++i){
        if (!nodeSubset[edgesArray[i].first] && !nodeSubset[edgesArray[i].second] ) return false;
    }
    return true;
}


bool WeightedVertexGraph::vertexCoverValidityNodesRemoved(bool* nodeSubset,const std::vector<uint>& nodesRemoved,const std::vector<uint>& nodesAdded){
    if (!nodeSubset) {
        throw std::invalid_argument("cost function on null pointer");
    }
    for (auto it = nodesRemoved.begin(); it != nodesRemoved.end(); it++) {
        for (auto adjIt = adjVector[*it].begin(); adjIt != adjVector[*it].end(); adjIt++) {
            if (nodeSubset[*adjIt] == false) {
                return false;
            }
        }
    }
    return true;
}

double WeightedVertexGraph::costFunction(bool* NodeSubset){
    if (!NodeSubset) {
        throw std::invalid_argument("cost function on null pointer");
    }
    double sum =0;
    for (uint i = 0 ; i < numberOfNodes ;i++ ) {
        if (NodeSubset[i]) {
            sum += nodeWeights[i];
        }
    }
    return sum;
}

double WeightedVertexGraph::getNodeWeight(uint node)const{
    return nodeWeights[node];
}

std::vector<uint> WeightedVertexGraph::getSharedAdjacentNodes(std::vector<uint>& nodes){
    std::unordered_set<uint> set = adjList[nodes[0]];
    for (auto it = nodes.begin()+1; it != nodes.end(); it++) {
        set = intersectionSet(set, adjList[nodes[0]]);
    }
    std::vector<uint> ret(set.begin(),set.end());
    return ret;
}

std::vector<uint> WeightedVertexGraph::getSwappablesIn(std::vector<uint>& nodes,bool* solution){
    std::vector<uint> candidateSwappables = getSharedAdjacentNodes(nodes);
    std::vector<uint> swappables;
    for (auto it = candidateSwappables.begin(); it!=candidateSwappables.end(); it++) {
        if(solution[*it]==false){
            swappables.push_back(*it);
        }
    }
    return swappables;
}
std::vector<uint> WeightedVertexGraph::getSwappablesOut(std::vector<uint>& nodes,bool* solution){
    std::vector<uint> candidateSwappables = getSharedAdjacentNodes(nodes);
    std::vector<uint> swappables;
    for (auto it = candidateSwappables.begin(); it!=candidateSwappables.end(); it++) {
        if(solution[*it]){
            swappables.push_back(*it);
        }
    }
    return swappables;
}


uint WeightedVertexGraph::getMaxDegree()const{
    uint max = 0;
    for (uint i = 0; i<numberOfNodes; i++) {
        if(degreeOfNode(i)>max)max=degreeOfNode(i);
    }
    return max;
}

double WeightedVertexGraph::getAverageDegree()const{
    double sum = 0;
    for (uint i = 0; i<numberOfNodes; i++) {
        sum += degreeOfNode(i);
    }
    return (sum / (2*numberOfNodes));
}