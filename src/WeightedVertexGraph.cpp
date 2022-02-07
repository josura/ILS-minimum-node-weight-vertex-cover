#include "WeightedVertexGraph.h"
#include <algorithm>
#include <iostream>
#include <iterator>
#include <ostream>
#include <string>
#include <utility>
#include <vector>


bool WeightedVertexGraph::adjNodes(uint node1, uint node2){
    return ( (adjList[node1].find(node2) != adjList[node1].end()) || (adjList[node2].find(node1) != adjList[node2].end())) ;
}


WeightedVertexGraph::WeightedVertexGraph(){
    numberOfNodes = 0;
    this->nodeWeights = nullptr;
    adjList = nullptr;
    
}

WeightedVertexGraph::WeightedVertexGraph(uint numNodes, double* nodeWeights){
    numberOfNodes = numNodes;
    this->nodeWeights = nodeWeights;
    adjList = new std::unordered_set<uint>[numNodes];
    //edgesVector = new std::vector<std::pair<uint, uint>>;
    for (int i = 0; i < numNodes; i++) {
        adjList[i] = std::unordered_set<uint>();

    }
    
}

WeightedVertexGraph::~WeightedVertexGraph(){
    //for (int i=0; i<numberOfNodes; i++) {
    //    delete this->adjList[i];
    //}
    delete [] this->adjList;
    
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
    for(int i = 0; i < numberOfEdges; ++i){
        if (!nodeSubset[edgesArray[i].first] && !nodeSubset[edgesArray[i].second] ) return false;
    }
    return true;
}