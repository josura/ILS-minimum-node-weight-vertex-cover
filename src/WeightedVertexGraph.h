#pragma once

#include <ostream>
#include <string>
#include <sys/types.h>
#include <unordered_set>
#include <iostream>
#include <vector>
#include <utility>

class WeightedVertexGraph{
    private:
        uint numberOfNodes;
        uint numberOfEdges=0;
        double* nodeWeights;
        std::unordered_set<uint>* adjList;
        std::vector<std::pair<uint, uint> > edgesVector;
        std::pair<uint, uint>* edgesArray;

    public:
        WeightedVertexGraph();

        WeightedVertexGraph(uint numNodes, double* nodeWeights);

        ~WeightedVertexGraph();

        WeightedVertexGraph* addEdge(uint node1, uint node2);

        std::pair<uint, uint>* makeEdgesArray();

        // accessory functions

        uint getNumNodes()const ;
        uint getNumEdges()const ;

        double* getNodeWeights()const;

        std::string getNodeWeightsStr()const;

        std::unordered_set<uint>* getAdjList(uint node)const;

        std::string getAdjListStr(uint node)const;

        bool adjNodes(uint node1, uint node2);

        std::vector<std::pair<uint, uint>> getEdgesVector()const;

        std::pair<uint, uint>* getEdgesArray()const;

        // optimization methods

        bool vertexCoverValidityEdgescheckBitArray(bool* nodeSubset);

};