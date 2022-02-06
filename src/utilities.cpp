
#include "utilities.h"
#include <cstddef>
#include <iostream>

using namespace std;

std::ostream& operator<< (std::ostream &out, WeightedVertexGraph const& data) {
            out << data.getNumNodes() << " " << data.getNumEdges() <<std::endl;
            string nodeweights = data.getNodeWeightsStr();
            out << nodeweights << std::endl;
            out << "Adj Lists" << std::endl;
            for(int i = 0; i<data.getNumNodes(); i++){
                out << "node " << i << " :" << data.getAdjListStr(i) << std::endl;
            }
            return out;
        }


std::ostream& operator<< (std::ostream &out, NodeBitList const& data) {
            for(auto it = data.begin() ; it!=data.end(); it++){
                if(*it){
                    out << " true ";
                } else{
                    out << " false ";
                }
            }
            return out;
        }


std::ostream& operator<< (std::ostream &out, NodeList const& data) {
            out << "(";
            for(auto it = data.begin() ; it!=data.end(); it++){
                out << *it << ",";
            }
            out << ")";
            return out;
        }


std::ostream& operator<< (std::ostream &out, NodeSet const& data) {
            out << "(";
            for(auto it = data.begin() ; it!=data.end(); it++){
                out << *it << ",";
            }
            out << ")";
            return out;
        }


NodeList* nodeBitArrayToList(NodeBitArray const& nodeArray, uint arraySize){
    NodeList* ret = new NodeList;
    for(uint i = 0 ; i < arraySize ; i++){
        if(nodeArray[i]){
            ret->push_back(i);
        }
    }
    return ret;
}

void printNodeBitArray(NodeBitArray nodeArray,uint size){
    for (uint i = 0; i < size; i++) {
        if(nodeArray[i]){
            cout << " true ";
        } else {
            cout << "false";
        }
    }
}