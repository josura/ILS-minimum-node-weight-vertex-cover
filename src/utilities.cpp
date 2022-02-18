
#include "utilities.h"
#include <cstddef>
#include <iostream>
#include <random>
#include <unordered_set>
#include <vector>

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

NodeSet* nodeBitArrayToSet(NodeBitArray const& nodeArray, uint arraySize){
    NodeSet* ret = new NodeSet;
    for(uint i = 0 ; i < arraySize ; i++){
        if(nodeArray[i]){
            ret->insert(i);
        }
    }
    return ret;
}

void printNodeBitArray(NodeBitArray nodeArray,uint size){
    for (uint i = 0; i < size; i++) {
        if(nodeArray[i]){
            cout << "1";
        } else {
            cout << "0";
        }
    }
    cout<<endl;
}


int randomNumber(int min, int max){
    // fixed seed because repeatability
    //unsigned seed = 777;
    std::random_device r;
    // range [min,max[
    std::default_random_engine e1(r());
    max = (max-1<0) ? 0 : max-1;
    std::uniform_int_distribution<int> uniform_dist(min, max);
    return uniform_dist(e1);
}

double randomRealNumber(double min, double max){
    // fixed seed because repeatability
    //unsigned seed = 777;
    std::random_device r;
    // range [min,max]
    std::default_random_engine e1(r());
    std::uniform_real_distribution<double> uniform_dist(min, max);
    return uniform_dist(e1);
}


std::vector<uint> randomVector(int min, int max , uint size){
    unsigned seed = std::chrono::system_clock::now()
                        .time_since_epoch()
                        .count();
    // fixed seed because repeatability
    //unsigned seed = 777;
    std::vector<uint> v(max-min,min);
    for (uint i = 0; i<max-min; i++) {
        v[i] += i;
    }
 
    std::shuffle(std::begin(v), std::end(v), std::default_random_engine(seed));
    size = (size<v.size()) ? size : v.size();
    std::vector<uint> ret(v.begin(),v.begin()+size);
    return ret;
}

NodeBitArray randomBooleanArray(uint size){
    NodeBitArray ret = new bool[size];
    for (uint i = 0; i< size; i++) {
        ret[i] = (randomNumber(0, 100) >=50) ? false : true;
    }

    return ret;
}


void printUsage(std::string execName){
    std::cout << "Usage : "<< execName << " <inputgraph>"<<std::endl;
}


std::string nodeBitArrayToString(NodeBitArray nodeArray,uint size){
    std::string ret = "";
    for (uint i = 0; i < size; i++) {
        if(nodeArray[i]){
             ret += "1";
        } else {
             ret += "0";
        }
    }
    return ret;
}


std::unordered_set<uint> intersectionSet(std::unordered_set<uint> set1,std::unordered_set<uint> set2){
    unordered_set<uint> ret;
    for (auto it = set1.begin(); it != set1.end(); it++) {
        if(set2.count(*it)>0){
            ret.insert(*it);
        }
    }
    return ret;
}