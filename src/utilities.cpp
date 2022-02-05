
#include "utilities.h"

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