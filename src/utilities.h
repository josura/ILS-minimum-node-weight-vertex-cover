#pragma  once

#include <cstddef>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <bitset>


#include "WeightedVertexGraph.h"

//custom types
typedef std::vector<uint> NodeList;
using NodeBitList = std::vector<bool>;


std::ostream& operator<< (std::ostream &out, WeightedVertexGraph const& data);

std::ostream& operator<< (std::ostream &out, NodeBitList const& data);
