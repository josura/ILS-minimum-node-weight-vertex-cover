#pragma  once

#include <array>
#include <cstddef>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>
#include <bitset>
#include <random>


#include "WeightedVertexGraph.h"

//custom types
typedef std::vector<uint> NodeList;
using NodeBitList = std::vector<bool>;
using NodeBitArray = bool*;
using NodeSet = std::unordered_set<uint>;


std::ostream& operator<< (std::ostream &out, WeightedVertexGraph const& data);

std::ostream& operator<< (std::ostream &out, NodeBitList const& data);
std::ostream& operator<< (std::ostream &out, NodeList const& data);
std::ostream& operator<< (std::ostream &out, NodeSet const& data);

NodeList* nodeBitArrayToList(NodeBitArray const& nodeArray,uint arraySize);
NodeSet* nodeBitArrayToSet(NodeBitArray const& nodeArray,uint arraySize);

void printNodeBitArray(NodeBitArray nodeArray,uint size);

uint randomNumber(uint min, uint max);

NodeBitArray randomBooleanArray(uint size);