//
// Created by Zeyu Wang on 2021/9/24.
//

#ifndef DUMPY_FULLARYTREENODE_H
#define DUMPY_FULLARYTREENODE_H
#include <string>
#include <vector>
#include <fstream>
#include <unordered_map>
#include "../Const.h"
#include "../Utils/SaxUtil.h"
#include "TimeSeries.h"
#include "PqItemSeries.h"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/unordered_map.hpp>

using namespace std;

class FullAryTreeNode {

    FullAryTreeNode(int _id, int _layer, const unsigned short *_sax, FullAryTreeNode* parent): id(_id), layer(_layer){
        offsets.clear();
        for(int i=0;i<Const::segmentNum;++i)
            sax[i] = _sax[i];
        SaxUtil::id2Sax2(_id, sax, Const::segmentNum);
        if(_layer == 1)
            file_id = to_string(_id);
        else if (_layer == 2){
            file_id = parent->file_id;
        }else{
            if(!parent->flag){
                parent->file_id += to_string(parent->id);
                parent->flag = true;
            }
            file_id = parent->file_id;
        }
    }

    FullAryTreeNode(){for(unsigned short & i : sax)   i = 0; offsets.clear();}

    void insert(int offset);

public:
    vector<int>offsets{};
    int id = -1;
    int partitionId = -1;
    void exactSearchKnn(int k, TimeSeries *queryTs, vector<PqItemSeries*>&heap) const;
    FullAryTreeNode *routeToTarget(unsigned short *asax);
    static FullAryTreeNode *BuildFullAryTree();

    int size = 0;
    bool flag = false;
    unsigned short sax[Const::segmentNum]{};
    int layer = 0;
    string file_id;
    unordered_map<int, FullAryTreeNode*>*children{};
    static unsigned short * saxes;
    static float *dataset;
    static int a,b,c;

    friend class boost::serialization::access;
    template<class Archive>
            void serialize(Archive &ar, const unsigned int version){
                ar & id; ar & offsets;
                ar & size;
                ar & sax;
                ar & layer;
                ar & children;
                ar & partitionId;
            }

    void save2Disk(const string& output){
        ofstream ofs(output.c_str(), ios::binary);
        boost::archive::binary_oarchive oa(ofs);
        oa << (*this);
        ofs.close();
    }

    static FullAryTreeNode *loadFromDisk(const string &idxfn);

    int getLeafNodeNumber(FullAryTreeNode *root);

    int getHeight(FullAryTreeNode *root);

    bool isLeafNode() const{return children == nullptr;}

    int getGraphNum(FullAryTreeNode *root);

    FullAryTreeNode *route(unsigned short *asax) const;

    static bool order(const FullAryTreeNode*a, const FullAryTreeNode *b){
        return a->size > b->size;
    }

    static void mergingTardis(FullAryTreeNode *root);

    void exactSearchKnnInMemory(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap) const;

    void exactSearchKnnInMemory(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap, int &threshold) const;

    int getTotalSize(FullAryTreeNode *root);

    int getNodeNumber(FullAryTreeNode *root);

    void getIndexStats();

    int getSumHeight(FullAryTreeNode *root) const;

    static void generateSaxTbl();

    void getLeafNodeSize(FullAryTreeNode *root, ofstream &f);
};


#endif //DUMPY_FULLARYTREENODE_H
