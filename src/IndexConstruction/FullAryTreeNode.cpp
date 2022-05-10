//
// Created by Zeyu Wang on 2021/9/24.
//

#include "../../include/DataStructures/FullAryTreeNode.h"
#include "../../include/Utils/FileUtil.h"
#include "../../include/Utils/TimeSeriesUtil.h"
#include "../../include/Const.h"
#include "../../include/Utils/MathUtil.h"
#include "../../include/DataStructures/DumpyNode.h"
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <chrono>
#include <cmath>
#include <iostream>

unsigned short * FullAryTreeNode::saxes = nullptr;
float *FullAryTreeNode::dataset = nullptr;
int FullAryTreeNode::a = MathUtil::nChooseK(Const::segmentNum, 1), FullAryTreeNode::b = MathUtil::nChooseK(Const::segmentNum, 2), FullAryTreeNode::c = MathUtil::nChooseK(Const::segmentNum, 3);

static long MAT1_TOTAL_TIME = 0, MAT1_READ_TIME = 0, MAT2_WRITE_TIME = 0, MAT2_TOTAL_TIME = 0, SAX_PAA_TOTAL_TIME = 0,
        MAT2_READ_TIME = 0, GROW_CPU_TIME = 0,  GROW_CPU_TIME_1st = 0, GROW_TOTAL_TIME = 0,
        RAND_READ_CNT = 0, RAND_WRITE_CNT = 0, SEQ_READ_CNT = 0, SEQ_WRITE_CNT = 0, SAX_PAA_READ_TIME, SAX_PAA_CPU_TIME;


FullAryTreeNode * FullAryTreeNode::BuildFullAryTree() {
    generateSaxTbl();

    auto* root = new FullAryTreeNode();
    root->children = new unordered_map<int, FullAryTreeNode*>();

    for(int i=0;i<Const::series_num;++i) {
        if(i%10000000 == 0)  cout << i << endl;
        root->insert(i);
    }

    return root;
}

void FullAryTreeNode::insert(int offset) {
    unsigned short *asax = saxes + (long)offset * Const::segmentNum;
    ++size;
    if(!isLeafNode()) {
        FullAryTreeNode *target = routeToTarget(asax);
        if(target->layer == 1 && !target->isLeafNode())
            target->offsets.push_back(offset);
        target->insert(offset);
    } else if(size <= Const::th || layer >= Const::bitsCardinality)   offsets.push_back(offset);
    else{   //split
        offsets.push_back(offset);
        children = new unordered_map<int, FullAryTreeNode*>();
        for(int off:offsets){
            unsigned short* _sax = saxes + (long)off * Const::segmentNum;;
            FullAryTreeNode* target = routeToTarget(_sax);
            target->insert(off);
        }
        if(layer > 1)
            vector<int>().swap(offsets);
    }
}

FullAryTreeNode* FullAryTreeNode::routeToTarget(unsigned short *asax){
    assert(!isLeafNode());
    int nav_id = SaxUtil::invSaxHeadkFromSax(asax, Const::bitsCardinality, Const::segmentNum, layer + 1);
    if(children->find(nav_id) == children->end()) {
        auto *target =  new FullAryTreeNode(nav_id, layer + 1, sax, this);
        (*children)[nav_id] = target;
    }
    return (*children)[nav_id];
}

FullAryTreeNode* FullAryTreeNode::route(unsigned short *asax) const{
    assert(!isLeafNode());
    int nav_id = SaxUtil::invSaxHeadkFromSax(asax, Const::bitsCardinality, Const::segmentNum, layer + 1);
    return (*children)[nav_id];
}

void FullAryTreeNode::exactSearchKnnInMemory(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap) const {
    if(isLeafNode()){
        double bsf = heap.size() < k? numeric_limits<double>::max() : heap[0]->dist;
        float *ts;
        for(int offset: offsets){
            ts = dataset + (long)offset * Const::tsLength;
            double dist = TimeSeriesUtil::euclideanDist(queryTs->ts, ts, Const::tsLength);
            if(heap.size() < k){
                heap.push_back(new PqItemSeries(ts, dist));
                push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            }else if(dist < bsf){
                pop_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
                delete heap.back();
                heap.pop_back();
                heap.push_back(new PqItemSeries(ts, dist));
                push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            }
            if(heap.size() >= k)    bsf = heap[0]->dist;
        }
    }
    else{
        for(auto & iter : *children)
            if(iter.second != nullptr)
                iter.second->exactSearchKnnInMemory(k, queryTs, heap);
    }
}

void FullAryTreeNode::exactSearchKnnInMemory(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap, int &threshold) const {
    if(isLeafNode()){
        double bsf = heap.size() < k? numeric_limits<double>::max() : heap[0]->dist;
        float *ts;
        int bound = min(threshold, size);
        for(int i=0;i<bound;++i){
            ts = dataset + (long)offsets[i] * Const::tsLength;
            double dist = TimeSeriesUtil::euclideanDist(queryTs->ts, ts, Const::tsLength);
            if(heap.size() < k){
                heap.push_back(new PqItemSeries(ts, dist));
                push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            }else if(dist < bsf){
                pop_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
                delete heap.back();
                heap.pop_back();
                heap.push_back(new PqItemSeries(ts, dist));
                push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            }
            if(heap.size() >= k)    bsf = heap[0]->dist;
        }
        threshold -=bound;
    }else{
        for(auto & iter : *children){
            if(threshold <=0) break;
            if(iter.second != nullptr)
                iter.second->exactSearchKnnInMemory(k, queryTs, heap, threshold);
        }
    }

}

void FullAryTreeNode::exactSearchKnn(int k, TimeSeries *queryTs, vector<PqItemSeries *> &heap) const {
    if(isLeafNode()){
        double bsf = heap.size() < k? numeric_limits<double>::max() : heap[0]->dist;
        FILE *f = fopen(Const::datafn.c_str(), "rb");
        float *ts;
        for(int offset: offsets){
            fseek(f, (long)offset * Const::tsLengthBytes, SEEK_SET);
            ts = new float[Const::tsLength];
            fread(ts, sizeof(float), Const::tsLength, f);
            double dist = TimeSeriesUtil::euclideanDist(queryTs->ts, ts, Const::tsLength);
            if(heap.size() < k){
                heap.push_back(new PqItemSeries(ts, dist, true));
                push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            }else if(dist < bsf){
                pop_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
                delete heap.back();
                heap.pop_back();
                heap.push_back(new PqItemSeries(ts, dist, true));
                push_heap(heap.begin(),  heap.end(), PqItemSeriesMaxHeap());
            } else delete[] ts;
            if(heap.size() >= k)    bsf = heap[0]->dist;
        }
        fclose(f);
    }
    else{
        for(auto & iter : *children)
            if(iter.second != nullptr)
                iter.second->exactSearchKnn(k, queryTs, heap);
    }
}

FullAryTreeNode *FullAryTreeNode::loadFromDisk(const string &idxfn) {
    long num = (long)Const::series_num * Const::tsLength;
    dataset = new float[num];
    FILE *f = fopen(Const::datafn.c_str(), "rb");
    fread(dataset, sizeof(float), num, f);
    fclose(f);
    ifstream ifs(idxfn.c_str(), ios::binary);
    boost::archive::binary_iarchive ia(ifs);
    auto *g = new FullAryTreeNode();
    ia >> (*g);
    ifs.close();
    return g;

}

void FullAryTreeNode::mergingTardis(FullAryTreeNode* root){
    if(root == nullptr || root->children == nullptr || root->size <= Const::th)   return;

    // do merge
    vector<FullAryTreeNode*>candidates;
    for(auto &iter:*root->children){
        FullAryTreeNode *tmp = iter.second;
        if(tmp != nullptr && tmp->isLeafNode() && tmp->size <= Const::th)
            candidates.push_back(tmp);
    }
    sort(candidates.begin(), candidates.end(), FullAryTreeNode::order);
    int pid = 0, cur_size = 0;
    for(FullAryTreeNode*_:candidates){
        if(cur_size + _->size <= Const::th) {
            _->partitionId = pid;
            cur_size += _->size;
        }else{
            ++pid;
            cur_size = _->size;
            _->partitionId = pid;
        }
    }

    // recursive
    for(auto &iter:*root->children){
        FullAryTreeNode *tmp = iter.second;
        if(tmp!= nullptr && !tmp->isLeafNode())
            mergingTardis(tmp);
    }
}

int FullAryTreeNode::getLeafNodeNumber(FullAryTreeNode * root){
    if(root == nullptr) return 0;
    if(root->size <= Const::th) return 1;
    int sum = 0;
    unordered_set<int>seen;
    for(auto& iter:*root->children){
        if(iter.second == nullptr) continue;
        if(iter.second->size <= Const::th){
            if(seen.find(iter.second->partitionId) != seen.end())   continue;
            sum++;
            seen.insert(iter.second->partitionId);
        }else{
            sum += getLeafNodeNumber(iter.second);
        }
    }
    return sum;
}

void FullAryTreeNode::getLeafNodeSize(FullAryTreeNode * root, ofstream &f){
    if(root == nullptr) return;
    if(root->size <= Const::th){
        int s = root->size;
        f << to_string(s) + ",";
        return;
    }
    int sum = 0;
    for(auto& iter:*root->children){
        if(iter.second == nullptr) continue;
        if(iter.second->size <= Const::th){
            f << to_string(iter.second->size) + ",";
        }else{
            getLeafNodeSize(iter.second, f);
        }
    }
    return;
}

int FullAryTreeNode::getNodeNumber(FullAryTreeNode * root){
    if(root == nullptr) return 0;
    if(root->children == nullptr) return 1;
    int sum = 0;
    unordered_set<int>seen;
    for(auto& iter:*root->children){
        if(iter.second == nullptr) continue;
        if(iter.second->size <= Const::th){
            if(seen.find(iter.second->partitionId) != seen.end())   continue;
            sum++;
            seen.insert(iter.second->partitionId);
        }else{
            sum += getNodeNumber(iter.second);
        }
    }
    return sum + 1;
}

int FullAryTreeNode::getTotalSize(FullAryTreeNode * root){
    if(root == nullptr) return 0;
    if(root->children == nullptr) return root->size;
    int sum = 0;
    unordered_set<int>seen;
    for(auto& iter:*root->children){
        if(iter.second == nullptr) continue;   // may not happen
        sum += getTotalSize(iter.second);
    }
    return sum;
}

int FullAryTreeNode::getHeight(FullAryTreeNode *root){
    if(root == nullptr) return 0;
    if(root->children == nullptr) return 1;
    int max = 0;
    for(auto& iter:*root->children){
        if(iter.second == nullptr) continue;
        int h = getHeight(iter.second);
        max = max >= h? max : h;
    }
    return max + 1;
}

int FullAryTreeNode::getSumHeight(FullAryTreeNode * root) const{
    if(root == nullptr) return 0;
    if(root->children == nullptr) return root->layer;
    int sum = 0;
    unordered_set<int>seen;
    for(auto& iter:*root->children){
        if(iter.second == nullptr || seen.find(iter.second->partitionId) != seen.end()) continue;   // may not happen
        sum += getSumHeight(iter.second);
        seen.insert(iter.second->partitionId);
    }
    return sum;
}

int FullAryTreeNode::getGraphNum(FullAryTreeNode* root){
    if(root == nullptr || root->children == nullptr || root->size <= Const::th)   return 0;
    int sum = 1;
    for(auto iter:*root->children){
        sum += getGraphNum(iter.second);
    }
    return sum;
}

void FullAryTreeNode::getIndexStats(){
    ofstream outfile;
    outfile.open("/mnt/c/codes/rand-leaf.out", ios::out | ios::app);//输入文件的路径
    getLeafNodeSize(this, outfile);
    int total_leaf_node_num = getLeafNodeNumber(this);
    int total_size = getTotalSize(this);
    cout << "Total size = " << total_size << endl;
    cout <<"Total nodes number = " << getNodeNumber(this) << endl;
    cout << "Leaf node number = " << total_leaf_node_num << endl;
    cout << "Max. height = " << getHeight(this) - 1 <<endl;
    cout << "Avg. Height = " << (double)getSumHeight(this) / (double) total_leaf_node_num << endl;
    cout <<"Avg. Filling Factor = "<< total_size / (double)total_leaf_node_num / Const::th << endl;
    outfile.close();
}

void FullAryTreeNode::generateSaxTbl(){
    string fn = Const::datafn;
    long fs = FileUtil::getFileSize(fn.c_str());
    long series_num = fs / Const::tsLengthBytes;
    if(series_num < Const::segmentNum) {
        cout << " File does not have enough data series!"<<endl;
        exit(-1);
    }
    series_num = Const::series_num;
    float ts[Const::tsLength];
    saxes = new unsigned short[series_num * Const::segmentNum];
    long rest = series_num, cur = 0;
    FILE *f = fopen(fn.c_str(), "rb");

    RAND_READ_CNT++;
    SEQ_READ_CNT += series_num;
    RAND_WRITE_CNT+=2;
    SEQ_WRITE_CNT += series_num;

    while(rest > 0){
        long num;
        if(rest > 4000000)    num = 4000000;
        else num = rest;
        auto *tss = new float[num * Const::tsLength];

        auto start = chrono::system_clock::now();
        fread(tss, sizeof(float),num * Const::tsLength,  f);
        auto end = chrono::system_clock::now();
        SAX_PAA_READ_TIME += chrono::duration_cast<chrono::microseconds>(end - start).count();

        for(long i=0;i<num;++i){
            if(isnan(tss[i * Const::tsLength])){
                for(int j = 0; j < Const::segmentNum; ++j)
                    saxes[i * Const::segmentNum + j] = 0;
                cout << "Dirty data: "<<i << "," <<endl;
            }
            else{
                SaxUtil::saxFromTs(tss + i * Const::tsLength, saxes + (cur+ i) * Const::segmentNum,
                                   Const::tsLengthPerSegment, Const::segmentNum, Const::cardinality);
            }
        }
        delete[] tss;
        rest -=num;
        cur+=num;
        start = chrono::system_clock::now();
        SAX_PAA_CPU_TIME += chrono::duration_cast<chrono::microseconds>(start - end).count();
    }

    fclose(f);
}