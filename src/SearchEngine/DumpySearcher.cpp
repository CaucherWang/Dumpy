//
// Created by pengwang5 on 2022/1/16.
//

#include <set>
#include <unordered_set>
#include <chrono>
#include <cmath>
#include "../../include/Searchers/DumpySearcher.h"
#include "../../include/DataStructures/PqItemSeries.h"
#include "../../include/DataStructures/TimeSeries.h"
#include "../../include/Utils/SaxUtil.h"
#include "../../include/Utils/MathUtil.h"
#include "../../include/Const.h"

extern long  LOADED_NODE_CNT = 0;

static DumpyNode* targetNode;
vector<PqItemSeries *> * DumpySearcher::approxSearch(DumpyNode *root, float *query, int k, vector<vector<int>> *g,
                                                     const string &index_dir) {
    auto* queryTs = new TimeSeries(query);
    auto*heap = new vector<PqItemSeries*>();

    unsigned short sax[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        sax[i] = (*(queryTs->sax))[i];

    int head = SaxUtil::invSaxHeadFromSax(sax, Const::bitsCardinality, Const::segmentNum);
    DumpyNode *cur = (root->children)[head];
    if(cur == nullptr){
        // Since nodes in the first layer are usually numerous, we adopt static Hamming-distance based search for efficiency
        DumpyNode *node = nullptr;
        for(int i=0; i < DumpyNode::a + DumpyNode::b + DumpyNode::c; ++i){
            if(root->children[(*g)[head][i]] != nullptr){
                node = root->children[(*g)[head][i]];
                break;
            }
        }
        assert(node!=nullptr);
        // we only concern whether the nearest node is a leaf or an internal node
        if(!node->isLeafNode()){
            approxSearchInterNode(node, queryTs, sax, k, heap, index_dir);
        }else { node->search(k, queryTs, *heap, index_dir); targetNode = node;}
    }else if(cur->isLeafNode()){
        { cur->search(k, queryTs, *heap, index_dir); targetNode = cur;}
    }else approxSearchInterNode(cur, queryTs, sax, k, heap, index_dir);
    
    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;
}

void DumpySearcher::approxSearchInterNode(DumpyNode *root, TimeSeries *queryTs, unsigned short *sax, int k,
                                          vector<PqItemSeries *> *heap, const string &index_dir) {
    DumpyNode *cur = root->route(sax);
    if(cur->isLeafNode()){
        cur->search(k, queryTs, *heap, index_dir);
        targetNode = cur;
        return;
    }

    double min_dist = numeric_limits<double>::max(), max_size = 0;
    DumpyNode *node;
    for(int i=0;i<cur->children.size();++i){
        if(cur->children[i] == nullptr)  continue;
        double dist;
        if(cur->children[i]->isLeafNode())  dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, cur->sax, cur->bits_cardinality, cur->children[i]->chosenSegments, i);
        else dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, cur->children[i]->sax, cur->children[i]->bits_cardinality);
        if(dist < min_dist){
            min_dist = dist;
            max_size = cur->children[i]->size;
            node = cur->children[i];
        }else if(dist == min_dist && cur->children[i]->size > max_size){
            max_size = cur->children[i]->size;
            node = cur->children[i];
        }
    }

    // we only concern whether the nearest node is a leaf or an internal node
    if(!node->isLeafNode()){
        approxSearchInterNode(node, queryTs, sax, k, heap, index_dir);
        return;
    }else { node->search(k, queryTs, *heap, index_dir); targetNode = node;}
}


struct PqItemDumpy{
    DumpyNode* node{};
    double dist{};

    PqItemDumpy(DumpyNode* n, double d){ node = n; dist = d;}
    PqItemDumpy(){ node = nullptr; dist = 0;}

    bool operator <(const PqItemDumpy & pt) const{
        if(node == pt.node) return false;
        else {
            if(dist != pt.dist) return dist < pt.dist;
            else return node < pt.node;
        }
    }
};

vector<PqItemSeries*>*DumpySearcher::exactSearch(DumpyNode* root, float *query, int k, vector<vector<int>> *g){

    vector<PqItemSeries*>* heap = approxSearch(root, query, k, g, Const::idxfn);
    make_heap(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    double bsf = (*heap)[0]->dist;
    auto *queryTs = new TimeSeries(query);

    set<PqItemDumpy>pq;
    pq.insert(PqItemDumpy(root, 0));

    PqItemDumpy cur;
    while(!pq.empty()){
        cur = *pq.begin();
        if(cur.dist > bsf)  break;
        pq.erase(pq.begin());
        if(!cur.node->isLeafNode()){
            unordered_set<DumpyNode*>inserted;
            for(DumpyNode* node:cur.node->children)
                if(node != nullptr && node != targetNode && inserted.find(node) == inserted.end()) {
                    inserted.insert(node);
                    double lb_dist = SaxUtil::LowerBound_Paa_iSax(queryTs->paa, node->sax, node->bits_cardinality);
                    if(lb_dist < bsf){
                        pq.insert(PqItemDumpy(node, lb_dist));
                    }
                }
            inserted.clear();
        }else{
            cur.node->search(k, queryTs, *heap, Const::idxfn);
            ++LOADED_NODE_CNT;
            bsf = (*heap)[0]->dist;
        }
    }
    pq.clear();
    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;
}

static float *t_paa;
bool comp_Dumpy(const DumpyNode* x, const DumpyNode* y){
    if(x == nullptr)    return false;
    if(y == nullptr)    return true;
    return SaxUtil::LowerBound_Paa_iSax(t_paa, x->sax, x->layer) < SaxUtil::LowerBound_Paa_iSax(t_paa, y->sax, y->layer);
}

vector<PqItemSeries *> * DumpySearcher::approxIncSearch(DumpyNode *root, float *query, int k, const string &index_dir,
                                                        int node_num) {
    auto* queryTs = new TimeSeries(query);
    t_paa = queryTs->paa;
    auto*heap = new vector<PqItemSeries*>();

    unsigned short sax[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        sax[i] = (*(queryTs->sax))[i];

    approxIncSearchInterNode(root, queryTs, sax, k, heap, index_dir, node_num);

    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;

}

void DumpySearcher::approxIncSearchInterNode(DumpyNode *root, TimeSeries *queryTs, unsigned short *sax, int k,
                                             vector<PqItemSeries *> *heap, const string &index_dir, int &node_num) {
    if(root->isLeafNode() || node_num <= 0)  return;
    DumpyNode *cur = root->route1step(sax), *parent = root;
    while (cur!= nullptr && !cur->isLeafNode() && cur->getLeafNodeNum() > node_num) {
        parent = cur;
        cur = cur->route(sax);
    }

    if(cur!= nullptr){
        if(cur->isLeafNode()){
            cur->search(k, queryTs, *heap, index_dir);
            --node_num;
        }else{
            approxIncSearchInterNode(cur, queryTs, sax, k, heap, index_dir, node_num);
        }
    }

    if(node_num <=0)    return;

    vector<DumpyNode*>candidates;
    unordered_set<DumpyNode*>cands;
    for(DumpyNode *node: parent->children)
        if(node != nullptr && node!=cur && cands.find(node) == cands.end()) {
            candidates.push_back(node);
            cands.insert(node);
        }
    cands.clear();
    sort(candidates.begin(), candidates.end(), comp_Dumpy);

    for(int i=0;i<candidates.size() && node_num > 0;++i){
        if(candidates[i]->isLeafNode()) {
            candidates[i]->search(k, queryTs, *heap, index_dir);
            --node_num;
        }
        else {
            approxIncSearchInterNode(candidates[i], queryTs, sax, k, heap, index_dir, node_num);
        }
    }

}


vector<PqItemSeries *> * DumpySearcher::approxIncSearchFuzzy(DumpyNode *root, float *query, int k, const string &index_dir,
                                                             int node_num) {
    auto* queryTs = new TimeSeries(query);
    t_paa = queryTs->paa;
    auto*heap = new vector<PqItemSeries*>();
    auto*hash_set = new unordered_set<float*, createhash, isEqual>();

    unsigned short sax[Const::segmentNum];
    for(int i=0;i<Const::segmentNum;++i)
        sax[i] = (*(queryTs->sax))[i];

    approxIncSearchInterNodeFuzzy(root, queryTs, sax, k, heap, index_dir, node_num, hash_set);

    delete queryTs;
    sort(heap->begin(), heap->end(), PqItemSeriesMaxHeap());
    return heap;

}


void DumpySearcher::approxIncSearchInterNodeFuzzy(DumpyNode *root, TimeSeries *queryTs, unsigned short *sax, int k,
                                                  vector<PqItemSeries *> *heap, const string &index_dir, int &node_num,
                                                  unordered_set<float*, createhash, isEqual>*hash_set) {
    if(root->isLeafNode() || node_num <= 0)  return;
    DumpyNode *cur = root->route1step(sax), *parent = root;
    while (cur!= nullptr && !cur->isLeafNode() && cur->getLeafNodeNum() > node_num) {
        parent = cur;
        cur = cur->route(sax);
    }

    if(cur!= nullptr){
        if(cur->isLeafNode()){
            cur->search(k, queryTs, *heap, index_dir);
            --node_num;
        }else{
            approxIncSearchInterNode(cur, queryTs, sax, k, heap, index_dir, node_num);
        }
    }

    if(node_num <=0)    return;

    vector<DumpyNode*>candidates;
    unordered_set<DumpyNode*>cands;
    for(DumpyNode *node: parent->children)
        if(node != nullptr && cands.find(node) == cands.end()) {
            candidates.push_back(node);
            cands.insert(node);
        }
    cands.clear();

    sort(candidates.begin(), candidates.end(), comp_Dumpy);

    for(int i=0;i<candidates.size() && node_num > 0;++i){
        if(candidates[i]->isLeafNode()) {
            candidates[i]->search(k, queryTs, *heap, index_dir, hash_set);
            --node_num;
        }
        else {
            approxIncSearchInterNodeFuzzy(candidates[i], queryTs, sax, k, heap, index_dir, node_num, hash_set);
        }
    }

}