#include <iostream>
#include <cstdio>
#include <vector>
#include <thread>
#include <chrono>
#include "../include/DataStructures/DumpyNode.h"
#include "../include/DataStructures/GraphConstruction.h"
#include "../include/Searchers/DumpySearcher.h"
#include "../include/Searchers/DumpyInMemorySearch.h"
#include "../include/DataStructures/FullAryTreeNode.h"
#include "../include/Utils/FileUtil.h"
#include "../include/Utils/MathUtil.h"
#include "../include/Utils/TimeSeriesUtil.h"
using namespace std;

vector<vector<int>>* loadGraphSkeleton(){
    int vd = 0;
    for(int i=1; i<=Const::bitsReserve;++i)
        vd += MathUtil::nChooseK(Const::segmentNum, i);
    auto nnList = new vector<vector<int>>(Const::vertexNum, vector<int>(vd, -1));

    if(!FileUtil::checkFileExists(Const::graphfn.c_str())){
        cout << "File not exists!" << Const::graphfn << endl;
        exit(-1);
    }
    FILE *f = fopen(Const::graphfn.c_str(), "rb");

    for(int i=1;i<Const::vertexNum;++i)
        fread(&((*nnList)[i][0]), sizeof(int), vd, f);

    return nnList;

}

void constructGraph(){
    GraphConstruction::buildAndSave2Disk();
}

void buildInMemoryIndexDumpy(){
    FullAryTreeNode* root = FullAryTreeNode::BuildFullAryTree();
    cout << "build finish" << endl;
    root->getIndexStats();
    root->save2Disk(Const::memoryidxfn);
}

void searchInMemory(){
    FullAryTreeNode* root = FullAryTreeNode::loadFromDisk(Const::memoryidxfn);
    cout << "load finish." << endl;
    float *queries = FileUtil::readQueries();
    cout << fixed <<setprecision(6);
    int search_num  = Const::series_num * Const::search_ratio;
    for(int i=0;i<Const::query_num;++i){
        Const::logPrint("Query " + to_string(i) +":");
        vector<PqItemSeries*> *approxKnn = DumpyInMemorySearch::approxSearch(root, queries + i * Const::tsLength, Const::k, search_num);
        Const::logPrint("Results:");
        for (int j = 0; j < approxKnn->size(); ++j) {
            cout << j + 1 << ": " << TimeSeriesUtil::timeSeries2Line((*approxKnn)[j]->ts) << endl;
        }
    }
}

void statMemoryDumpy(){
    FullAryTreeNode* root = FullAryTreeNode::loadFromDisk(Const::memoryidxfn);
    cout << "load finish" << endl;
    root->getIndexStats();
}

void buildDumpy(){
    auto g = loadGraphSkeleton();
    DumpyNode* root = DumpyNode::BuildIndex(Const::datafn, Const::saxfn, Const::paafn, g);
    root->save2Disk(Const::idxfn + "root.idx");
}

void approxSearchOneNode() {
    DumpyNode *root = DumpyNode::loadFromDisk(Const::saxfn, Const::idxfn + "root.idx", false);
    auto *g = loadGraphSkeleton();
    float *queries = FileUtil::readQueries();
    for (int i = 0; i < Const::query_num; ++i) {
        Const::logPrint("Query " + to_string(i) + ":");
        vector<PqItemSeries *> *approxKnn = DumpySearcher::approxSearch(root, queries + i * Const::tsLength, Const::k,
                                                                        g, Const::idxfn);
        Const::logPrint("Results:");
        for (int j = 0; j < approxKnn->size(); ++j) {
            cout << j + 1 << ": " << TimeSeriesUtil::timeSeries2Line((*approxKnn)[j]->ts) << endl;
        }
    }
}
void approxSearchMoreNode() {
    DumpyNode *root = DumpyNode::loadFromDisk(Const::saxfn, Const::idxfn + "root.idx", false);
    float *queries = FileUtil::readQueries();
    for (int i = 0; i < Const::query_num; ++i) {
        Const::logPrint("Query " + to_string(i) + ":");
        vector<PqItemSeries *> *approxKnn = DumpySearcher::approxIncSearch(root, queries + i * Const::tsLength,
                                                                           Const::k, Const::idxfn,
                                                                           Const::visited_node_num);
        Const::logPrint("Results:");
        for (int j = 0; j < approxKnn->size(); ++j)
            cout << j + 1 << ": " << TimeSeriesUtil::timeSeries2Line((*approxKnn)[j]->ts) << endl;
    }
}

void buildDumpyFuzzy(){
    auto g = loadGraphSkeleton();
    int bound = Const::fuzzy_f * 100;
    Const::fuzzyidxfn += "/" + to_string(bound) + "-" + to_string(Const::delta) + "/";
    DumpyNode* root = DumpyNode::BuildIndexFuzzy(Const::datafn, Const::saxfn, Const::paafn, g);
    root->save2Disk(Const::fuzzyidxfn + "root.idx");
}

void approxSearchOneNodeFuzzy() {
    int bound = Const::fuzzy_f * 100;
    Const::fuzzyidxfn += "/" + to_string(bound) + "-" + to_string(Const::delta) + "/";
    DumpyNode *root = DumpyNode::loadFromDisk(Const::saxfn, Const::fuzzyidxfn + "root.idx", false);
    auto *g = loadGraphSkeleton();
    float *queries = FileUtil::readQueries();
    for (int i = 0; i < Const::query_num; ++i) {
        Const::logPrint("Query " + to_string(i) + ":");
        vector<PqItemSeries *> *approxKnn = DumpySearcher::approxSearch(root, queries + i * Const::tsLength, Const::k,
                                                                        g, Const::fuzzyidxfn);
        Const::logPrint("Results:");
        for (int j = 0; j < approxKnn->size(); ++j)
            cout << j + 1 << ": " << TimeSeriesUtil::timeSeries2Line((*approxKnn)[j]->ts) << endl;
    }
}
void approxSearchMoreNodeFuzzy() {
    int bound = Const::fuzzy_f * 100;
    Const::fuzzyidxfn += "/" + to_string(bound) + "-" + to_string(Const::delta) + "/";
    DumpyNode *root = DumpyNode::loadFromDisk(Const::saxfn, Const::fuzzyidxfn + "root.idx", false);
    float *queries = FileUtil::readQueries();
    for (int i = 0; i < Const::query_num; ++i) {
        Const::logPrint("Query " + to_string(i) + ":");
        auto start = chrono::system_clock::now();
        vector<PqItemSeries *> *approxKnn = DumpySearcher::approxIncSearchFuzzy(root, queries + i * Const::tsLength,
                                                                                Const::k, Const::fuzzyidxfn,
                                                                                Const::visited_node_num);
        Const::logPrint("Results:");
        for (int j = 0; j < approxKnn->size(); ++j)
            cout << j + 1 << ": " << TimeSeriesUtil::timeSeries2Line((*approxKnn)[j]->ts) << endl;
    }
}

void exactExprDumpy() {
    DumpyNode *root = DumpyNode::loadFromDisk(Const::saxfn, Const::idxfn + "root.idx", false);
    auto *g = loadGraphSkeleton();
    float *queries = FileUtil::readQueries();
    for (int i = 0; i < Const::query_num; ++i) {
        Const::logPrint("Query " + to_string(i) + ":");
        auto start = chrono::system_clock::now();
        vector<PqItemSeries *> *exactKnn = DumpySearcher::exactSearch(root, queries + i * Const::tsLength, Const::k, g);
        Const::logPrint("Results:");
        for (int j = 0; j < exactKnn->size(); ++j)
            cout << j + 1 << ": " << TimeSeriesUtil::timeSeries2Line((*exactKnn)[j]->ts) << endl;
    }
}

void statIndexDumpy(){
    DumpyNode* root = DumpyNode::loadFromDisk(Const::saxfn, Const::idxfn + "root.idx", false);
    root->getIndexStats();
}

void statIndexDumpyFuzzy(){
    int bound = Const::fuzzy_f * 100;
    Const::fuzzyidxfn += "/" + to_string(bound) + "-" + to_string(Const::delta) + "/";
    DumpyNode* root = DumpyNode::loadFromDisk(Const::saxfn, Const::fuzzyidxfn + "root.idx", false);
    root->getIndexStats();
}

int main() {
    Const::readConfig();

    switch (Const::index) {
        case 0:
            constructGraph();
            break;
        case 1:
            switch (Const::ops) {
                case 0:
                    buildDumpy();
                    break;
                case 1:
                    approxSearchOneNode();
                    break;
                case 2:
                    exactExprDumpy();
                    break;
                case 3:
                    statIndexDumpy();
                    break;
                case 4:
                    approxSearchMoreNode();
                default:
                    break;
            }
            break;
        case 2:
            if(Const::ops == 0){
                buildDumpyFuzzy();
            }
            else if(Const::ops == 1){
                approxSearchOneNodeFuzzy();
            }else if(Const::ops == 3){
                statIndexDumpyFuzzy();
            }else if(Const::ops == 4){
                approxSearchMoreNodeFuzzy();
            }
            break;
        case 3:
            if(Const::ops == 0) buildInMemoryIndexDumpy();
            else if(Const::ops == 1) searchInMemory();
            else if(Const::ops == 4)    statMemoryDumpy();
            break;
        default:    break;
    }
}
