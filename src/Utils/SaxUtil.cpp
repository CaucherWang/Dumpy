//
// Created by Zeyu Wang on 2021/8/7.
//

#include "../../include/Utils/SaxUtil.h"
#include "../../include/Utils/FileUtil.h"
#include "../../include/DataStructures/TimeSeries.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cassert>
using namespace std;


double * SaxUtil::breakpoints = readDoubleFromFileAtOnce();

std::vector<std::string> split(std::string& str, std::string pattern)
{
    std::string::size_type pos;
    std::vector<std::string> result;
    str += pattern;
    int size = str.size();
    for (int i = 0; i < size; i++)
    {
        pos = str.find(pattern, i);
        if (pos < size)
        {
            std::string s = str.substr(i, pos - i);
            result.push_back(s);
            i = pos + pattern.size() - 1;
        }
    }
    return result;
}


double* SaxUtil::readFromstring(string str) {
    vector<string> strings = split(str, ",");
    double* ret = new double[strings.size()];

    for (int i = 0; i < strings.size(); i++) {
        string s = strings[i];
        string::size_type size_type;
        if (s.length() > 0)
            ret[i] = stod(s, &size_type);
    }

    return ret;
}

double * SaxUtil::readDoubleFromFileAtOnce() {

    INIReader reader("../config.ini");

    if (reader.ParseError() < 0) {
        cout << "Can't load '.ini'\n";
        exit(-1);
    }
    string breakpointsfn = reader.Get("other", "breakpointsfn","");
    cout << "breakpointsfn: " << breakpointsfn <<endl;



    string line;
    std::ifstream fin(breakpointsfn);
    getline(fin, line);

    return readFromstring(line);
}

void SaxUtil::id2Sax2(int id, unsigned short *sax, int segment_num){
    for(int i=0;i<segment_num;++i){
        int temp = id >> 1;
        sax[segment_num - i - 1] = (sax[segment_num - i - 1] << 1) + ((temp + temp == id) ? 0 : 1);
        id = temp;
    }
}

float * SaxUtil::paaFromTs(const float* ts, int tsLengthPerSegment, int segmentNum){
    // Create PAA representation
    auto* paa = new float [segmentNum];

    int s, i;
    for (s=0; s<segmentNum; s++) {
        paa[s] = 0;
        for (i=0; i<tsLengthPerSegment; i++) {
            paa[s] += ts[(s * tsLengthPerSegment)+i];
        }
        paa[s] /= tsLengthPerSegment;
    }
    return paa;
}

void SaxUtil::paaFromTs(const float* ts, float *paa, int tsLengthPerSegment, int segmentNum){
    // Create PAA representation
    double tmp[segmentNum];
    int s, i;
    for (s=0; s<segmentNum; s++) {
        tmp[s] = 0;
        for (i=0; i<tsLengthPerSegment; i++) {
            tmp[s] += ts[(s * tsLengthPerSegment)+i];
        }
        tmp[s] /= tsLengthPerSegment;
        paa[s] = tmp[s];
    }
}

void SaxUtil::paaAndSaxFromTs(const float* ts, float *paa, unsigned short *sax, int tsLengthPerSegment, int segmentNum, int cardinality){
    // Create PAA representation
    int offset = ((cardinality - 1) * (cardinality - 2)) / 2;
    int s, i;
    for (s=0; s<segmentNum; s++) {
        paa[s] = 0;
        for (i=0; i<tsLengthPerSegment; i++) {
            paa[s] += ts[(s * tsLengthPerSegment)+i];
        }
        paa[s] /= tsLengthPerSegment;
        int index = findFirstGE(breakpoints, offset, cardinality -1, paa[s]);

        if(index >= 0)
            sax[s] = (index - offset);
        else
            cout<<"ERROR!!!!!!!";
    }
}

vector<int> * SaxUtil::saxFromTs(float*ts, int tsLengthPerSegment, int segmentNum, int cardinality)
{
    // Create PAA representation
    float * paa = paaFromTs(ts, tsLengthPerSegment, segmentNum);

    // Convert PAA to SAX
    // Note: Each cardinality has cardinality - 1 break points if c is cardinality
    //       the breakpoints can be found in the following array positions:
    //       FROM (c - 1) * (c - 2) / 2
    //       TO   (c - 1) * (c - 2) / 2 + c - 1
    int offset = ((cardinality - 1) * (cardinality - 2)) / 2;
    //printf("FROM %lf TO %lf\n", sax_breakpoints[offset], sax_breakpoints[offset + cardinality - 2]);

    auto*sax = new vector<int>(segmentNum);
//    int* sax = new int[segmentNum];
    int si;
    for (si=0; si<segmentNum; si++) {
        (*sax)[si] = 0;

        // First object = sax_breakpoints[offset]
        // Last object = sax_breakpoints[offset + cardinality - 2]
        // Size of sub-array = cardinality - 1
        int index = findFirstGE(breakpoints, offset, cardinality -1, paa[si]);

        if(index >= 0)
            (*sax)[si] = (index - offset);
        else
            cout<<"ERROR!!!!!!!";

    }

    delete[] paa;
    //sax_print(sax_out, segments, cardinality);
    return sax;
}

void SaxUtil::saxFromTs(const float*ts, unsigned short *sax, int tsLengthPerSegment, int segmentNum, int cardinality)
{
    // Create PAA representation
    double paa[segmentNum];
    int s, i;
    for (s=0; s<segmentNum; s++) {
        paa[s] = 0;
        for (i=0; i<tsLengthPerSegment; i++) {
            paa[s] += ts[(s * tsLengthPerSegment)+i];
        }
        paa[s] /= tsLengthPerSegment;
    }

    // Convert PAA to SAX
    // Note: Each cardinality has cardinality - 1 break points if c is cardinality
    //       the breakpoints can be found in the following array positions:
    //       FROM (c - 1) * (c - 2) / 2
    //       TO   (c - 1) * (c - 2) / 2 + c - 1
    int offset = ((cardinality - 1) * (cardinality - 2)) / 2;
    //printf("FROM %lf TO %lf\n", sax_breakpoints[offset], sax_breakpoints[offset + cardinality - 2]);

    int si;
    for (si=0; si<segmentNum; si++) {
        // First object = sax_breakpoints[offset]
        // Last object = sax_breakpoints[offset + cardinality - 2]
        // Size of sub-array = cardinality - 1
        int index = findFirstGE(breakpoints, offset, cardinality -1, paa[si]);

        if(index >= 0)
            sax[si] = (index - offset);
        else
            cout<<"ERROR!!!!!!!";

    }
    //sax_print(sax_out, segments, cardinality);
}

vector<unsigned short> * SaxUtil::saxFromPaa(float *paa, int segmentNum, int cardinality)
{

    // Convert PAA to SAX
    // Note: Each cardinality has cardinality - 1 break points if c is cardinality
    //       the breakpoints can be found in the following array positions:
    //       FROM (c - 1) * (c - 2) / 2
    //       TO   (c - 1) * (c - 2) / 2 + c - 1
    int offset = ((cardinality - 1) * (cardinality - 2)) / 2;
    //printf("FROM %lf TO %lf\n", sax_breakpoints[offset], sax_breakpoints[offset + cardinality - 2]);

    auto*sax = new vector<unsigned short >(segmentNum);
    //    int* sax = new int[segmentNum];
    int si;
    for (si=0; si<segmentNum; si++) {
        (*sax)[si] = 0;

        // First object = sax_breakpoints[offset]
        // Last object = sax_breakpoints[offset + cardinality - 2]
        // Size of sub-array = cardinality - 1
        int index = findFirstGE(breakpoints, offset, cardinality -1, paa[si]);

        if(index >= 0)
            (*sax)[si] = (index - offset);
        else
            cout<<"ERROR!!!!!!!";

    }
    //sax_print(sax_out, segments, cardinality);
    return sax;
}

int SaxUtil::findFirstGE(const double* array, int start, int length, double target) // satisfy condition: array[?] >= target  and the first one
{
    int end = start + length - 1;
    while (start <= end) {
        int mid = (start + end) / 2;
        if (array[mid] < target)
            start = mid + 1;
        else if (array[mid] >= target)
            end = mid - 1;
    }
    if (end == start + length - 1)
        return -1;
    return end + 1;

}

int SaxUtil::findFirstGE(const int* array, int start, int length, int target) // satisfy condition: array[?] >= target  and the first one
{
    int end = start + length - 1;
    while (start <= end) {
        int mid = (start + end) / 2;
        if (array[mid] < target)
            start = mid + 1;
        else if (array[mid] >= target)
            end = mid - 1;
    }
    if (end == start + length - 1)
        return -1;
    return end + 1;

}

int SaxUtil::invSaxHeadFromSax(vector<int> *sax, int bitsCardinality, int segmentNum)
{
    int i=bitsCardinality-1,s = 0,n;
    for (int j=0; j < segmentNum; j++)
    {
        n = (*sax)[j];
        n >>= i;
        s |= (n % 2) << (segmentNum - 1 - j);
    }
    return s;
}

int SaxUtil::invSaxHeadFromSax(const unsigned short *sax, int bitsCardinality, int segmentNum)
{
    int i=bitsCardinality-1,s = 0,n;
    for (int j=0; j < segmentNum; j++)
    {
        n = sax[j];
        n >>= i;
        s |= (n % 2) << (segmentNum - 1 - j);
    }
    return s;
}

int SaxUtil::invSaxHeadkFromSax(const unsigned short *sax, int bitsCardinality, int segmentNum, int k)
{
    int i=bitsCardinality-k,s = 0,n;
    for (int j=0; j < segmentNum; j++)
    {
        n = sax[j];
        n >>= i;
        s |= (n % 2) << (segmentNum - 1 - j);
    }
    return s;
}

int SaxUtil::invSaxHeadFromPaa(const float *paa, int tsLengthPerSegment, int segmentNum) {
    int res = 0;
    for (int s=0; s<segmentNum; s++) {
        res <<= 1;
        if(paa[s] >= 0)
            res++;
    }
    return res;
}

double SaxUtil::LowerBound_Paa_iSax(const float *paa, const unsigned short *sax){
    double frontCoef = (double)Const::tsLength / Const::segmentNum; // n / w
    int offset = Const::offset, cardinality = Const::cardinality;
    double paaValue, lb, ub;
    unsigned short saxValue;
    double sum = 0;

    for (int i = 0; i < Const::segmentNum; i++) {
        paaValue = paa[i];
        saxValue = sax[i];
        if(saxValue == 0){
            lb = -numeric_limits<double>::max();
            ub = breakpoints[offset];
        }
        else if(saxValue == cardinality - 1){
            lb = breakpoints[offset + cardinality - 2];
            ub = numeric_limits<double>::max();
        }
        else {
            lb = breakpoints[offset + saxValue - 1];
            ub = breakpoints[offset + saxValue];
        }

        if(paaValue < lb){
            sum += (lb - paaValue) * (lb - paaValue);
        }
        else if(paaValue > ub){
            sum += (paaValue - ub) * (paaValue - ub);
        }
    } // for

    return frontCoef * sum;
}

double SaxUtil::LowerBound_Paa_iSax(const float *paa, const unsigned short *sax, const int* bits_cardinality, vector<int>&chosen_segs, int new_id){
    double frontCoef = (double)Const::tsLength / Const::segmentNum; // n / w
    double paaValue, lb, ub;
    double sum = 0;
    int saxValue, cur = chosen_segs.size() - 1, bc;

    for (int i = 0; i < Const::segmentNum; i++) {
        paaValue = paa[i];
        if(chosen_segs[cur] == i){
            saxValue = (sax[i] << 1) + (new_id % 2);
            new_id >>= 1;
            bc = bits_cardinality[i] + 1;
        } else {
            saxValue = sax[i];
            bc = bits_cardinality[i];
        }
        getValueRange(saxValue, bc, &lb, &ub);

        if(paaValue < lb){
            sum += (lb - paaValue) * (lb - paaValue);
        }
        else if(paaValue > ub){
            sum += (paaValue - ub) * (paaValue - ub);
        }
    } // for

    return frontCoef * sum;
}

double SaxUtil::LowerBound_Paa_iSax(const float *paa, const unsigned short *sax, int bits_cardinality) {
    double frontCoef = (double)Const::tsLength / Const::segmentNum; // n / w
    long offset; int cardinality;
    if(bits_cardinality == Const::bitsCardinality)
        offset = Const::offset, cardinality = Const::cardinality;
    else
        cardinality = 1 << bits_cardinality, offset = ((long )(cardinality - 1) * (cardinality - 2)) / 2;
    double paaValue, lb, ub;
    int saxValue;
    double sum = 0;

    for (int i = 0; i < Const::segmentNum; i++) {
        paaValue = paa[i];
        saxValue = sax[i];
        if(saxValue == 0){
            lb = -numeric_limits<double>::max();
            ub = breakpoints[offset];
        }
        else if(saxValue == cardinality - 1){
            lb = breakpoints[offset + cardinality - 2];
            ub = numeric_limits<double>::max();
        }
        else {
            lb = breakpoints[offset + saxValue - 1];
            ub = breakpoints[offset + saxValue];
        }

        if(paaValue < lb){
            sum += (lb - paaValue) * (lb - paaValue);
        }
        else if(paaValue > ub){
            sum += (paaValue - ub) * (paaValue - ub);
        }
    } // for

    return frontCoef * sum;
}

double SaxUtil::LowerBound_Paa_iSax(const float *paa, const unsigned short *sax, const int* bits_cardinality) {
    double frontCoef = (double)Const::tsLength / Const::segmentNum; // n / w
    double paaValue, lb, ub;
    double sum = 0;

    for (int i = 0; i < Const::segmentNum; i++) {
        paaValue = paa[i];
        if(bits_cardinality[i] == 0)    continue;
        getValueRange(sax[i], bits_cardinality[i], &lb, &ub);

        if(paaValue < lb){
            sum += (lb - paaValue) * (lb - paaValue);
        }
        else if(paaValue > ub){
            sum += (paaValue - ub) * (paaValue - ub);
        }
    } // for

    return frontCoef * sum;
}

void SaxUtil::getValueRange(int sax_single, int bits_cardinality, double *lb, double *ub){
    int cardinality = 1 << bits_cardinality;
    int offset = ((cardinality - 1) * (cardinality - 2)) / 2;
    if(sax_single == 0){
        *lb = -numeric_limits<double>::max();
        *ub = breakpoints[offset];
    }else if(sax_single == cardinality - 1){
        *lb = breakpoints[offset + sax_single - 1];
        *ub = numeric_limits<double>::max();
    }else{
        *lb = breakpoints[offset + sax_single - 1];
        *ub = breakpoints[offset + sax_single];
    }
}

// return the new id
int SaxUtil::extendSax(float *paa, const int *bits_cardinality, vector<int> &segments) {
    int res = 0, cardinality, sw, offset;
    for(int segment:segments){
        cardinality = 1 << (bits_cardinality[segment] + 1);
        offset = ((cardinality - 1) * (cardinality - 2)) / 2;
        int index = findFirstGE(breakpoints, offset, cardinality -1, paa[segment]);
        if(index >= 0)
            sw = (index - offset);
        else
            cout<<"ERROR!!!!!!!";
        res = (res << 1) + (sw % 2);
    }
    return res;
}

// return the new id
int SaxUtil::extendSax(const unsigned short *sax, const int *bits_cardinality, vector<int> &segments) {
    int res = 0, sw;
    for(int segment:segments){
        sw = sax[segment] >> (Const::bitsCardinality - bits_cardinality[segment] - 1);
        res = (res << 1) + (sw % 2);
    }
    return res;
}

// return the new id
int SaxUtil::extendSax(const unsigned short *sax, const int *bits_cardinality, vector<int> &segments,
                       const unsigned short *parent_sax) {
    int res = 0, sw;
    for(int segment:segments){
        sw = sax[segment] >> (Const::bitsCardinality - bits_cardinality[segment] - 1);
        if((sw >> 1) == parent_sax[segment])    res = (res << 1) + (sw % 2);
        else if((sw >> 1) > parent_sax[segment])    res = (res << 1) + 1;
        else res <<= 1;
    }
    return res;
}

// return the new id
int SaxUtil::extendSax(const unsigned short *sax, const int *bits_cardinality) {
    int res = 0, sw;
    for(int segment = 0;segment < Const::segmentNum; ++segment){
        sw = sax[segment] >> (Const::bitsCardinality - bits_cardinality[segment] - 1);
        res = (res << 1) + (sw % 2);
    }
    return res;
}

int SaxUtil::getNewId(const float *paa, const float *split_line){
    int ret = 0;
    for(int i=0;i<Const::segmentNum;++i){
        ret <<= 1;
        if(paa[i] > split_line[i])  ret += 1;
    }
    return  ret;
}

int SaxUtil::getNewId(const float *paa, const float *split_line, vector<int>&segments){
    int ret = 0;
    for(int i:segments){
        ret <<= 1;
        if(paa[i] > split_line[i])  ret += 1;
    }
    return  ret;
}

void SaxUtil::saxPrint(int* sax, int bits_cardinality, int segment_num)
{
    int i;
    for (i=0; i < segment_num; i++) {
        cout<<i << "\t";
        printBinary(sax[i], bits_cardinality);
        cout<<endl;
    }
    cout<<endl;
}

void SaxUtil::printBinary(long n, int size) {
    char b[size + 1];
    b[size] = 0;
    int i;
    for (i=0; i<size; i++) {
        b[i] = '0';
    }

    for (i=0; i<size; i++, n=n/2)
        if (n%2 == 1) b[size-1-i] = '1';

    cout<<b;
}

void SaxUtil::generateSaxFile(const string &fn , const string &output){
    long fs = FileUtil::getFileSize(fn.c_str());
    int num = fs / Const::tsLengthBytes;
    cout << "Total Number is "<< num <<endl;
    float ts[Const::tsLength];
    unsigned short sax[Const::segmentNum];
    FILE *f = fopen(fn.c_str(), "rb");
    FILE *of = fopen(output.c_str(), "wb");
    for(int i=0;i<num;++i){
        if(i % 1000000 == 0)    cout << i << endl;
        if(i == 74625529)
            cout <<"hereh"<<endl;
        fread(ts, sizeof(float ), Const::tsLength, f);
        if(isnan(ts[0])){
            for(auto &t:sax)    t = 0;
            cout << i << "," <<endl;
        }
        else    saxFromTs(ts, sax, Const::tsLengthPerSegment, Const::segmentNum, Const::cardinality);
        fwrite(sax, sizeof(unsigned short), Const::segmentNum, of);
    }
    fclose(f);
    fclose(of);
}

void SaxUtil::generatePaaFile(const string &fn, const string &output) {
    long fs = FileUtil::getFileSize(fn.c_str());
    int num = fs / Const::tsLengthBytes;
    cout << "Total number is " << num << endl;
    float ts[Const::tsLength];
    float paa[Const::segmentNum];
    FILE *f = fopen(fn.c_str(), "rb");
    FILE *of = fopen(output.c_str(), "wb");
    for(int i=0;i<num;++i){
        if(i % 1000000 == 0)    cout << i << endl;
        if(i == 74625529)
            cout << i << endl;
        fread(ts, sizeof(float ), Const::tsLength, f);
        if(isnan(ts[0])){
            for(auto &t:paa)    t = 0;
            cout << i << "," <<endl;
        }
        else    paaFromTs(ts, paa, Const::tsLengthPerSegment, Const::segmentNum);
        fwrite(paa, sizeof(float ), Const::segmentNum, of);
    }
    fclose(f);
    fclose(of);
}