#ifndef GENERAL_H
#define GENERAL_H
#include "TChain.h"
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include "TSystem.h"
#include <algorithm>

using namespace std;

namespace ImagingAlgorithm
{
double ErrorDouble();
int ErrorInt();

// Generate Chain with treename, keywords vector, already exsited TChain*(optional), folder name(optional) as argument
TChain *GenerateChain(std::string treeName, const std::vector<std::string> &keyWord, TChain *chain = NULL, string sFolder = "./", std::vector<std::string> notInclude = {});

// Find ele in vector<T>, return index, -1 if not found
template <typename T>
int FindInVector(const vector<T> &Array, const T &Element);

// Arrange vector<T> array, return map vector
template <typename T>
std::vector<pair<int, T>> EncodeArray(const vector<T> &array);
template <typename T>
std::vector<int> DecodeArray(const std::vector<pair<int, T>> &arrayEncode);
template <typename T>
std::vector<int> SortVectorDescend(const vector<T> &array);
template <typename T>
std::vector<int> SortVectorAscend(const vector<T> &array);

template <typename T>
void TryPrintVector(const vector<T> &array, ostream &os = std::cout);

template <typename T>
const T &MaxTemplate(const T &a, const T &b);
template <typename T>
const T &MinTemplate(const T &a, const T &b);

} // namespace ImagingAlgorithm

template <typename T>
int ImagingAlgorithm::FindInVector(const vector<T> &array, const T &ele)
{
    for (int i = 0; i < array.size(); i++)
    {
        if (array[i] == ele)
        {
            return i;
        }
    }
    return -1;
}

template <typename T>
std::vector<pair<int, T>> ImagingAlgorithm::EncodeArray(const vector<T> &array)
{
    vector<pair<int, T>> encodeArray;
    for (int i = 0; i < array.size(); i++)
    {
        encodeArray.push_back({i, array[i]});
    }
    return encodeArray;
}

template <typename T>
std::vector<int> ImagingAlgorithm::DecodeArray(const std::vector<pair<int, T>> &arrayEncode)
{
    vector<int> arrayDecode;
    for (int i = 0; i < arrayEncode.size(); i++)
    {
        arrayDecode.push_back(arrayEncode[i].first);
    }
    return arrayDecode;
}

template <typename T>
std::vector<int> ImagingAlgorithm::SortVectorDescend(const vector<T> &array)
{
    auto encodeArray = ImagingAlgorithm::EncodeArray(array);
    std::sort(encodeArray.begin(), encodeArray.end(), [](pair<int, T> a, pair<int, T> b) { return b.second < a.second; });
    auto decodeArray = ImagingAlgorithm::DecodeArray(encodeArray);
    return decodeArray;
}

template <typename T>
std::vector<int> ImagingAlgorithm::SortVectorAscend(const vector<T> &array)
{
    auto encodeArray = ImagingAlgorithm::EncodeArray(array);
    std::sort(encodeArray.begin(), encodeArray.end(), [](pair<int, T> a, pair<int, T> b) { return b.second > a.second; });
    auto decodeArray = ImagingAlgorithm::DecodeArray(encodeArray);
    return decodeArray;
}

template <typename T>
void ImagingAlgorithm::TryPrintVector(const vector<T> &array, ostream &os)
{
    for (int i = 0; i < array.size(); i++)
    {
        os << array[i] << '\t';
    }
}

template <typename T>
const T &ImagingAlgorithm::MaxTemplate(const T &a, const T &b)
{
    if (a > b)
        return a;
    else
        return b;
}

template <typename T>
const T &ImagingAlgorithm::MinTemplate(const T &a, const T &b)
{
    if (a < b)
        return a;
    else
        return b;
}

#define gErrorDouble (ImagingAlgorithm::ErrorDouble())
#define gErrorInt (ImagingAlgorithm::ErrorInt())

#endif