#ifndef PARALLEL_H
#define PARALLEL_H

#include "junction.hpp"

#include <string>
#include <iostream>
#include <stdio.h>
#include <fstream>

#include <vector>
#include <tuple>
#include <map>

#include <thread>
#include <future>

typedef std::tuple<double, std::map<std::string, double>> fnRes;
class ComputeUtil
{
public:
    static void customResultMap(std::map<std::string, std::vector<double>> resultMap, std::string filename)
    {

        std::ofstream saveFile;

        saveFile.open(filename);
        int logLength = 0;
        for (const auto [key, _] : resultMap)
        {
            if (!logLength)
                logLength = resultMap[key].size();
            saveFile << key << ";";
        }
        saveFile << "\n";
        for (unsigned int i = 0; i < logLength; i++)
        {
            for (const auto [_, value] : resultMap)
            {
                saveFile << value[i] << ";";
            }
            saveFile << "\n";
        }
        saveFile.close();
    }

    static std::map<std::string, std::vector<double>> parallelFieldScan(Junction &mtj, double minField, double maxField,
                                                                        int numberOfPoints,
                                                                        int threadNumber,
                                                                        const std::function<fnRes(Junction &mtj, const double scanningParam)>
                                                                            runnableFunction)
    {
 
        int threadNum = std::thread::hardware_concurrency();
        if (threadNumber > 0)
        {
            threadNum = threadNumber <= threadNum ? threadNumber : threadNum;
        }
        std::vector<std::future<std::vector<fnRes>>> threadResults;
        threadResults.reserve(threadNum);
        std::cout << "Using: " << threadNum << " threads." << std::endl;
        const int pointsPerThread = (numberOfPoints / threadNum) + 1;
        const double spacing = (maxField - minField) / numberOfPoints;
        for (int i = 0; i < threadNum; i++)
        {
            const double threadMinField = pointsPerThread * i * spacing;
            const double threadMaxField = pointsPerThread * (i + 1) * spacing;
            threadResults.emplace_back(std::async([mtj, runnableFunction, threadMinField, threadMaxField, spacing]() mutable {
                std::vector<fnRes>
                    resAcc;
                for (double field = threadMinField; field < threadMaxField;
                     field += spacing)
                {
                    auto subRes = runnableFunction(mtj, field);
                    resAcc.push_back(subRes);
                }
                return resAcc;
            }));
        }

        std::map<std::string, std::vector<double>> finalRes;
        for (auto &result : threadResults)
        {
            for (const auto [field, map] : result.get())
            {
                finalRes["H"].push_back(field);
                for (const auto [key, value] : map)
                {
                    finalRes[key].push_back(value);
                }
            }
        }
        return finalRes;
    }
};

#endif
