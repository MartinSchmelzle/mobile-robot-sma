#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <string>
#include <algorithm>
#include <filesystem>

using h_map=std::unordered_map<std::string, std::vector<double>>;

// Function to convert a string containing numbers to a vector of doubles
std::vector<double> stringToDoubleVector(const std::string& inputString) {
    std::vector<double> result;
    std::istringstream iss(inputString);
    std::string token;

    // Check if the string contains multiple numbers separated by spaces
    while (iss >> token) {
        double num;
        std::istringstream(token) >> num;
        result.push_back(num);
    }

    return result;
}

//function for splitting a string between left of the comma and right of the comma
void splitString(const std::string& inputString, std::string& leftPart, std::string& rightPart) {
    // Find the position of the comma in the string
    size_t commaPos = inputString.find(',');

    if (commaPos != std::string::npos) {
        // Extract left part (substring before comma)
        leftPart = inputString.substr(0, commaPos);

        // Extract right part (substring after comma)
        rightPart = inputString.substr(commaPos + 1);
    } else {
        // If comma is not found, set both parts as empty strings
        leftPart = "";
        rightPart = "";
    }
}

//function that unpacks h_data from the thing MATLAB returned into h
h_map h_unpack() {
    std::ifstream file("h_data.csv");
    h_map hdata;
    if (!file.is_open()) {
        std::cerr << "Unable to open file!" << std::endl; //error handling
        // Handle file opening error here if needed
    }

    std::vector<std::string> upper = {"c100", "c10", "c1", "c20", "c2", "c3", "c4", "c6", "s", "c0", "c40", "c5", "c7", "c8", "c9", "dock", "updock", "b", "R", "Accel", "Slow"};//all the structs directly below h
    std::vector<std::string> upper3 = {"c100", "c10", "c1", "c20", "c2", "c3", "c4", "c6", "s", "c0", "c40", "c5", "c7", "c8", "c9"};//structs with three layers (has substructs)
    std::vector<std::string> upper2 = {"dock", "updock", "b", "R", "Accel", "Slow"};//structs with two layers (don't have other structs beneath them)
    std::vector<std::string> middle = {"Acc","Slw","Con"}; //middle layer

    std::string line;
    std::string cur_upper;
    std::string cur_mid;
    std::string cur_low;
    std::string leftOfComma;
    std::string rightOfComma;
    std::vector<double> cur_value;
    std::string key;

    while (std::getline(file, line)) { //until we run out of lines
        if (line.empty()) {//skip empty lines
            continue;
        }

        auto k1 = std::find(upper.begin(), upper.end(), line); //find whether line is in upper
        auto k2 = std::find(upper2.begin(), upper2.end(), cur_upper); //find whether current upper is in two- or three-layered upper
        auto k3 = std::find(upper3.begin(), upper3.end(), cur_upper);

        if (k1 != upper.end()) { //if line is in upper, update current upper
            cur_upper = line;
            //std::cout << cur_upper << std::endl;
        } else if (k3 != upper3.end()) { //if there are three layers
            auto k4 = std::find(middle.begin(), middle.end(), line); //find out if line is in middle
            if(k4!=middle.end()){ //if line is in middle struct, update cur_mid
                cur_mid=line;
                //std::cout<<cur_upper<<"."<<cur_mid<<std::endl;
            }
            else{
                splitString(line,leftOfComma,rightOfComma);//separate the two sides of the comma
                cur_low=leftOfComma;
                key = cur_upper +"."+ cur_mid +"."+ cur_low;
                //std::cout<<cur_upper<<"."<<cur_mid<<"."<<cur_low<<std::endl;
                cur_value=stringToDoubleVector(rightOfComma);
                hdata[key] = cur_value;

            }
        } else if (k2 != upper2.end()) {
            splitString(line,leftOfComma,rightOfComma);//separate the two sides of the comma
            cur_low=leftOfComma;
            cur_value=stringToDoubleVector(rightOfComma);
            key = cur_upper +"."+ cur_low;
            hdata[key] = cur_value;
            //std::cout<<cur_upper<<"."<<cur_low<<std::endl;
        }
    }

    return hdata;
}

//function that writes the unordered_map back into another csv file h_Data (for better runtime)
void exportToCSV(const h_map& myMap, const std::string& filename) {
    std::ofstream csvFile(filename);

    if (csvFile.is_open()) {
        // Write headers
        csvFile << "Key,Values" << std::endl;

        // Iterate through the map and write key-value pairs to the CSV file
        for (const auto& pair : myMap) {
            csvFile << pair.first << ",";

            // Write values from the vector
            for (size_t i = 0; i < pair.second.size(); ++i) {
                csvFile << pair.second[i];
                if (i != pair.second.size() - 1) {
                    csvFile << " "; // Separate values by space
                }
            }

            csvFile << std::endl;
        }

        std::cout << "Data has been exported to " << filename << " successfully." << std::endl;
        csvFile.close();
    } else {
        std::cout << "Unable to open the file " << filename << " for writing." << std::endl;
    }
}

//struct for key + value
struct KeyValue {
    std::string key;
    std::string value;
};

// Function to find a string in a CSV file by key
//Returns a vector of doubles => if you just want a double, add [0]
//examplary usage: double g=findinh("c5.Acc.L")[0];
std::vector<double> findinh(const std::string& searchString) {
    const std::string& filename = "SMA_conf.csv";
    std::ifstream csvFile(filename);

    if (!csvFile.is_open()) {
        throw std::runtime_error("Unable to open the csv file " + filename);
    }

    std::string foundValue;
    std::vector<KeyValue> keyValuePairs;

    // Read each line in the CSV file
    std::string line;
    while (std::getline(csvFile, line)) {
        std::istringstream iss(line);
        std::string key, value;
        if (std::getline(iss, key, ',') && std::getline(iss, value)) {
            // Skip lines containing a comma in the value
            if (value.find(',') != std::string::npos) {
                continue; // Skip this line
            }

            // Store key-value pairs in a vector for future reference
            keyValuePairs.push_back({key, value});

            // Check if the key matches the search string
            if (key == searchString) {
                if (!foundValue.empty()) {
                    throw std::runtime_error("Error: Multiple results found for the search string.");
                }
                foundValue = value; // Store the found value
            }
        }
    }

    // Check if the search string was found
    if (foundValue.empty()) {
        throw std::runtime_error("Error: Search string not found in the CSV file.");
    }
    //transform to double vector
    std::vector<double> valdouble = stringToDoubleVector(foundValue);
    // Return the value corresponding to the found key
    return valdouble;
}
/*
int main()
{
    if(!std::filesystem::exists("h_Data_2.csv")) //If h_Data_2.csv doesn't exist, go through the computationally expensive step of creating it
    {h_map h = h_unpack();
    exportToCSV(h, "h_Data_2.csv");}
    std::vector<double> result = findinh("c5.Con.L");

    //print out result vector
    for (const auto& element : result) {
        std::cout << ", " << element;
    }
    std::cout << std::endl;
    return 0;
}
*/