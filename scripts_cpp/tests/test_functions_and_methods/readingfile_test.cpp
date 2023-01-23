#include <fstream>
#include <iostream>
#include <string>
#include <vector>

int main(){

    std::cout<<"\n===========TEST READING LINES==============\n";
    // Store the contents into a vector of strings
    std::vector<std::string> outputs;

    std::cout << "Reading from readingfile_test.txt....\n";

    // Create the file object (input)
    std::ifstream infile("readingfile_test.txt");

    // Temporary buffer
    std::string temp;

    // Get the input from the input file until EOF
    while (std::getline(infile, temp)) {
        // Add to the list of output strings
        outputs.push_back(temp);
    }

    // Check output vector
    for (const auto& i : outputs)
        std::cout << i << std::endl;

    std::cout<<"\n===========TEST READING NON EXISTING FILE==============\n";
    // Store the contents into a vector of strings
    std::vector<std::string> outputs2;

    // Create the file object (input)
    std::ifstream infile2("readingfile2_test.txt");

    // Temporary buffer
    std::string temp2;

    // Get the input from the input file until EOF
    while (std::getline(infile2, temp2)) {
        // Add to the list of output strings
        outputs2.push_back(temp2);
    }

    // Check output vector
    for (const auto& i : outputs2)
        std::cout << i << std::endl;
    
    std::cout<<"Length of outputs2 is "<<outputs2.size();


    return 0;
}