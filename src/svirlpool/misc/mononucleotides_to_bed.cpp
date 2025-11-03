// this is a small program to write a bed file that 
// contains all mononucleotide stretches of a minimum given length
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>

bool inAlphabet(char c){
    return (c == 'A' || c == 'C' || c == 'G' || c == 'T' ||
            c == 'a' || c == 'c' || c == 'g' || c == 't');
}

int main(int argc, char *argv[]) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <input_fasta> <output_bed> <N>\nWhere N is the minimum length of the stretch\n";
        return 1;
    }

    std::ifstream input(argv[1]);
    if (!input.is_open()) {
        std::cerr << "Error opening input file\n";
        return 1;
    }

    std::ofstream output(argv[2]);
    if (!output.is_open()) {
        std::cerr << "Error opening output file\n";
        input.close();
        return 1;
    }

    int N = std::stoi(argv[3]);
    char last_char = '\0';
    int stretch = 0;
    std::string chrName;
    int i = 0;

    std::string line;
    float totalLines = 0.0;
    while (std::getline(input, line)) {
        totalLines++;
    }

    // Return to the beginning of the file
    input.clear();
    input.seekg(0);

    while (std::getline(input, line)) {
        if (!line.empty() && line[0] == '>') {
            // let chrName be the first word of the line withput the first character
            std::istringstream iss(line);
            iss >> chrName;
            chrName = chrName.substr(1);
            // print name of chromosome
            std::cout << "working on chromsome " << chrName << std::endl;
            i = 0;
            last_char = '\0';
            stretch = 0;
        } else {
            for (char current_char : line) {
                if (current_char == last_char && inAlphabet(current_char)) {
                    stretch++;
                } else {
                    if (stretch >= N) {
                        output << chrName << '\t' << i - stretch << '\t' << i << '\t' << stretch << '\n';
                    }
                    stretch = 1;
                    last_char = current_char;
                }
                i++;
            }
        }
    }

    // Check if there's a final stretch longer than N
    if (stretch >= N) {
        output << chrName << '\t' << i - stretch << '\t' << i << '\t' << stretch << '\n';
    }

    // Close files
    input.close();
    output.close();

    std::cout << std::endl; // Print a newline to clear the progress bar
    return 0;
}