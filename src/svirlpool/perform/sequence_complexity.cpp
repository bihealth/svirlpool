#include <iostream>
#include <vector>
#include <cstdlib>
#include <string>
#include <cmath>
#include <set>
#include <algorithm>

size_t hash_dna5(char c)
{
    switch (c)
    {
        case 'A':
            return 0;
        case 'a':
            return 0;
        case 'C':
            return 1;
        case 'c':
            return 1;
        case 'G':
            return 2;
        case 'g':
            return 2;
        case 'T':
            return 3;
        case 't':
            return 3;
        default:
            return 4;
    }
}

// this function recieves a start iterator and an end iterator
// and returns a hash value for the sequence between the two iterators
size_t kmer_to_hash(std::string::iterator it_start, std::string::iterator it_end, size_t base=3)
{
    size_t hashvalue = 0;
    for (auto it = it_start; it != it_end; it++)
    {
        hashvalue = (hashvalue << base) + hash_dna5(*it);
    }
    return hashvalue;
}

// this function is ran once at initialization
// it computes the maximum number of unique hashes for each kmer
// in a given sequence of characters
std::vector<size_t> sequence_to_kmer_hashes(
    std::string::iterator it_start,
    std::string::iterator it_end,
    size_t k)
{
    // loop over the sequence and calculate the hash value for each kmer
    std::vector<size_t> hashvalues{};
    for (auto it = it_start; it != it_end-k+1; it++)
    {
        hashvalues.push_back(kmer_to_hash(it,it+k));
    }
    return hashvalues;
}

size_t extend_hash(size_t kmer_hash, size_t letter_dna5 , size_t k, size_t base=3)
{
    return ((kmer_hash << base) + letter_dna5) & ((1 << (k*base))-1);
}


// this function computes the product of the element wise division of two arrays.
// the arrays must have the same size. The result is a float.
float product(std::vector<size_t>& a, std::vector<size_t>& b, size_t size)
{
    float result(1.0);
    for (size_t i = 0; i < size; i++)
    {
        result *= (float(a[i]) / float(b[i]));
    }
    return result;
}

size_t count_equal_kmer_hashes(std::vector<size_t>& kmer_hashes, size_t start, size_t end, size_t value)
{
    if (start >= kmer_hashes.size() || end > kmer_hashes.size()) {
        throw std::runtime_error("Error: start or end index out of bounds in count_equal_kmer_hashes.");
    }
    size_t count(0);
    for (size_t i = start; i < end; i++)
    {
        if (kmer_hashes[i] == value)
        {
            count++;
        }
    }
    return count;
}

size_t calc_position(size_t w, size_t k, size_t i, size_t index)
{
    size_t alpha = w-k+1;
    return index+((alpha+i-1)%alpha);
}
size_t calc_last_position(size_t w, size_t k, size_t i, size_t index)
{
    size_t alpha = w-k+1;
    return index+((alpha+i-2)%alpha);
}



std::vector<float> run_program(
        size_t wsize,
        std::vector<uint8_t> kmers,
        std::string dna)
{
    // check if the window size is larger than the maximum kmer size
    if (wsize < *std::max_element(kmers.begin(), kmers.end()))
    {
        throw std::runtime_error("The window size must be larger than the maximum kmer size.");
    }
    // check if wsize is odd and 2 < wsize < 21.
    if (wsize % 2 == 0 || wsize < 2 || wsize > 21)
    {
        throw std::runtime_error("The window size must be an odd number between 2 and 21.");
    }
    // check if there are more than 0 kmers and if their size ranges between 3 and wsize
    if (kmers.size() == 0)
    {
        throw std::runtime_error("At least one kmer size must be provided.");
    }
    if (*std::min_element(kmers.begin(), kmers.end()) < 3 || *std::max_element(kmers.begin(), kmers.end()) > wsize)
    {
        throw std::runtime_error("The kmer sizes must be between 3 and the window size.");
    }
    // check if all k-mers are in monotonic ascending order
    for (size_t i = 1; i < kmers.size(); i++)
    {
        if (kmers[i-1] < kmers[i])
        {
            throw std::runtime_error("The kmer sizes must be in strict ascending order.");
        }
    }
    // check if the dna sequence is not empty
    if (dna.size() == 0)
    {
        throw std::runtime_error("The DNA sequence must not be empty.");
    }

    size_t nk(kmers.size());
    size_t kmers_array[nk] = {0};
    std::copy(kmers.begin(), kmers.end(), kmers_array);
    std::vector<size_t> MAX_UNIQUE_HASHES(nk, 0);
    // compute the maximum number of unique hashes for each k
    for (size_t i = 0; i < nk; i++){
        MAX_UNIQUE_HASHES[i] = std::min(wsize-kmers_array[i]+1,size_t(pow(4,kmers_array[i]))+1);
    }
    // replace the vector of vectors with an array. every array according to a kmer is saved in
    // a certain interval on the array. this way, we can use the same array for all kmers.
    std::vector<size_t> kmer_hashes_index(nk, 0);
    for (size_t i = 1; i < nk; i++)
    {
        kmer_hashes_index[i] = kmer_hashes_index[i-1]
            + wsize - kmers_array[i-1]+1;
    }
    size_t kmer_hashes_size = kmer_hashes_index[nk-1]+wsize-kmers_array[nk-1]+1;
    
    
    
    std::vector<size_t> kmer_hashes(kmer_hashes_size, 0);
    
    
    
    // std::vector<std::vector<size_t>> kmer_hashes(kmers.size());
    std::vector<size_t> current_unique_n_hashes(nk,0);

    for (size_t ki = 0; ki < nk; ki++)
    {
        std::vector<size_t> kmer_hashes_init = sequence_to_kmer_hashes(dna.begin(),dna.begin()+wsize,kmers_array[ki]);
        // i is the index of the kmer in the window - offset
        // the offset is given in the index
        // std::cout << "calc first hashes with k = " << kmers_array[ki] << std::endl;
        for (size_t i = 0; i < wsize-kmers_array[ki]+1; i++)
        {
            kmer_hashes[kmer_hashes_index[ki] + i] = kmer_hashes_init[i];
        }
    }

    // Starting unique k-mer counts
    for (size_t ki = 0; ki < nk; ki++)
    {
        size_t interval_start = kmer_hashes_index[ki];
        size_t interval_end = interval_start + wsize - kmers_array[ki] + 1;
        if (interval_end > kmer_hashes.size()) {
            std::cerr << "Error: interval_end out of bounds!" << std::endl;
            exit(1);
        }
        std::set<size_t> unique_hashes(kmer_hashes.begin() + interval_start, 
                                    kmer_hashes.begin() + interval_end);
        current_unique_n_hashes[ki] = unique_hashes.size();
    }

    std::vector<float> results(dna.size(),0.0);

    float p = product(current_unique_n_hashes,MAX_UNIQUE_HASHES,nk);
    results[0] = p;

    
    // main loop
    // while loop through the rest of standard input
    // create the return vector<float>
    
    // iterate over lines instead of std::getline
    // iterate over dna from position wsize/2 to end using iterators
    size_t i(0);
    for (auto it = dna.begin()+wsize/2; it != dna.end(); it++)
    {
        i++;
        size_t h = hash_dna5(*it);
        for (size_t j = 0; j < nk; j++)
        {
            size_t k = kmers_array[j];
            size_t position = calc_position(wsize,k,i,kmer_hashes_index[j]);
            if (position >= kmer_hashes_size) {throw std::runtime_error("Position exceeds kmer_hashes_size.");}
            size_t last_position=calc_last_position(wsize,k,i,kmer_hashes_index[j]);
            if (last_position >= kmer_hashes_size) {throw std::runtime_error("Last position exceeds kmer_hashes_size.");}
            size_t old_hash = kmer_hashes[last_position];
            bool last_hash_unique = count_equal_kmer_hashes(kmer_hashes,kmer_hashes_index[j],kmer_hashes_index[j]+wsize-k+1,kmer_hashes[position]) == 1;
            current_unique_n_hashes[j] -= size_t(last_hash_unique);
            size_t new_hash = extend_hash(old_hash,h,k);
            kmer_hashes[position]=new_hash;
            // update count of unique elements in k-mers
            bool new_hash_unique = count_equal_kmer_hashes(kmer_hashes,kmer_hashes_index[j],kmer_hashes_index[j]+wsize-k+1,new_hash) == 1;
            current_unique_n_hashes[j] += int(new_hash_unique);
        }
        // write result to vector
        float p = product(current_unique_n_hashes,MAX_UNIQUE_HASHES,nk);
        results[i-1] = p;
    }
    // padding at the end
    //p = product(current_unique_n_hashes,MAX_UNIQUE_HASHES,nk);
    while (i < dna.size())
    {
        results[i] = p;
        i++;
    }
    return results;
}


struct cmd_arguments {
    int w;
    std::vector<uint8_t> k_values;
    std::string dna;
    bool verbose;
};

void print_help() {
    std::cout << "Usage: program_name [-w <odd_integer>] [-k <ascending_integers>]\n"
              << "Options:\n"
              << "  -w   Set an odd integer between 5 and 21 (default: 21)\n"
              << "  -k   Set multiple ascending integers between 2 and w/2 (default: 2 3 4 5 6 7 8 9 10)\n"
              << "  -s   DNA sequence or several sequences, e.g. 'ACGTCGCTGCAT'\n"
              << "  -v   verbosity: print DNA letter, its position and its hash results value\n";
}

void parse_arguments(int argc, char **argv, cmd_arguments &args) {
    args.w = 21;
    args.k_values = {2, 3, 4, 5, 6, 7, 8, 9, 10};
    args.verbose = false;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-w") {
            if (i + 1 < argc) {
                args.w = std::atoi(argv[++i]);
                if (args.w < 5 || args.w > 21 || args.w % 2 == 0) {
                    std::cerr << "Error: Invalid value for -w. Please provide an odd integer between 5 and 21.\n";
                    print_help();
                    std::exit(EXIT_FAILURE);
                }
            } else {
                std::cerr << "Error: -w option requires an argument.\n";
                print_help();
                std::exit(EXIT_FAILURE);
            }
        } else if (arg == "-k") {
            // erase default values
            args.k_values.clear();
            while (i + 1 < argc && argv[i + 1][0] != '-') {
                int k_value = std::atoi(argv[++i]);
                if (k_value < 2 || k_value > args.w / 2) {
                    std::cerr << "Error: Invalid value for -k. Please provide ascending integers between 2 and w/2.\n";
                    print_help();
                    std::exit(EXIT_FAILURE);
                }
                args.k_values.push_back(k_value);
                // check if k_values is in ascending order
                if (args.k_values.size() > 1 && args.k_values[args.k_values.size() - 2] >= args.k_values[args.k_values.size() - 1]) {
                    std::cerr << "Error: Invalid value for -k. Please provide ascending integers between 2 and w/2.\n";
                    print_help();
                    std::exit(EXIT_FAILURE);
                }
            }
        } else if (arg == "-s") {
            if (i + 1 < argc) {
                args.dna = argv[++i];
            } else {
                std::cerr << "Error: -s option requires an argument.\n";
                print_help();
                std::exit(EXIT_FAILURE);
            }
            // concatenate all remaining arguments if they dont start with '-'
            while (i + 1 < argc && argv[i + 1][0] != '-') {
                args.dna += argv[++i];
            }
        } else if (arg == "-v") {
            args.verbose = true;
        } else if (arg == "-h") {
            print_help();
            std::exit(EXIT_SUCCESS);
        } else {
            std::cerr << "Error: Unknown option '" << arg << "'.\n";
            print_help();
            std::exit(EXIT_FAILURE);
        }
    }
}

int main(int argc, char **argv) {
    cmd_arguments args;
    parse_arguments(argc, argv, args);
    std::vector<float> result = run_program(args.w, args.k_values, args.dna);
    for (size_t i = 0; i < result.size(); i++)
    {
        if (!args.verbose){
            std::cout << result[i] << std::endl;
        }
        else{
            std::cout << args.dna[i] << '\t' << i << '\t' << result[i] << std::endl;
        }
    }

    return 0;
}