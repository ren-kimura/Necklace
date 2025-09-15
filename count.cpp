#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>

using INT = uint64_t;

INT count_bases(const std::string& filename) {
    std::ifstream file(filename, std::ios::in | std::ios::binary);
    if (!file) {
        std::cerr << "Error: Cannot open file " << filename << "\n";
        return -1;
    }

    const size_t BUF_SIZE = 64 * 1024 * 1024;
    std::vector<char> buffer(BUF_SIZE);
    file.rdbuf()->pubsetbuf(buffer.data(), BUF_SIZE);

    std::string line, read;
    INT count = 0;
    bool is_new_read = true;

    while (std::getline(file, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') {
            if (!read.empty()) {
                count += read.size();
                read.clear();
            }
            is_new_read = true;
        } else {
            if (is_new_read) {
                read = line;
                is_new_read = false;
            } else {
                read += line;
            }
        }
    }
    if (!read.empty()) count += read.size();

    file.close();
    return count;
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <FASTA file>" << "\n";
        return 1;
    }

    INT base_count = count_bases(argv[1]);
    if (base_count >= 0) {
        std::cout << "Total base count: " << base_count << "\n";
    }

    return 0;
}