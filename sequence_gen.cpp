#include <iostream>
#include <fstream>
#include <string>
#include <random>

void generateRandomDNA(const std::string& filename, int64_t length) {
    // DNAの文字（A, C, G, T）
    const char dnaChars[] = {'A', 'C', 'G', 'T'};
    const int numChars = sizeof(dnaChars) / sizeof(dnaChars[0]);

    // ランダム数生成器の初期化
    std::random_device rd;  // シードの取得
    std::mt19937 gen(rd()); // メルセンヌ・ツイスタ法による乱数生成器
    std::uniform_int_distribution<> dis(0, numChars - 1);

    // ファイルのオープン
    std::ofstream file(filename);
    if (!file) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // 1行目にヘッダーを書き込む
    file << "> TEST DATA" << std::endl;

    // ランダムDNA配列を生成して書き込む
    for (int64_t i = 0; i < length; ++i) {
        char randomChar = dnaChars[dis(gen)];
        file << randomChar;
    }
    file << std::endl; // 最後に改行を追加

    // ファイルを閉じる
    file.close();
}

int main(int argc, char *argv[]) {
    std::string filename;
    int64_t x;

    if (argc >= 2) {
        filename = argv[1];
        if (argc >= 3) {
            x = atoi(argv[2]);
        }
    } else {
        // ファイル名と文字列長を入力
        std::cout << "Enter the filename (e.g., filename.fa): ";
        std::cin >> filename;
        std::cout << "Enter the length of the random DNA sequence: ";
        std::cin >> x;
    }

    // ランダムDNA配列を生成
    generateRandomDNA(filename, x);

    std::cout << "File " << filename << " has been created with " << x << " random DNA characters." << std::endl;

    return 0;
}
