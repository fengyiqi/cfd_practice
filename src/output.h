#ifndef OUTPUT_H
#define OUTPUT_H

#include <string>
#include "block.h"

class OutputWriter {
    std::string file_name_;
    const Block& block_;
public:
    OutputWriter() = delete;
    OutputWriter(const Block& block, std::string file_name);
    void WriteCSV();
};

#endif
