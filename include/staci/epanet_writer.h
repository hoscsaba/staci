#ifndef STACI_EPANET_WRITER_H
#define STACI_EPANET_WRITER_H

#include <string>
#include <vector>

class Agelem;
class Csomopont;
class EpanetDocument;

class EpanetWriter {
public:
    static void write(const std::string &filename,
                      const std::vector<Csomopont *> &nodes,
                      const std::vector<Agelem *> &edges,
                      const std::string &friction_model);

    static void copy(const std::string &input_filename,
                     const std::string &output_filename);

    static void write_document(const std::string &output_filename,
                               const EpanetDocument &document);

    static void write_modified_copy(const std::string &input_filename,
                                    const std::string &output_filename,
                                    const std::string &element_id,
                                    const std::string &property,
                                    double value_si);
};

#endif
