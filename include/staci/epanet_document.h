#ifndef STACI_EPANET_DOCUMENT_H
#define STACI_EPANET_DOCUMENT_H

#include <cstddef>
#include <string>
#include <vector>

struct EpanetDocumentLine {
    std::size_t line_number = 0;
    std::string raw;
    std::string section;
    std::string content;
    std::string inline_comment;
    std::string record_id;
    std::string object_type;
    bool section_header = false;
    bool data_record = false;
};

// Lossless EPANET INP representation. raw_text() is retained byte-for-byte for
// export, while lines() provides section, association, comment, and ordering
// information to STACI features that need structured access.
class EpanetDocument {
public:
    EpanetDocument() = default;

    static EpanetDocument read(const std::string &filename);
    void write(const std::string &filename) const;

    const std::string &raw_text() const { return raw_text_; }
    const std::vector<EpanetDocumentLine> &lines() const { return lines_; }
    bool empty() const { return raw_text_.empty() && lines_.empty(); }

    std::vector<const EpanetDocumentLine *> section_lines(
        const std::string &section) const;

private:
    std::string raw_text_;
    std::vector<EpanetDocumentLine> lines_;

    void parse_lines();
};

#endif
