#include "epanet_document.h"

#include <iostream>
#include <string>
#include <vector>

namespace {

bool has_associated_comment(const EpanetDocument &document,
                            const std::string &section,
                            const std::string &record_id,
                            const std::string &comment_text) {
    for (const EpanetDocumentLine *line : document.section_lines(section))
        if (line->record_id == record_id &&
            line->inline_comment.find(comment_text) != std::string::npos)
            return true;
    return false;
}

bool has_record(const EpanetDocument &document,
                const std::string &section,
                const std::string &record_id,
                const std::string &object_type) {
    for (const EpanetDocumentLine *line : document.section_lines(section))
        if (line->data_record && line->record_id == record_id &&
            line->object_type == object_type)
            return true;
    return false;
}

} // namespace

int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: test_epanet_document fixture.inp\n";
        return 2;
    }

    try {
        const EpanetDocument document = EpanetDocument::read(argv[1]);
        if (!has_associated_comment(document, "JUNCTIONS", "J1",
                                    "junction demand comment") ||
            !has_associated_comment(document, "PIPES", "P1",
                                    "pipe association comment")) {
            std::cerr << "Inline comment association was not retained\n";
            return 1;
        }
        if (!has_record(document, "TAGS", "J1", "NODE") ||
            !has_record(document, "TAGS", "P1", "LINK")) {
            std::cerr << "Node/link tag association was not retained\n";
            return 1;
        }
        std::vector<const EpanetDocumentLine *> vertices;
        for (const EpanetDocumentLine *line : document.section_lines("VERTICES"))
            if (line->data_record)
                vertices.push_back(line);
        if (vertices.size() != 2 || vertices[0]->record_id != "P1" ||
            vertices[1]->record_id != "P1" ||
            vertices[0]->line_number >= vertices[1]->line_number) {
            std::cerr << "Vertex association or ordering was not retained\n";
            return 1;
        }
        for (const std::string &section : {
                 "TITLE", "COORDINATES", "LABELS", "BACKDROP",
                 "TIMES", "REPORT", "ENERGY", "OPTIONS"}) {
            if (document.section_lines(section).empty()) {
                std::cerr << "Section was not retained: " << section << '\n';
                return 1;
            }
        }
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }

    return 0;
}
