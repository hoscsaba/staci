#include "epanet_document.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iterator>
#include <sstream>
#include <stdexcept>

namespace {

std::string trim(const std::string &value) {
    const std::string whitespace = " \t\r\n";
    const std::size_t begin = value.find_first_not_of(whitespace);
    if (begin == std::string::npos)
        return "";
    const std::size_t end = value.find_last_not_of(whitespace);
    return value.substr(begin, end - begin + 1);
}

std::string upper(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    return value;
}

std::vector<std::string> fields(const std::string &value) {
    std::istringstream input(value);
    std::vector<std::string> result;
    std::string field;
    while (input >> field)
        result.push_back(field);
    return result;
}

} // namespace

EpanetDocument EpanetDocument::read(const std::string &filename) {
    std::ifstream input(filename.c_str(), std::ios::binary);
    if (!input)
        throw std::runtime_error("Cannot open EPANET input file: " + filename);

    EpanetDocument document;
    document.raw_text_.assign(std::istreambuf_iterator<char>(input),
                              std::istreambuf_iterator<char>());
    document.parse_lines();
    return document;
}

void EpanetDocument::write(const std::string &filename) const {
    std::ofstream output(filename.c_str(), std::ios::binary | std::ios::trunc);
    if (!output)
        throw std::runtime_error("Cannot create EPANET input file: " + filename);
    output.write(raw_text_.data(), static_cast<std::streamsize>(raw_text_.size()));
    if (!output)
        throw std::runtime_error("Cannot write EPANET input file: " + filename);
}

std::vector<const EpanetDocumentLine *> EpanetDocument::section_lines(
    const std::string &section) const {
    const std::string requested = upper(section);
    std::vector<const EpanetDocumentLine *> result;
    for (const EpanetDocumentLine &line : lines_)
        if (line.section == requested)
            result.push_back(&line);
    return result;
}

void EpanetDocument::parse_lines() {
    lines_.clear();
    std::string section;
    std::size_t position = 0;
    std::size_t line_number = 0;
    while (position < raw_text_.size()) {
        const std::size_t newline = raw_text_.find('\n', position);
        const std::size_t end = newline == std::string::npos ? raw_text_.size() : newline;
        std::string raw = raw_text_.substr(position, end - position);
        ++line_number;

        std::string parsed = raw;
        if (!parsed.empty() && parsed.back() == '\r')
            parsed.pop_back();
        if (line_number == 1 && parsed.size() >= 3 &&
            static_cast<unsigned char>(parsed[0]) == 0xef &&
            static_cast<unsigned char>(parsed[1]) == 0xbb &&
            static_cast<unsigned char>(parsed[2]) == 0xbf)
            parsed.erase(0, 3);

        EpanetDocumentLine line;
        line.line_number = line_number;
        line.raw = raw;
        const std::size_t comment = parsed.find(';');
        if (comment != std::string::npos)
            line.inline_comment = parsed.substr(comment);
        line.content = trim(parsed.substr(0, comment));

        if (!line.content.empty() && line.content.front() == '[' &&
            line.content.back() == ']') {
            section = upper(trim(line.content.substr(1, line.content.size() - 2)));
            line.section = section;
            line.section_header = true;
        } else {
            line.section = section;
            line.data_record = !line.content.empty();
            if (line.data_record) {
                const std::vector<std::string> row = fields(line.content);
                if (!row.empty()) {
                    if (section == "TAGS" && row.size() > 1) {
                        line.object_type = upper(row[0]);
                        line.record_id = row[1];
                    } else {
                        line.record_id = row[0];
                        if (section == "JUNCTIONS" || section == "RESERVOIRS" ||
                            section == "TANKS" || section == "DEMANDS" ||
                            section == "QUALITY" || section == "EMITTERS" ||
                            section == "COORDINATES")
                            line.object_type = "NODE";
                        else if (section == "PIPES" || section == "PUMPS" ||
                                 section == "VALVES" || section == "VERTICES" ||
                                 section == "STATUS")
                            line.object_type = "LINK";
                    }
                }
            }
        }
        lines_.push_back(line);
        if (newline == std::string::npos)
            break;
        position = newline + 1;
    }
}
