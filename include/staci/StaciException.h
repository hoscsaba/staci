#ifndef STACI_EXCEPTION_H
#define STACI_EXCEPTION_H

#include <stdexcept>
#include <string>

class StaciException : public std::runtime_error
{
public:
    explicit StaciException(const std::string &message)
        : std::runtime_error(message) {}

    std::string getDescription() const { return what(); }
};

#endif
