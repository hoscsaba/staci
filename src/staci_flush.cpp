#include "flushing.h"
#include "diagnostics.h"
int main(int argc, char **argv) {
    return diagnostics::run("staci_flush", argc, argv, flushing::run);
}
