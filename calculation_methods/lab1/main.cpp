#include <iostream>
#include "Matrix.h"

int main() {
    Matrix matrix(5, 1000.23123);
    matrix.make_ones_on_main_diag();
    matrix.set_console_text_colour(YELLOW);
    std::cout << matrix;
    return 0;
}
