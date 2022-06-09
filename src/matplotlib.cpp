#include <iostream>
#include <Eigen/Dense>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

int main() {
    int n = 5000;
    Eigen::VectorXd x(n), y(n), z(n), w = Eigen::VectorXd::Ones(n);
    for (int i = 0; i < n; ++i) {
    double value = (1.0 + i) / n;
    x(i) = value;
    y(i) = value * value;
    z(i) = value * value * value;
    }
    plt::figure();
    plt::plot(x, y, {{"label", "$x^2$"}}); // f(x) = x^2
    plt::plot(x, w, "r--",{{"label", "$1$"}}); // f(x) = 1, red dashed line
    plt::plot(x, z, "g:", {{"label", "$x^3$"}}); // f(x) = x^3, green dots + label
    plt::legend();
    plt::xlim(-0.2, 1.2);
    plt::title("中文Some functions of HipY$_{set}$"); // add a title
    plt::savefig("../matplotlib.png");
    plt::show();
}