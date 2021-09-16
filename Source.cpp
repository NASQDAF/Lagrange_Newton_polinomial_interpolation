#include <vector>
#include <iostream>
#include <iomanip>

size_t n = 50, another_n = 10; 
constexpr double   a = -2, b = 2;
std::vector<std::pair<double, double>> pointsData;
inline double f(double x) { return -x * (1 - exp(cos(x)) + 2 * x * cos(x)); } //y = f(x)


inline float h(size_t k = n) { return (b - a) / k; }
inline double offset(size_t i, double k = h()) { return a + i * k; }

void PointsWeKnow_FromSteps() {
    for (size_t i = 0; i < n; i++)
        pointsData.push_back({ offset(i), f(offset(i)) });
}

double Lagrange(double x) {
    double l=1,approx=0;
    for (size_t i = 0; i < pointsData.size(); l = 1.0f, i++) {
        for (size_t j = 0; j < pointsData.size(); j++)
            if (i != j) 
                l *= (x - pointsData[j].first) / (pointsData[i].first - pointsData[j].first);
        approx += l * pointsData[i].second;
    }
    return approx;
}

int main() {
    PointsWeKnow_FromSteps();
    
    std::cout.setf(std::ios::left, std::ios::adjustfield);
    std::cout << "\`-. |   |  .-.  Function: f(x) = -x (cos(x)2x - exp(cos(x)) + 1)\n\___`-.__|,`_|_  Count of Data Points: " << pointsData.size() << "\n    | `-`       Count of prediction points: " << another_n << "\n\n\tLagrange method:\n\n" << std::setw(7) << "" << "x:" << std::setw(13) << "" << "y = f(x):" << std::endl;
    for (size_t i = 0; i < pointsData.size(); i++)
        std::cout << "[" << std::setw(2) << i + 1 << "]:  " << std::setw(15) << pointsData[i].first << std::setw(15) << pointsData[i].second << std::endl;

    n = another_n;
    std::cout << "\n\n" << std::setw(8) << "" << "f(x): " << std::setw(5) << "" << "Ln(x):" << std::setw(7) << "" << "Residual:" << std::endl;
    double tempLag, tempFunc;
    for (size_t i = 0; i < n; i++) {
        tempLag = Lagrange(offset(i));
        tempFunc = f(offset(i));
        std::cout << "[" << std::setw(2) << i + 1 << "]:  " << std::setw(12) << tempFunc << std::setw(12) << tempLag << std::setw(15) << abs(tempFunc) - abs(tempLag) << std::endl;
    }
    system("pause");
}