#include <vector>
#include <iostream>
#include <iomanip>
using std::vector;
#pragma region Data

size_t n = 25, another_n = 10;
constexpr double   a = -2, b = 2;
vector<std::pair<double, double>> pointsData;
vector<double> NewtonFwPoints;
inline double f(double x) { return -x * (1 - exp(cos(x)) + 2 * x * cos(x)); } //y = f(x)


class Offset{ // [x] steps frequency
public:
    inline double f(size_t i) {
        switch (tumb) {
        case 0:
            return a + i * ((b - a) / n);
        case 1:
            return 0.5 * (b + a - (b - a) * cos(((2 * i + 1) * pi) / (2 * n + 2)));
        }
    }
    void operator=(size_t i) { tumb = i; }
private:
    const double pi = acos(-1);
    static int tumb;
}offset;
int Offset::tumb = 0;


#pragma endregion

template<typename T>
void PointsWeKnow_FromSteps(T f) {
    pointsData.clear();
    for (size_t i = 0; i < n; i++)
        pointsData.push_back({ offset.f(i), f(offset.f(i)) });
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

#pragma region Newton
void insertFirstY(vector<double>& a) {
    for (auto i : pointsData)
        a.push_back(i.second);
}
void Newton_forward_points() {
    vector<vector<double>> vvd(n);
    NewtonFwPoints.push_back(pointsData[0].second);
    insertFirstY(vvd[0]);
    for (size_t i = 1; i < n; ++i) {
        for (size_t j = 0; j < n - i; ++j)
            vvd[i].push_back((vvd[i - 1][j + 1] - vvd[i - 1][j]) / (offset.f(j + i) - offset.f(j)));
        NewtonFwPoints.push_back(vvd[i][0]);
    }
}
double Newton_forward(double x) {
    double d = 0, f = NewtonFwPoints[0];
    for (size_t i = 1; i < n; ++i) {
        d = NewtonFwPoints[i];
        for (size_t j = 0; j < i; ++j)
            d *= (x - offset.f(j));
        f += d;
    }
    return f;
}
#pragma endregion

template<typename T>
inline void show(T Func) {
    double tempInterpolation, tempFunc;
    for (size_t i = 0; i < n; i++) {
        tempInterpolation = Func(offset.f(i));
        tempFunc = f(offset.f(i));
        std::cout << "[" << std::setw(2) << i + 1 << "]:  " 
            << std::setw(12) << tempFunc 
            << std::setw(12) << tempInterpolation 
            << std::setw(15) << abs(tempFunc) - abs(tempInterpolation) 
            << std::endl;
    }
}

inline void ioman(bool r) {
    if(r)std::cout.setf(std::ios::right, std::ios::adjustfield); 
    else std::cout.setf(std::ios::left, std::ios::adjustfield);
}

int main() {
    using namespace std;
    PointsWeKnow_FromSteps(f);

    cout << "\`-. |   |  .-.  Function: f(x) = -x (cos(x)2x - exp(cos(x)) + 1)\n" 
        "\___`-.__|,`_|_  Count of Data Points: " << pointsData.size() 
        << "\n    | `-`       Count of prediction points: " << another_n 
        << "\n\n"<< setw(10) << "x" << setw(15) << "f(x):\n" << endl;
    ioman(0);
    for (size_t i = 0; i < pointsData.size(); i++)
        cout << "[" << setw(2) << i + 1 << "]:  " 
        << setw(10) << pointsData[i].first 
        << pointsData[i].second 
        << endl;

    n = another_n;

    ioman(1);
    cout << "\n\n\tLagrange method:\n\n" 
        << setw(15) << "f(x): "
        << setw(11) << "Ln(x):"
        << setw(15) << "Residual:"
        << std::endl;
    ioman(0);
    show(Lagrange);

    ioman(1);
    cout << "\n\n\tNewton method:\n\n" 
        << setw(15) << "f(x): "
        << setw(11) << "Nw(x):" 
        << setw(15) << "Residual:" 
        << std::endl;
    ioman(0);
    Newton_forward_points();
    show(Newton_forward);

    offset = 1; //Algorithm step change
    PointsWeKnow_FromSteps(f); //Recalc
    ioman(1);
    cout << "\n\n\tLagrange method:\n\n" 
        << setw(15) << "f(x): "
        << setw(11) << "Ln(x):"
        << setw(15) << "Residual:"
        << std::endl;
    ioman(0);
	n = 15; //Changing amount of points again
	show(Lagrange); 
    system("pause");
}