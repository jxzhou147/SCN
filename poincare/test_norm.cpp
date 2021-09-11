#include <iostream>
#include <random>
using namespace std;

int main()
{
    normal_distribution<double> norm{1, 0.05};
    default_random_engine rng;
    for (size_t i = 0; i < 10; i++)
        cout << norm(rng) << '\t';
    cout << endl;

    for (size_t i = 0; i < 10; i++)
        cout << norm(rng) << '\t';
    cout << endl;

    return 0;
}