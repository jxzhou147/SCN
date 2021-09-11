#include <iostream>

using namespace std;

void matrix(double (*a)[3])
{
    for (size_t i = 0; i < 3; i++)
    {
        for (size_t j = 0; j < 3; j++)
            a[i][j]++;
    }
}

void array(double p[3])
{
    for (size_t i = 0; i < 3; i++)
        p[i]++;
}

int main()
{
    double a[3][3] = {0};
    double p[3] = {0};
    for (size_t i = 0; i < 3; i++)
    {
        /*for (size_t j = 0; j < 3; j++)
        {
            cout << a[i][j] << '\t';
        }
        cout << endl;*/

        cout << p[i] << endl;
    }

    //matrix(a);
    array(p);

    for (size_t i = 0; i < 3; i++)
    { /*
        for (size_t j = 0; j < 3; j++)
        {
            cout << a[i][j] << '\t';
        }
        cout << endl;*/

        cout << p[i] << endl;
    }

    return 0;
}