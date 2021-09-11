#include <iostream>
#include <random>

int main()
{
    int rand_num;
    for (int i = 0; i < 10; i++)
    {
        rand_num = (rand() % 4);
        std::cout << rand_num << std::endl;
    }

    return 0;
}