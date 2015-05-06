#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

int main()
{
    std::ifstream in("inputP.txt");
    std::string line;

    float v[2][3];
    int i = 0, k = 0;

    while (std::getline(in, line))
    {
        float value;
        int k = 0;
        std::stringstream ss(line);

        while (ss >> value)
        {
            v[i][k] = value;
            ++k;
        }
        ++i;
    }
    for(int i=0;i<2;i++)
    {
    	for(int j=0;j<3;j++)
    	{
    		std::cout << v[i][j] << " ";
		}
		std::cout << endl;
	}
}
