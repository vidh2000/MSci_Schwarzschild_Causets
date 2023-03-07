#include <iostream>
#include <vector>

using namespace std;

class MyClass {
public:
    vector<int> x, y, z;

    MyClass()
    {
        x = {1,1,1,1,1};
        y = {1,1,1,1,1, 1,1,1,1,1};
        z = {10001, 10001, 10001, 10001, 10001};
    }
};

int main(){
    vector<int> a = {1,1,1,1,1};
    vector<int> b = {1,1,1,1,1, 1,1,1,1,1};
    vector<int> c = {10001, 10001, 10001, 10001, 10001};

    cout<<sizeof(a)<<endl;
    int A = 0;
    for (auto ai : a)
    {A += sizeof(ai);}
    cout<<A<<endl;
    cout<<sizeof(b)<<endl;
    int B = 0;
    for (auto bi : b)
    {B += sizeof(bi);}
    cout<<B<<endl;
    cout<<sizeof(c)<<endl;
    int C = 0;
    for (auto ci : c)
    {C += sizeof(ci);}
    cout<<C<<endl;

    MyClass XYZ;
    cout<<sizeof(XYZ)<<endl;
    cout<<sizeof(XYZ.x)<<endl;
    cout<<sizeof(XYZ.y)<<endl;
    cout<<sizeof(XYZ.z)<<endl;

    return 0;
}