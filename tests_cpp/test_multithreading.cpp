#include<iostream>
#include<thread>



void myFunc()
{
    std::cout << "So, we're in!" << std::endl;
}

int main()
{
    //myFunc();
    std::thread t1(myFunc); // t1 starts running

    //t1.join(); //main thread waits for t1 to finish           
    //t1.detach(); // t1 will run freely (deamon process)
    
    //!! can only detach/join the specific thread only ones
    if (t1.joinable())
    {
        t1.join(); //crash
    }

    return 0;
}



















