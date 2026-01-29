 #include<iostream>
 #include<cmath>
 #include"myComplex.hpp"
 
 
 
    myComplex::myComplex()
    {
        std::cout << "Default Constructor Called" << std::endl;
    }
    myComplex::myComplex(float x, float y)
    {
        a = x;
        b = y;
        std::cout << "Parameterized Constructor with two parameters Called" << std::endl;
    }
    myComplex::myComplex(float x, float y, float z)
    {
        a = x;
        b = y;
        c = z;
        std::cout << "Parameterized Constructor with three parameters Called" << std::endl;
    }
    myComplex myComplex::add(myComplex obj)
    {
        myComplex temp;
        temp.a = a + obj.a;
        temp.b = b + obj.b;
        return temp;
    }
    myComplex myComplex::subtract(myComplex obj)
    {
        myComplex temp;
        temp.a = a - obj.a;
        temp.b = b - obj.b;
        return temp;
    }
    myComplex myComplex::multiply(myComplex obj)
    {
        myComplex temp;
        temp.a = (a * obj.a- b * obj.b);
        temp.b = (a * obj.b + b * obj.a);
        return temp;
    }
    myComplex myComplex::myComplexConjugate(myComplex obj)
    {
        myComplex temp;
        temp.a = obj.a;
        temp.b = -obj.b;
        return temp;
    }   
    float myComplex::myComplexNorm(myComplex obj)
{
     float temp;
    temp = sqrt(obj.a * obj.a + obj.b * obj.b);
    return temp;
}

double myComplex::myComplexDivision(myComplex obj)
{
     double temp;
    temp = (a * obj.a + b * obj.b) / (obj.a * obj.a + obj.b * obj.b);
    temp = (b * obj.a - a * obj.b) / (obj.a * obj.a + obj.b * obj.b);
    return temp;
}
