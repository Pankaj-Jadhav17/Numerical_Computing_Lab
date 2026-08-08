#include<iostream>
#include"myComplex.hpp"


int main() {
    myComplex c1; // Default Constructor
    myComplex c2(3.0f, 4.0f); // Parameterized Constructor with two parameters
    myComplex c3(1.0f, 2.0f, 5.0f); // Parameterized Constructor with three parameters
    float cnorm;
    cnorm = c2.myComplexNorm(c2); // Norm Function
    std::cout << "Norm of c2: " << cnorm << std::endl;
    double cdiv;
    cdiv = c2.myComplexDivision(c3); // Division Function
    std::cout << "Result of Division: " << cdiv<< " + " << cdiv << "i" << std::endl;

    myComplex c4;
    c4 = c2.add(c3); // Addition Function
    std::cout<< "Result of Addition: " << c4.a << " + " << c4.b << "i" << std::endl;

    myComplex c5;
    c5 = c3.subtract(c2); // Subtraction Function
    std::cout << "Result of Subtraction: " << c5.a << " + " << c5.b << "i" << std::endl;

    myComplex c6;
    c6 = c2.multiply(c3); // Multiplication Function
    std::cout << "Result of Multiplication: " << c6.a << " + " <<c6.b << "i" << std::endl;

    myComplex c7;
    c7 = c3.myComplexConjugate(c2); // Conjugate Function
    std::cout << "Result of Conjugate: " << c7.a << " + " << c7.b << "i" << std::endl; 
    
    
    return 0;
}
