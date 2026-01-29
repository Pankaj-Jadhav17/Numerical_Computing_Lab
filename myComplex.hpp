
class myComplex {
public:
    float a, b, c;
    myComplex();
    myComplex(float x, float y);
    myComplex(float x, float y, float z);
    myComplex add(myComplex obj);
    myComplex subtract(myComplex obj);
    myComplex multiply(myComplex obj);
    myComplex myComplexConjugate(myComplex obj);
    float myComplexNorm(myComplex obj);
    double myComplexDivision(myComplex obj);

};
   