// functions for calculating linear indices.
#ifndef MAIN
#define INDEX
#include "myheader.cpp"
#endif

namespace index {
int index2(int i1,int i2,int N1,int N2){
    return i2 + i1*N2;
}
int index3(int i1,int i2,int i3,int N1,int N2, int N3){
    return i3 + (i2 + i1*N2)*N3;
}
int index4(int i1,int i2,int i3,int i4,int N1,int N2, int N3, int N4){
    return i4 + (i3 + (i2 + i1*N2)*N3)*N4;
}

int index8(int i1,int i2,int i3,int i4,int i5,int i6,int i7,int i8,int N1,int N2, int N3, int N4,int N5,int N6, int N7, int N8){
    return i8 + (i7 + (i6 + (i5 + (i4 + (i3 + (i2 + i1*N2)*N3)*N4)*N5)*N6)*N7)*N8;
}
}