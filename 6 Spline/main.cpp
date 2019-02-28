 #include "pch.h"
#include "spline.h"


int main() {
	Spline S;
	S.load();
	S.calcSpline();
	S.output();
}
