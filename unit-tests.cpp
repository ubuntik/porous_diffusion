#include "lalgebra.h"

int main()
{
	matrix a(2);
	matrix b(2);
	matrix c(2);

	a(0,0) = 1; a(1,0) = 2; a(0,1) = 3; a(1,1) = 4;
	b(0,0) = 8; b(1,0) = 7; b(0,1) = 6; b(1,1) = 5;

	std::cout << "*|1 2|*|8 7| = |20 17|************" << std::endl;
	std::cout << "*|3 4|*|6 5| = |48 41|************" << std::endl;

	c = a * b;
	c.print();
	std::cout << std::endl;

	vector f(2);
	f(0) = 5; f(1) = 6;

	std::cout << "*|1 2|*|5| = |17|*****************" << std::endl;
	std::cout << "*|3 4|*|6| = |39|*****************" << std::endl;

	vector d(2);
	d = a * f;
	d.print();
	std::cout << std::endl;

	std::cout << "*2*|1 2| = |2 4|******************" << std::endl;
	std::cout << "*2*|3 4| = |6 8|******************" << std::endl;

	(a * 2).print();
	std::cout << std::endl;

	std::cout << "*|1 2|+|8 7| = |9 9|**************" << std::endl;
	std::cout << "*|3 4|+|6 5| = |9 9|**************" << std::endl;

	c = a + b;
	c.print();
	std::cout << std::endl;

	std::cout << "*|1 2|-|8 7| = |-7 -5|************" << std::endl;
	std::cout << "*|3 4|-|6 5| = |-3 -1|************" << std::endl;

	c = a - b;
	c.print();
	std::cout << std::endl;

	std::cout << "*|5|+2*|17| = |39|****************" << std::endl;
	std::cout << "*|6|+2*|39| = |84|****************" << std::endl;

	(f + (d * 2)).print();
	std::cout << std::endl;

	std::cout << "*2*|5|-|17| = |-7|****************" << std::endl;
	std::cout << "*2*|6|-|39| = |-27|***************" << std::endl;

	((f * 2) - d).print();
	std::cout << std::endl;

	matrix x(3);
	x(0,0) = 1; x(1,0) = 2; x(2,0) = 3;
	x(0,1) = 3; x(1,1) = 2; x(2,1) = 1;
	x(0,2) = 2; x(1,2) = 1; x(2,2) = 3;

	std::cout << "*|1 2 3| = coFactor = |5 -7 -1|***" << std::endl;
	std::cout << "*|3 2 1| = coFactor = |-3 -3 3|***" << std::endl;
	std::cout << "*|2 1 3| = coFactor = |-4 8 -4|***" << std::endl;

	x.coFactor().print();
	std::cout << std::endl;

	std::cout << "*|1 2 3| = transpose = |1 3 2|****" << std::endl;
	std::cout << "*|3 2 1| = transpose = |2 2 1|****" << std::endl;
	std::cout << "*|2 1 3| = transpose = |3 1 3|****" << std::endl;

	x.transpose().print();
	std::cout << std::endl;

	std::cout << "*|1 2 3| = inverse = |-0.41(6) 0.25 0.(3)|****" << std::endl;
	std::cout << "*|3 2 1| = inverse = |0.58(3) 0.25 -0.(6)|****" << std::endl;
	std::cout << "*|2 1 3| = inverse = |0.08(3) -0.25 0.(3)|****" << std::endl;

	x.inverse().print();
	std::cout << std::endl;

	return 0;
}
