#include "optimalRMQ.hpp"
#include <iostream>
#include <ctime>
using namespace std;

int main(int argc, char* argv[]) {

	srand((unsigned)time(NULL));
	const DTidx n = 16;               //DTidx --  unsigned long
	DT* a = new DT[n];           //DT   --  int

	//for (DTidx i = 0; i < n; i++) {
	//	a[i] = (1 + rand()) % 1000;         //   [0,999]
	//////	if (i < 100) cout << a[i] << ",";
	//}

	a[0] = 5;
	a[1] = 3;
	a[2] = 4;
	a[3] = 3;
	a[4] = 4;
	a[5] = 5;
	a[6] = 1;
	a[7] = 3;
	a[8] = 2;
	a[9] = 4;
	a[10] = 2;
	a[11] = 5;
	a[12] = 3;
	a[13] = 5;
	a[14] = 5;
	a[15] = 4;

	cout << endl;
	RMQ_optimal RMQ(a, n);
	cout << "structure occupies " << RMQ.getSize() << " bytes of memory.\n";

	//DTidx nr_tests = 100000;
	//cout << "now testing " << nr_tests << " random RMQs...\n";
 
	DTidx i = 6;
	DTidx j = 10;

	clock_t start = clock();

	//for (DTidx x = 0; x < nr_tests; x++) {
	//	i = (DTidx)(n * ((float)rand() / (RAND_MAX + 1.0)));    // n * [ 0, 1)  =  [ 0, n)
	//	//    j = i+2; if (j >= n) j = n-1; // small queries
	//	j = (DTidx)(n * ((float)rand() / (RAND_MAX + 1.0)));
	//	//		if (j >= n) j = n-1;
	//	DTidx min = RMQ.query(i, j);

	//	//uncomment this if you want to check results (naively) for correctness:
	//	//for (DTidx k = i; k <= j; k++) {
	//	//	if ((a[k] < a[min]) || (a[k]==a[min] && k < min)) {    
	//	             //论文中提到，多个最值出现时，返回的最值位置是leftmost
	//	// 		cout << "ERROR: " << i << "," << j << endl;
	//	//		cout << i << "," << j << endl;
	//	// 		exit(-1);
	//	//  }
	//	//}
	//}

	DTidx max = RMQ.query(i, j);
	cout << "数组范围[ " << i << " , " << j << "]中，最大值max ： "<<max<<"， value为："<< a[max] << endl;

	clock_t end = clock();

	cout << "time : "<<end - start<<"  ms "<< "\n";
	cout << "...done. Bye!\n";

	system("pause");
	return 0;
}
