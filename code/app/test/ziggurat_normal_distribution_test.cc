#include <gtest/gtest.h>

#include <random>

#include "ipic3d/app/ziggurat_normal_distribution.h"


namespace ipic3d {

	TEST(Cell, initCellsUniform) {

		ziggurat_normal_distribution dist;

		int N = 100;
		std::vector<int> hist(N);

		for(int i=0; i<5000; i++) {

			int p = dist() * N/8 + N/2;
			if (0 <= p && p < N) {
				hist[p]++;
			}
		}

		for(int i=0; i<N; i++) {
			std::cout<< "i=" << i << ": ";
			for(int j=0; j<hist[i];j++) {
				std::cout << "*";
			}
			std::cout << "\n";
		}
		std::cout << "\n";
	}

	int N = 10000000;

	TEST(Norm,Bench_ziggurat) {

		ziggurat_normal_distribution dist;
		float sum = 0;
		for(int i=0; i<N; i++) {
			sum += dist();
		}
		std::cout << sum << "\n";
	}

	TEST(Norm,Bench_std_double) {

		std::minstd_rand rd;
		std::normal_distribution<> norm;

		ziggurat_normal_distribution dist;
		float sum = 0;
		for(int i=0; i<N; i++) {
			sum += norm(rd);
		}
		std::cout << sum << "\n";
	}

	TEST(Norm,Bench_std_float) {

		std::minstd_rand rd;
		std::normal_distribution<float> norm;

		ziggurat_normal_distribution dist;
		float sum = 0;
		for(int i=0; i<N; i++) {
			sum += norm(rd);
		}
		std::cout << sum << "\n";
	}
}
