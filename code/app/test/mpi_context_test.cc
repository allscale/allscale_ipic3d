#include <gtest/gtest.h>

#include "ipic3d/app/mpi_context.h"

namespace ipic3d {

	TEST(MPI_Context, PartitionTest) {

		using Size = allscale::utils::Vector<int64_t,3>;

		// create partitions for various number of processes
		for(int s=1; s<3000; s++) {	// the group size
			// the grid size
			for(auto gs : { Size(4), Size(16), Size(27), Size(32), Size(64) }) {

				// partition the space
				auto dist = partitioning::partitionSpace(s,gs);

				// check size
				EXPECT_EQ(dist.size(),gs.x*gs.y*gs.z);

				// check that all entries are valid
				std::vector<int> count(s);
				for(const auto& cur : dist) {
					EXPECT_LE(0,cur);
					EXPECT_LT(cur,s);
					count[cur]++;
				}

				// compute min/max
				int min = std::numeric_limits<int>::max();
				int max = std::numeric_limits<int>::min();
				int sum = 0;
				for(const auto& cur : count) {
					min = (min < cur) ? min : cur;
					max = (max > cur) ? max : cur;
					sum += cur;
				}

//				std::cout << "Load-Inbalance: " << s << " / " << gs << ": " << min << " vs. " << max << " - " << (sum/s) << "\n";
				EXPECT_LE(max-min,1);

			}
		}

	}

} // end namespace ipic3d
