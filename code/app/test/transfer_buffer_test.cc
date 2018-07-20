#include <gtest/gtest.h>

#include "ipic3d/app/transfer_buffer.h"

namespace ipic3d {

	TEST(TransferBuffers, Creation) {

		using size_t = TransferBuffers::grid_size_t;

		size_t size(10,12,14);

		TransferBuffers buffers(size);

		std::set<void*> all_buffers;

		// check that all buffers have been properly initialized
		for(int x = 0; x<10; x++) {
			for(int y = 0; y<12; y++) {
				for(int z = 0; z<14; z++) {

					// create the source position
					size_t pos(x,y,z);

					// in each direction
					for(int i=0; i<3; i++) {
						for(int j=0; j<3; j++) {
							for(int k=0; k<3; k++) {

								// create the direction
								TransferDirection dir(i,j,k);

								auto& buffer = buffers.getBuffer(pos,dir);

								// check that the buffer is empty
								EXPECT_TRUE(buffer.empty());

								// check that this is a buffer not seen before
								EXPECT_TRUE(all_buffers.insert(&buffer).second);

							}
						}
					}

				}
			}
		}

		EXPECT_EQ(10*12*14*3*3*3, all_buffers.size());

		// test that particles can be added
		auto& buffer = buffers.getBuffer({5,6,7},TransferDirection(0,1,2));

		EXPECT_TRUE(all_buffers.find(&buffer) != all_buffers.end());

		EXPECT_EQ(0,buffer.size());

		Particle p;
		buffer.push_back(p);

		EXPECT_EQ(1,buffer.size());

		buffer.clear();

		EXPECT_EQ(0,buffer.size());

	}

}
