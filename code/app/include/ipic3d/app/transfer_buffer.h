#pragma once

#include <array>
#include <vector>

#include "allscale/api/user/data/grid.h"

#include "ipic3d/app/particle.h"
#include "ipic3d/app/universe_properties.h"

namespace ipic3d {

	/**
	 * A class modeling the transfer direction of particles.
	 */
	class TransferDirection {

		friend class TransferBuffers;

		unsigned direction;

	public:

		// The type used to specify directions (int to make writing loops easier)
		using Direction = int;

		// some constants to specify directions in each dimension
		constexpr static Direction Predecessor = 0;
		constexpr static Direction Center      = 1;
		constexpr static Direction Successor   = 2;

		/**
		 * Creates a new transfer direction instance pointing in the specified direction.
		 */
		TransferDirection(Direction x, Direction y, Direction z) : direction((x << 4) | (y << 2) | z) {
			assert_true(0 <= x && x < 3) << x;
			assert_true(0 <= y && y < 3) << y;
			assert_true(0 <= z && z < 3) << z;
		}

		TransferDirection inverse() const {
			return TransferDirection(2 - getDirectionX(), 2 - getDirectionY(), 2 - getDirectionZ());
		}

		Direction getDirectionX() const {
			return (direction >> 4) & 0x3;
		}

		Direction getDirectionY() const {
			return (direction >> 2) & 0x3;
		}

		Direction getDirectionZ() const {
			return direction & 0x3;
		}

		int getDeltaX() const {
			return getDirectionX() - 1;
		}

		int getDeltaY() const {
			return getDirectionY() - 1;
		}

		int getDeltaZ() const {
			return getDirectionZ() - 1;
		}

		Vector3<allscale::api::user::data::coordinate_type> getDelta() const {
			return { getDeltaX(), getDeltaY(), getDeltaZ() };
		}

		coordinate_type step(const coordinate_type& position, const coordinate_type& boundaries) const {
			return ((position + getDelta()) + boundaries) % boundaries;
		}

		template<typename Body>
		static void forEach(const Body& body) {
			for(int x = 0; x < 3; ++x) {
				for(int y = 0; y < 3; ++y) {
					for(int z = 0; z < 3; ++z) {
						body(TransferDirection(x, y, z));
					}
				}
			}
		}

	};

	/**
	 * A class organizing particle-transfer buffers within a simulation.
	 * Internally, this class maintains particle transfer lists connecting
	 * each pair of adjacent cells in a given 3D grid. In each direction,
	 * two distinct buffers are maintained.
	 */
	class TransferBuffers {

		using particle_list = std::vector<Particle>;

		using buffer_grid = allscale::api::user::data::Grid<particle_list,3>;

		using d1 = std::array<buffer_grid,3>;
		using d2 = std::array<d1,3>;
		using d3 = std::array<d2,3>;

		d3 buffers;

	public:

		using grid_size_t = typename buffer_grid::coordinate_type;
		using grid_pos_t = grid_size_t;

		/**
		 * Creates a transfere buffer for a grid of the given size.
		 */
		TransferBuffers(const grid_size_t& size)
			: buffers(d3{
				d2{
					d1{buffer_grid(size),buffer_grid(size),buffer_grid(size)},
					d1{buffer_grid(size),buffer_grid(size),buffer_grid(size)},
					d1{buffer_grid(size),buffer_grid(size),buffer_grid(size)}
				},
				d2{
					d1{buffer_grid(size),buffer_grid(size),buffer_grid(size)},
					d1{buffer_grid(size),buffer_grid(size),buffer_grid(size)},
					d1{buffer_grid(size),buffer_grid(size),buffer_grid(size)}
				},
				d2{
					d1{buffer_grid(size),buffer_grid(size),buffer_grid(size)},
					d1{buffer_grid(size),buffer_grid(size),buffer_grid(size)},
					d1{buffer_grid(size),buffer_grid(size),buffer_grid(size)}
				}
			}) {}

		/**
		 * Obtains a transfer buffer for a source cell at the given position targeting
		 * the given direction.
		 */
		particle_list& getBuffer(const grid_pos_t& src, const TransferDirection& dir) {
			return buffers[dir.getDirectionX()][dir.getDirectionY()][dir.getDirectionZ()][src];
		}

	};

} // end namespace ipic3d

