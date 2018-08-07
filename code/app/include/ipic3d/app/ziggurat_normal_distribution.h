#include <math.h>

namespace ipic3d {

	/**
	 * An approximated normal distribution value generator.
	 */
	class ziggurat_normal_distribution {

		uint64_t jz, jsr = 123456789;

		int64_t hz;
		uint64_t iz, kn[128]/*, ke[256]*/;
		float wn[128], fn[128]/*, we[256], fe[256]*/;

		uint64_t shr3() {
			jz=jsr;
			jsr^=(jsr<<13);
			jsr^=(jsr>>17);
			jsr^=(jsr<<5);
			return jz+jsr;
		}

		float uni() {
			return (float)(.5 + (signed) shr3()*.2328306e-9);
		}

		float rnor() {
			hz = (int64_t)shr3();
			iz = hz & 127;
			return (fabs(hz) < kn[iz]) ? hz * wn[iz] : nfix();
		}

	public:

		explicit ziggurat_normal_distribution(uint64_t seed = 0) {
			zigset(seed);
		}

		float operator()() {
			return nfix();
		}

	private:

		float nfix() {
			int hz;
			uint32_t iz;
			const float r = 3.442620f;
			float value;
			float x;
			float y;

			hz = (int) shr3();
			iz = (hz & 127);

			if (fabs(hz) < kn[iz]) {
				value = (float) (hz) * wn[iz];
			} else {
				for (;;) {
					if (iz == 0) {
						for (;;) {
							x = (float)(-0.2904764 * log(uni()));
							y = -log(uni());
							if (x * x <= y + y) {
								break;
							}
						}

						if (hz <= 0) {
							value = -r - x;
						} else {
							value = +r + x;
						}
						break;
					}

					x = (float) (hz) * wn[iz];

					if (fn[iz] + uni() * (fn[iz - 1] - fn[iz]) < exp(-0.5 * x * x)) {
						value = x;
						break;
					}

					hz = (int) shr3();
					iz = (hz & 127);

					if (fabs(hz) < kn[iz]) {
						value = (float) (hz) * wn[iz];
						break;
					}
				}
			}

			return value;

		}

		/*--------This procedure sets the seed and creates the tables------*/

		void zigset(uint64_t jsrseed) {
			const double m1 = 2147483648.0;
			double dn = 3.442619855899, tn = dn, vn = 9.91256303526217e-3, q;
			int i;
			jsr ^= jsrseed;

			/* Set up tables for RNOR */
			q = vn / exp(-.5 * dn * dn);
			kn[0] = (uint64_t)((dn / q) * m1);
			kn[1] = 0;

			wn[0] = (float)(q / m1);
			wn[127] = (float)(dn / m1);

			fn[0] = 1.;
			fn[127] = (float)(exp(-.5 * dn * dn));

			for (i = 126; i >= 1; i--) {
				dn = sqrt(-2. * log(vn / dn + exp(-.5 * dn * dn)));
				kn[i + 1] = (uint64_t)((dn / tn) * m1);
				tn = dn;
				fn[i] = (float)(exp(-.5 * dn * dn));
				wn[i] = (float)(dn / m1);
			}
		}

	};

} // end namespace ipic3d
