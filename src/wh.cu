#include "types.cuh"
#include "wh.cuh"
#include "convert.h"
#include "util.cuh"

namespace sr
{
namespace wh
{
	using namespace sr::data;

	const size_t MAXKEP = 5;
	const float64_t TOLKEP = 1E-14;

	struct MVSKernel
	{
		const float64_t* planet_m;
		const uint32_t* planet_id;
		const float64_t mu;
		const f64_3* planet_h0_log;
		const f64_3* planet_r_log;
		const uint32_t planet_n;
		const uint32_t tbsize;
		const float64_t dt;

		const float64_t outer_r;
		const float64_t* planet_rh;

		MVSKernel(const DevicePlanetPhaseSpace& planets, const Dvf64_3& h0_log, const Dvf64& _planet_rh, double _outer_r, uint32_t _tbsize, float64_t _dt) :
			planet_m(planets.m.data().get()),
			planet_id(planets.id.data().get()),
			mu(planets.m[0]),
			planet_h0_log(h0_log.data().get()),
			planet_r_log(planets.r_log.data().get()),
			planet_n(static_cast<uint32_t>(planets.n_alive)),
			tbsize(_tbsize),
			dt(_dt),
			outer_r(_outer_r),
			planet_rh(_planet_rh.data().get())
		{ }

		__host__ __device__
		static void kepeq(double dM, double ecosEo, double esinEo, double* dE, double* sindE, double* cosdE, uint16_t& flags)
		{
			double f, fp, delta;

			*sindE = sin(*dE);
			*cosdE = cos(*dE);

			for (size_t i = 0; i < MAXKEP; i++)
			{
				f = *dE - ecosEo * (*sindE) + esinEo * (1. - *cosdE) - dM;
				fp = 1. - ecosEo * (*cosdE) + esinEo * (*sindE);
				delta = -f / fp;

				*dE += delta;
				*sindE = sin(*dE);
				*cosdE = cos(*dE);

#ifdef CUDA_KEPEQ_CHECK_CONVERGENCE
				if (fabs(delta) < TOLKEP)
				{
					goto done;
				}
#endif
			}

			flags = static_cast<uint16_t>(flags | ((fabs(delta) > TOLKEP) << 3));
#ifdef CUDA_KEPEQ_CHECK_CONVERGENCE
done: ;
#endif
		}

		__host__ __device__
		static void drift(f64_3& r, f64_3& v, uint16_t& flags, double dt, double mu)
		{
			float64_t dist = sqrt(r.lensq());
			float64_t vdotr = v.x * r.x + v.y * r.y + v.z * r.z;

			float64_t energy = v.lensq() * 0.5 - mu / dist;

			flags = static_cast<uint16_t>(flags | ((energy >= 0) << 2));

			float64_t a = -0.5 * mu / energy;
			float64_t n_ = sqrt(mu / (a * a * a));
			float64_t ecosEo = 1.0 - dist / a;
			float64_t esinEo = vdotr / (n_ * a * a);
			float64_t e = sqrt(ecosEo * ecosEo + esinEo * esinEo);

			// subtract off an integer multiple of complete orbits
			float64_t dM = dt * n_ - M_2PI * (int) (dt * n_ / M_2PI);

			// remaining time to advance
			float64_t _dt = dM / n_;

			// call kepler equation solver with initial guess in dE already
			float64_t dE = dM - esinEo + esinEo * cos(dM) + ecosEo * sin(dM);
			float64_t sindE, cosdE;
			kepeq(dM, ecosEo, esinEo, &dE, &sindE, &cosdE, flags);

			float64_t fp = 1.0 - ecosEo * cosdE + esinEo * sindE;
			float64_t f = 1.0 + a * (cosdE - 1.0) / dist;
			float64_t g = _dt + (sindE - dE) / n_;
			float64_t fdot = -n_ * sindE * a / (dist * fp);
			float64_t gdot = 1.0 + (cosdE - 1.0) / fp;

			f64_3 r0 = r;
			r = r0 * f + v * g;
			v = r0 * fdot + v * gdot;
		}

		__host__ __device__
		static void step_forward(
				f64_3& r,
				f64_3& v,
				uint16_t& flags,
				f64_3& a,
				uint32_t& deathtime_index,
				uint32_t _tbsize,
				uint32_t planet_n,
				const f64_3* h0_log,
				const f64_3* r_log,
				const float64_t* m,
				const uint32_t* pl_id,
				const float64_t* rh,
				double outer_radius,
				float64_t dt,
				float64_t mu)
		{
			deathtime_index = 0;

			for (uint32_t step = 0; step < static_cast<uint32_t>(_tbsize); step++)
			{
				if (flags == 0)
				{
					// kick
					// if step = 0, the acceleration is preloaded - this comes from 
					v = v + a * (dt / 2);

					drift(r, v, flags, dt, mu);

					a = h0_log[step];

					// planet 0 is not counted
					for (uint32_t i = 1; i < static_cast<uint32_t>(planet_n); i++)
					{
						f64_3 dr = r - r_log[step * (planet_n - 1) + i - 1];

						float64_t rad = dr.lensq();

						if (rad < rh[i] * rh[i] && flags == 0)
						{
							flags = flags & 0x00FF;
							flags = static_cast<uint16_t>(flags | (pl_id[i] << 8) | 0x0001);
						}

						float64_t inv3 = 1. / (rad * sqrt(rad));
						float64_t fac = m[i] * inv3;

						a -= dr * fac;
					}

					float64_t rad = r.lensq();
					if (rad < rh[0] * rh[0])
					{
						flags = flags & 0x00FF;
						flags = flags | 0x0001;
					}
					if (rad > outer_radius * outer_radius)
					{
						flags = flags | 0x0002;
					}


					v = v + a * (dt / 2);

					deathtime_index = step + 1;
				}
			}
		}

		template<typename Tuple>
		__host__ __device__
		void operator()(Tuple args) const
		{
			uint32_t _tbsize = this->tbsize;
			const f64_3* h0_log = this->planet_h0_log;
			const f64_3* r_log = this->planet_r_log;
			const float64_t* m = this->planet_m;
			const uint32_t* ids = this->planet_id;
			const float64_t* rh = this->planet_rh;

			f64_3 r = thrust::get<0>(thrust::get<0>(args));
			f64_3 v = thrust::get<1>(thrust::get<0>(args));
			uint16_t flags = thrust::get<2>(thrust::get<0>(args));
			uint32_t deathtime_index = 0;
			f64_3 a = thrust::get<1>(args);

			step_forward(r, v, flags, a, deathtime_index, _tbsize,
				planet_n, h0_log, r_log, m, ids, rh, this->outer_r, this->dt, this->mu);

			thrust::get<0>(thrust::get<0>(args)) = r;
			thrust::get<1>(thrust::get<0>(args)) = v;
			thrust::get<2>(thrust::get<0>(args)) = flags;
			thrust::get<3>(thrust::get<0>(args)) = deathtime_index;

			thrust::get<1>(args) = a;
		}
	};

	__global__
	void MVSKernel_(
			f64_3* r,
			f64_3* v,
			uint16_t* flags,
			f64_3* a,
			uint32_t* deathtime_index,
			uint32_t n,
			uint32_t tbsize,
			uint32_t planet_n,
			const f64_3* h0_log,
			const f64_3* r_log,
			const float64_t* m,
			const uint32_t* pl_id,
			const float64_t* rh,
			double outer_r,
			float64_t dt,
			float64_t mu)
	{
		// max timeblock size: 384
		__shared__ f64_3 h0_log_shared[384];
		__shared__ f64_3 r_log_shared[1536];

		for (int i = threadIdx.x; 
				i < tbsize;
				i += blockDim.x)
		{
			h0_log_shared[i] = h0_log[i];
		}

		for (int i = threadIdx.x; 
				i < tbsize * (planet_n - 1);
				i += blockDim.x)
		{
			r_log_shared[i] = r_log[i];
		}

		__syncthreads();

		for (int i = blockIdx.x * blockDim.x + threadIdx.x; 
				i < n; 
				i += blockDim.x * gridDim.x) 
		{
			f64_3 ri = r[i];
			f64_3 vi = v[i];
			uint16_t flagsi = flags[i];
			f64_3 ai = a[i];
			uint32_t deathtime_indexi;
		
			MVSKernel::step_forward(ri, vi, flagsi, ai, deathtime_indexi, tbsize,
					planet_n, h0_log_shared, r_log_shared, m, pl_id, rh, outer_r, dt, mu);

			r[i] = ri;
			v[i] = vi;
			flags[i] = flagsi;
			a[i] = ai;
			deathtime_index[i] = deathtime_indexi;
		}
	}
	

	WHCudaIntegrator::WHCudaIntegrator() { }

	WHCudaIntegrator::WHCudaIntegrator(HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, const Configuration& config)
		: base(pl, pa, config)
	{
		device_h0_log_0 = Dvf64_3(config.tbsize);
		device_h0_log_1 = Dvf64_3(config.tbsize);
		device_particle_a = Dvf64_3(pa.n());

		device_planet_rh = Dvf64(pl.n());

		memcpy_htd(device_planet_rh, base.planet_rh, 0);
		cudaStreamSynchronize(0);
	}

	Dvf64_3& WHCudaIntegrator::device_h0_log(size_t planet_data_id)
	{
		return planet_data_id % 2 ? device_h0_log_1 : device_h0_log_0;
	}

	void WHCudaIntegrator::helio_acc_particles(
			const HostPlanetPhaseSpace& pl,
			HostParticlePhaseSpace& pa,
			size_t begin,
			size_t len,
			float64_t time,
			size_t timestep_index,
			bool old)
	{
		base.helio_acc_particles(pl, pa, begin, len, time, timestep_index, old);
	}

	void WHCudaIntegrator::upload_planet_log_cuda(cudaStream_t stream, size_t planet_data_id)
	{
		memcpy_htd(device_h0_log(planet_data_id), base.planet_h0_log.log, stream);
		cudaStreamSynchronize(stream);
	}

	void WHCudaIntegrator::gather_particles(const std::vector<size_t>& indices, size_t begin, size_t length)
	{
		base.gather_particles(indices, begin, length);
	}

	void WHCudaIntegrator::integrate_planets_timeblock(HostPlanetPhaseSpace& pl, size_t nstep, float64_t t, double dt)
	{
		base.integrate_planets_timeblock(pl, nstep, t, dt);
	}

	void WHCudaIntegrator::load_h0(const HostPlanetPhaseSpace& pl)
	{
		base.load_h0(pl);
	}

	void WHCudaIntegrator::integrate_particles_timeblock(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, size_t begin, size_t length, float64_t t, double dt)
	{
		base.integrate_particles_timeblock(pl, pa, begin, length, t, dt);
	}

	void WHCudaIntegrator::swap_logs()
	{
		base.swap_logs();
	}

	void WHCudaIntegrator::upload_data_cuda(cudaStream_t stream, size_t begin, size_t length)
	{
		memcpy_htd(device_particle_a, base.particle_a, stream, begin, begin, length);
		cudaStreamSynchronize(stream);
	}

	void WHCudaIntegrator::integrate_particles_timeblock_cuda(cudaStream_t stream, size_t planet_data_id, const DevicePlanetPhaseSpace& pl, DeviceParticlePhaseSpace& pa, double dt)
	{
#ifndef CUDA_USE_SHARED_MEM_CACHE
		auto it = thrust::make_zip_iterator(thrust::make_tuple(pa.begin(), device_begin()));
		thrust::for_each(
				thrust::cuda::par.on(stream),
				it,
				it + pa.n_alive,
				MVSKernel(
					pl,
					device_h0_log(planet_data_id),
					device_planet_rh,
					base.outer_radius,
					static_cast<uint32_t>(pl.log_len),
					dt
				)
		);
#else
		cudaDeviceProp prop;
		cudaGetDeviceProperties(&prop, 0);
		uint32_t block_size, grid_size, shared_mem = 0;
		if (pa.n_alive < static_cast<uint32_t>(1024 * prop.multiProcessorCount))
		{
			block_size = static_cast<uint32_t>((pa.n_alive / prop.multiProcessorCount + 31) / 32 * 32);
			if (block_size == 0) block_size = 32;
		}
		else
		{
			block_size = 1024;
		}
		grid_size = static_cast<uint32_t>((pa.n_alive + block_size - 1) / block_size);

		MVSKernel_<<<grid_size, block_size, shared_mem, stream>>>(
				pa.r.data().get(),
				pa.v.data().get(),
				pa.deathflags.data().get(),
				device_particle_a.data().get(),
				pa.deathtime_index.data().get(),
				static_cast<uint32_t>(pa.n_alive),
				static_cast<uint32_t>(pl.log_len),
				static_cast<uint32_t>(pl.n_alive),
				device_h0_log(planet_data_id).data().get(),
				pl.r_log.data().get(),
				pl.m.data().get(),
				pl.id.data().get(),
				device_planet_rh.data().get(),
				base.outer_radius,
				dt,
				pl.m[0]
		);
#endif
	}
}
}
