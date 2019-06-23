#pragma once
#include "data.h"
#include "util.h"

#include <unordered_map>

namespace sr
{
namespace wh
{
	using namespace sr::data;

	bool kepeq(double dM, double ecosEo, double esinEo, double* dE, double* sindE, double* cosdE, uint32_t* iterations);
	bool kepeq_fixed(double dM, double ecosEo, double esinEo, double* dE, double* sindE, double* cosdE, uint32_t iterations);

	void kepmd(double dm, double es, double ec, double* x, double* s, double* c);
	bool kepu(double dt, double r0, double mu, double alpha, double u, double* fp, double* c1, double* c2, double* c3);

	class WHIntegrator
	{
	public:
		Vf64 planet_inverse_helio_cubed, planet_inverse_jacobi_cubed;

		Vf64 particle_dist, particle_energy, particle_vdotr;
		Vf64 planet_dist, planet_energy, planet_vdotr;

		Vu8 particle_mask;
		Vu8 planet_mask;

		Vf64 planet_mu, planet_eta;
		Vf64 particle_mu;

		Vf64_3 planet_rj, planet_vj;
		Vf64_3 planet_a, particle_a;

		sr::util::History<Vf64_3> planet_h0_log;

		Vf64 planet_rh;

		double outer_radius, cull_radius, encounter_sphere_factor;

		WHIntegrator();
		WHIntegrator(HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, const Configuration& config);

		void recalculate_rh(const HostPlanetPhaseSpace& pl);

		void swap_logs();

		void integrate_planets_timeblock(HostPlanetPhaseSpace& pl, size_t nsteps, float64_t t, float64_t dt);

		// integrates for as many steps as are in the planetary log: i.e. nsteps = pl.r_log().len
		void integrate_particles_timeblock(
				const HostPlanetPhaseSpace& pl,
				HostParticlePhaseSpace& pa,
				size_t begin,
				size_t length,
				float64_t t,
				float64_t dt
		);

		void gather_particles(const Vs& indices, size_t begin, size_t length);

		void load_h0(const HostPlanetPhaseSpace& pl);

		void step_planets(HostPlanetPhaseSpace& pl, float64_t t, double dt, size_t timestep_index);

		void step_particles(
				const HostPlanetPhaseSpace& pl,
				HostParticlePhaseSpace& pa,
				size_t begin,
				size_t length,
				float64_t t,
				double dt,
				size_t timestep_index
		);

		// this is just the public interface to the other private function
		void helio_acc_particles(
				const HostPlanetPhaseSpace& pl,
				HostParticlePhaseSpace& pa,
				size_t begin,
				size_t len,
				float64_t time,
				size_t timestep_index,
				bool old
		);

		/*
		   Calculates accelerations for planets, and stores h0 into
		   the planetary log at the given index

		   Side effects:
		   planet_inverse_helio_cubed
		   planet_inverse_jacobi_cubed
		   planet_a
		   planet_h0_log
	   	*/
		void helio_acc_planets(HostPlanetPhaseSpace& p, size_t index);

	private:
		template<bool danby>
		static bool drift_single(float64_t t, float64_t mu, f64_3* r, f64_3* v);

		static bool drift_single_hp(float64_t t, float64_t mu, f64_3* r, f64_3* v);

		template<bool fixedit>
		static void drift(
				float64_t t,
				Vf64_3& r,
				Vf64_3& v,
				size_t start,
				size_t n,
				Vf64& dist,
				Vf64& energy,
				Vf64& vdotr,
				Vf64& mu,
				Vu8& mask
		);

		template<bool old>
		void helio_acc_particle(
				const HostPlanetPhaseSpace& pl,
				HostParticlePhaseSpace& pa,
				size_t particle_index,
				float64_t time,
				size_t timestep_index
		);

		/*
		   Calculates the particle accelerations using planet_h0_log and planet_r_log
		   at the given timestep index. Also detects encounters.

		   Side effects:
		   particle_a
		   particles.deathflags
		   particles.deathtime
		 */
		template<bool old>
		void helio_acc_particles(
				const HostPlanetPhaseSpace& pl,
				HostParticlePhaseSpace& p,
				size_t begin,
				size_t length,
				float64_t time,
				size_t timestep_index
		);
	};

	void calculate_planet_metrics(const HostPlanetPhaseSpace& p, double* energy, f64_3* l);
}
}
