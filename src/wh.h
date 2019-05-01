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

	struct WHIntegratorEncounterContinuationContext
	{
		std::unordered_map<uint32_t, uint8_t> id_to_enc_level;
	};

	class WHIntegrator
	{
	public:
		WHIntegratorEncounterContinuationContext ecc;
		Vf64 planet_inverse_helio_cubed, planet_inverse_jacobi_cubed;

		Vf64 particle_dist, particle_energy, particle_vdotr;
		Vf64 planet_dist, planet_energy, planet_vdotr;

		Vu8 particle_mask;
		Vu8 planet_mask;

		Vf64 planet_mu, planet_eta;
		Vf64 particle_mu;

		Vf64_3 planet_rj, planet_vj;
		Vf64_3 planet_a, particle_a;

		sr::util::LogQuartet<Vf64_3> planet_h0_log;

		Vf64 planet_rh;

		size_t encounter_n1, encounter_n2;
		double encounter_r1, encounter_r2;

		bool resolve_encounters;

		WHIntegrator();
		WHIntegrator(HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, const Configuration& config);

		void swap_logs();

		void integrate_planets_timeblock(HostPlanetPhaseSpace& pl, size_t nsteps, float64_t t, float64_t dt);
		// integrates for as many steps as are in the planet log
		void integrate_particles_timeblock(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, size_t begin, size_t length, float64_t t, float64_t dt);
		void gather_particles(const std::vector<size_t>& indices, size_t begin, size_t length);

		void load_h0(const HostPlanetPhaseSpace& pl);
		void integrate_encounter_particle_catchup(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, size_t particle_index, size_t particle_deathtime_index, double t, double dt);

		template<bool old>
		size_t integrate_encounter_particle_step(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, size_t particle_index, size_t timestep_index, uint8_t* encounter_level, double t, double dt);

		void step_planets(HostPlanetPhaseSpace& pl, float64_t t, double dt, size_t timestep_index);
		void step_particles(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, size_t begin, size_t length, float64_t t, double dt, size_t timestep_index);

		template<bool danby>
		static bool drift_single(float64_t t, float64_t mu, f64_3* r, f64_3* v);
		static bool drift_single_hp(float64_t t, float64_t mu, f64_3* r, f64_3* v);
		template<bool fixedit>
		static void drift(float64_t t, Vf64_3& r, Vf64_3& v, size_t start, size_t n, Vf64& dist, Vf64& energy, Vf64& vdotr, Vf64& mu, Vu8& mask);

		template<bool old>
		void nonhelio_acc_encounter_particle(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& p, size_t particle_index, float64_t time, size_t timestep_index, size_t central_planet_index);

		template<bool encounter, bool old>
		uint8_t helio_acc_particle(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& pa, size_t particle_index, float64_t time, size_t timestep_index);

		template<bool encounter, bool old>
		void helio_acc_particles(const HostPlanetPhaseSpace& pl, HostParticlePhaseSpace& p, size_t begin, size_t length, float64_t time, size_t timestep_index);

		template<bool slow>
		void helio_acc_planets(HostPlanetPhaseSpace& p, size_t index);

		static uint8_t detect_encounter(float64_t r_rel_sq, float64_t rh, double r1, double r2);
	};

	void calculate_planet_metrics(const HostPlanetPhaseSpace& p, double* energy, f64_3* l);
}
}
