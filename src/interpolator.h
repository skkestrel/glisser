namespace sr
{
	class Interpolator
	{
		Interpolator(std::string file);

		private:
		std::ostream input;
		Vf64_3 aei0, aei1;
		Vf64_3 oof0, oof1;
		float64_t t0, t1;
	}
}
