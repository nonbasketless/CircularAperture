#include <vector>


/* This class is intended for simulating the passage of light from a point at (0, 0, -ExtremelyFar) onto a sensor plane at z = FocalLength
after passage through a circular aperture (infinitely thin material) at z = 0.

The math is based on info from https://tomroelandts.com/articles/the-psf-of-a-pinhole-camera. Here's some Maxima code I used:

	k*z/(2*%pi*%i)*exp(%i*k*r)/r^2*(1 - 1/(%i*k*r));
	rectform(%);

	[realpart(%), imagpart(%)];

Normally, Fourier transforms are used; however here integrals are written in polar form (reducing complexity) and evaluated with midpoint rule. This gives much
better accuracy than Fourier methods and no edge effects, but only works for circular apertures. So you can simulate pinholes, photon sieves, zone plates, etc. */
class CircularAperture{
public:
	//SensorWidth is used to limit number of calculations done. Resolution controls fineness of output (not integration size), in multiples of Wavelength, smaller is finer/slower.
	CircularAperture(double Resolution, double SensorWidth, double Radius, double FocalLength, double Wavelength);

	//Retrieve integral (convolution point) at (x, y), adds to (or subtracts froms) sum. Requires 2*max(|x|, |y|) <= SensorWidth.
	bool GetIntegral(double x, double y, double *sum_re, double *sum_im, bool Subtract = false);

	const double Radius;
	const double FocalLength;
	const double Wavelength;

	//Computed convolutions.
	const double *I_re;
	const double *I_im;
	const double rStep, rMax;

private:
	size_t nI;
	std::vector<double> vI_re, vI_im;
};

