#include "CircularAperture.h"

#include <cmath>
#include <cstdio>
#include <algorithm>
#include <assert.h>

using namespace std;
constexpr double pi = 3.14159265358979324;


//Evaluates h*dr*r for IntegrateCircle at all integration points, to save repeated work. Returns r step size.
static double CircularAperture_IntegrateCircle_hrdr(double rMax, double FocalLength, double Wavelength, vector<double> &hrdr_re, vector<double> &hrdr_im){
	//Wave number.
	const double k = 2.0*pi/Wavelength;

	//Integration step size in r direction. Provided light Wavelength is the worst case wavelength of what's being integrated, so a small fraction of that is perfect and always applicable.
	const double dr = Wavelength/50.0;

	//Loop over points.
	const size_t np = (size_t)(rMax/dr) + 2;
	hrdr_re.resize(np);
	hrdr_im.resize(np);

	for(size_t p = 0; p != np; p++){
		double r = dr*p;

		//Get R = |(x, y, z)|.
		double R = hypot(r, FocalLength);

		//Evaluate h(r).
		double kR = k*R;
		double zRR = FocalLength/(R*R*Wavelength);

		double s = sin(kR);
		double c = cos(kR);

		double hre = ( s + c/kR)*zRR;
		double him = (-c + s/kR)*zRR;

		hrdr_re[p] = hre*dr*r;
		hrdr_im[p] = him*dr*r;
	}

	return dr;
}

//This numerically integrates E0*h where E0 is the indicator function of a circle at a distance of AxisOffset from the 3D line x = y = 0. Uses results from IntegrateCircle_hrdr.
//This integral may be understood in different ways (eg, gives emission from circular function as a function of r when r = AxisOffset), but is intended as a building block integral.
static void CircularAperture_IntegrateCircle(double ApertureRadius, double FocalLength, double Wavelength, double AxisOffset, double *sum_re, double *sum_im, const double *hrdr_re, const double *hrdr_im, const double dr){
	assert(AxisOffset >= 0.0);

	//Limits of integration in r direction, corresponding indices.
	double rs = max(AxisOffset - ApertureRadius, 0.0) + dr/4.0;
	double re =     AxisOffset + ApertureRadius       - dr/4.0;

	const size_t ps = (size_t)ceil( rs/dr);
	const size_t pe = (size_t)floor(re/dr);

	//Loop over r. Oh, and to be clear, r = |(x, y)|, ie, 2D distance from optic axis.
	*sum_re = 0.0;
	*sum_im = 0.0;

	for(size_t p = ps; p <= pe; p++){
		/* Craft the integration size in tangential direction. Requires intersection of circle at origin with radius r, and circle a distance of AxisOffset at radius Radius.
		That problem can be translated into a triangle with known side lengths and unknown angles. Cosine rule reveals angle we're interested in. */
		double dt = 2.0*pi;

		if(AxisOffset > 0.0 && p != 0){
			double r = p*dr;
			double c = (r*r + AxisOffset*AxisOffset - ApertureRadius*ApertureRadius)/(2.0*r*AxisOffset);

			if(abs(c) <= 1.0)
				dt = 2.0*acos(c);
		}

		//Add, finally.
		*sum_re += hrdr_re[p]*dt;
		*sum_im += hrdr_im[p]*dt;
	}
}


CircularAperture::CircularAperture(double Resolution, double SensorWidth, double Radius, double FocalLength, double Wavelength) : 
		Radius(Radius), FocalLength(FocalLength), Wavelength(Wavelength), 
		rStep(Resolution*Wavelength), rMax(1.01*sqrt(2.0)*SensorWidth*SensorWidth + 4.0*Wavelength){ //Translate sensor size to a limit on distance from optic axis.

	//Size the memory.
	nI = (size_t)ceil(rMax/rStep) + 1;

	vI_re.resize(nI);
	vI_im.resize(nI);

	I_re = vI_re.data();
	I_im = vI_im.data();

	//Like above, but not const.
	auto Ire = vI_re.data();
	auto Iim = vI_im.data();

	//Store h, note need for extending rMax.
	vector<double> hrdr_re, hrdr_im;
	double dr = CircularAperture_IntegrateCircle_hrdr(rMax + Radius + 10.0*Wavelength, FocalLength, Wavelength, hrdr_re, hrdr_im);

	//Convolve.
	for(size_t i = 0; i != nI; i++){
		double r = i*rStep;

		double sum_re, sum_im;
		CircularAperture_IntegrateCircle(Radius, FocalLength, Wavelength, r, &sum_re, &sum_im, hrdr_re.data(), hrdr_im.data(), dr);

		Ire[i] = sum_re;
		Iim[i] = sum_im;
	}
}

bool CircularAperture::GetIntegral(double x, double y, double *sum_re, double *sum_im, bool Subtract){
	double r = hypot(x, y);

	//Linear interpolate with respect to r.
	r = r/rStep;
	size_t i = (size_t)r;
	r -= (double)i;

	if(i >= nI)
		return false;

	if(!Subtract){
		*sum_re += I_re[i + 1]*r + I_re[i]*(1.0 - r);
		*sum_im += I_im[i + 1]*r + I_im[i]*(1.0 - r);
	}else{
		*sum_re -= I_re[i + 1]*r + I_re[i]*(1.0 - r);
		*sum_im -= I_im[i + 1]*r + I_im[i]*(1.0 - r);
	}

	return true;
}

