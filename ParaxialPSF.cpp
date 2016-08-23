#include "CircularAperture.h"

using namespace std;


//From Wikipedia.
double OptimalZonePlaneTransition(double FocalLength, double Wavelength, size_t n){
	return sqrt(n*Wavelength*FocalLength + n*n*Wavelength*Wavelength/4.0);
}

int main(int ArgumentCount, char **ArgumentVector){
	//Description of imaging system. These are inputs. All length measurements are millimeters.
	static const double FocalLength    = 14.0;
	static const double PinHoleRadius  = 0.1;
	static const double SensorRadius   = 1.0;
	static const double Wavelength     = 0.5*1e-3;	//Greenish light.


	//Simulate the point spread function due to four pinholes next to each other, arranged in a cross/plus sign of size Separation.
	//Note that since the individual holes are the same size, we need only one CircularAperture. Just a slice through the x axis.
	CircularAperture ca(20.0, SensorRadius*2.0, PinHoleRadius, FocalLength, Wavelength);

	double Separation = 0.3;

	for(double x = -SensorRadius; x < SensorRadius; x += 0.01){
		double sum_re = 0.0, sum_im = 0.0;

		//By linearity of integration, the desired integral is just sum of effects from the four holes.
		ca.GetIntegral(x - Separation - PinHoleRadius, 0.0, &sum_re, &sum_im);
		ca.GetIntegral(x + Separation + PinHoleRadius, 0.0, &sum_re, &sum_im);
		ca.GetIntegral(x, -Separation - PinHoleRadius, &sum_re, &sum_im);
		ca.GetIntegral(x,  Separation + PinHoleRadius, &sum_re, &sum_im);

		printf("%f\n", hypot(sum_re, sum_im));
	}


/*
	//Zone plate.
	const size_t nZones = 9;

	const size_t N = 2*nZones - 1;
	size_t i;

	vector<CircularAperture *> CA;

	for(i = 0; i != N; i++)
		CA.push_back(new CircularAperture(5.0, SensorRadius*2.0, OptimalZonePlaneTransition(FocalLength, Wavelength, i + 1), FocalLength, Wavelength));

	for(double r = 0.0; r < SensorRadius; r += 0.001){
		//Compute sums.
		double sum_re = 0.0, sum_im = 0.0;

		for(i = 0; i != N; i++){
			bool negative = (i & 1) != 0;
			CA[i]->GetIntegral(r, 0.0, &sum_re, &sum_im, negative);

			//Print magnitude for every positive ring, so that we're outputting multiple graphs corresponding to different numbers of rings.
			double m = hypot(sum_re, sum_im);

			printf("%.3e,", m);
		}

		printf("\n");
	}

	for(i = 0; i != N; i++)
		delete CA[i];
*/


	return EXIT_SUCCESS;
}




