#include "Density.h"

 
int main(int argc, char* argv[])
{
	Density den;
	den.readXYZR(argv[1], stod(argv[2]), stod(argv[3]), stod(argv[4]));
	den.createGrid();
	den.fillDensity();
	den.writeDensity();
	den.createLevelSetSurface();
	den.writeMesh();

	return 1;
}