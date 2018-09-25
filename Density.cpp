#include "Density.h"

int Density::readXYZR(string fn, double gs, double ex, double iso)
{
	ifstream in(fn.c_str());

	string str;
	while (getline(in, str)) {
		stringstream ss(str);

		double x, y, z, r;
		ss >> x >> y >> z >> r;
		//cout << "Atoms: " << x << " " << y << " " << z << " " << r << endl;
		vec3 v(x, y, z);
		
		atoms.push_back(v);
		radius.push_back(r);
	}

	gridsize = gs;
	extension = ex;
	isovalue = iso;

	return 1;
}

int Density::createGrid()
{
	// find bounding box;
	xmin = atoms[0].x();
	ymin = atoms[0].y();
	zmin = atoms[0].z();
	xmax = atoms[0].x();
	ymax = atoms[0].y();
	zmax = atoms[0].z();

	for (int i = 1; i < atoms.size(); ++i) {
		if (xmin > atoms[i].x())
			xmin = atoms[i].x();
		if (ymin > atoms[i].y())
			ymin = atoms[i].y();
		if (zmin > atoms[i].z())
			zmin = atoms[i].z();

		if (xmax < atoms[i].x())
			xmax = atoms[i].x();
		if (ymax < atoms[i].y())
			ymax = atoms[i].y();
		if (zmax < atoms[i].z())
			zmax = atoms[i].z();
	}

	xmin -= extension;
	ymin -= extension;
	zmin -= extension;
	xmax += extension;
	ymax += extension;
	zmax += extension;
	cout << xmin << " " << ymin << " " << zmin << endl;
	cout << xmax << " " << ymax << " " << zmax << endl;
	cout << gridsize << endl;

	cout << xmax - xmin << endl;
	xdim = ceil((xmax - xmin) / gridsize);
	ydim = ceil((ymax - ymin) / gridsize);
	zdim = ceil((zmax - zmin) / gridsize);
	cout << xdim << " " << ydim << " " << zdim << endl;

	xmax = xmin + xdim * gridsize;
	ymax = ymin + ydim * gridsize;
	zmax = zmin + zdim * gridsize;

	voxs.resize(xdim*ydim*zdim, 0);

	return 1;
}

int Density::fillDensity()
{

	for (int l = 0; l < atoms.size(); ++l) {
		for (int k = 0; k < zdim; ++k) {
			for (int j = 0; j < ydim; ++j) {
				for (int i = 0; i < xdim; ++i) {
					vec3 pos(xmin + (i + 0.5)*gridsize, ymin + (j + 0.5)*gridsize, zmin + (k + 0.5)*gridsize);
					voxs[i + j * xdim + k * xdim*ydim] += getDensity(pos, atoms[l], radius[l]);
					//cout << voxs[i + j * xdim + k * xdim*ydim] << endl;
				}
			}
		}
	}

	return 1;
}

int Density::writeDensity()
{
	ofstream out("density.vox");

	out << xmin << " " << ymin << " " << zmin << endl;
	out << xmax << " " << ymax << " " << zmax << endl;
	out << xdim << " " << ydim << " " << zdim << endl;

	for (int i = 0; i < voxs.size(); ++i) {
		if (voxs[i] > 1)
			voxs[i] = 1;
		if (voxs[i] < 0)
			voxs[i] = 0;

		out << voxs[i] << endl;
	}

	return 1;
}

int Density::createLevelSetSurface()
{
	openvdb::initialize();
	//cout << "11111111111111"<< endl;
	grid = openvdb::DoubleGrid::create();
	//cout << "11111111111111" << endl;
	openvdb::DoubleGrid::Accessor accessor = grid->getAccessor();
	//cout << "11111111111111" << endl;
	for (int k = 0; k < zdim; ++k) {
		for (int j = 0; j < ydim; ++j) {
			for (int i = 0; i < xdim; ++i) {
				openvdb::Coord xyz(i, j, k);
				//cout << i<<" " << j << " " << k<< endl;
				//cout << xdim<<" " << ydim << " " << zdim << endl;
				//cout << voxs[i + j*xdim + k*xdim*ydim] << endl;
				accessor.setValue(xyz, voxs[i + j*xdim + k*xdim*ydim]);
			}
		}
	}

	openvdb::tools::volumeToMesh<openvdb::DoubleGrid>(*grid, vertices, quads, isovalue);
}

int Density::writeMesh()
{
	ofstream out("mesh.obj");

	for (int i = 0; i < vertices.size(); ++i) {
		out << "v " << vertices[i].x() << " " << vertices[i].y() << " " << vertices[i].z() << endl;
	}

	for (int i = 0; i < quads.size(); ++i) {
		out << "f " << quads[i].x() + 1 << " " << quads[i].y() + 1 << " " << quads[i].z() + 1 << " " << quads[i].w() + 1 << endl;
	}

	return 1;
}

double Density::getDensity(vec3 pos, vec3 anchor, double radius)
{
	vec3 d = pos - anchor;

	double val = exp(-pow((sqrt(d.dot(d))) / radius, 2));

	return val;
}