#include <iostream>
#include <algorithm>
#include <cstring>
#include <cmath>
#include <utility>

#include "GridSolver.h"

using std::vector;
using std::swap;
using std::deque;
using std::shared_ptr;
using std::make_shared;

float nv = 1.0f, kappa = 0.1f; // diffuse rate. nv is for velocity field and kapp is for heat transfer.
float dt = 0.01f; // simulation step size
float simTime;
float dx = 1; // size of a cell in simulation
float velScale = 10.0f; // scaling for display
float buoyancy = 0.5f; // buoyancy force coefficient
float refinementThreshold = 20.0f;
int iterations = 20; // number of iterations for iterative solvers (Gauss–Seidel, Conjugate Gradients, etc.)

/*
 * Quantities in the grid.
 * Each has 2 entries (e.g., vx[0] and vy[1] for storing the current/next state)
 * But you had better use pointers like cur_vx, next_vx, etc. defined below.
 * Note that rows are horizontal and columns are vertical.
 * e.g., cur_vx[3][4] is the x component of the velocity at position (4, 3).
 * Also note that the first rows/columns and last rows/columns are for the boundary.
 */
// velocity
double vx[2][rows + 2][cols + 2];
double vy[2][rows + 2][cols + 2];

// temperature
double tp[2][rows + 2][cols + 2];

// pointers for the corresponding quantities in the current and the next states
double (*cur_vx)[cols + 2] = vx[0], (*next_vx)[cols + 2] = vx[1];
double (*cur_vy)[cols + 2] = vy[0], (*next_vy)[cols + 2] = vy[1];
double (*cur_tp)[cols + 2] = tp[0], (*next_tp)[cols + 2] = tp[1];

// flow source
double vxSrc[rows + 2][cols + 2];
double vySrc[rows + 2][cols + 2];
// heat source
double tpSrc[rows + 2][cols + 2];

//filament
deque<shared_ptr<vector<Vertex>>> filaments;
deque<double> ages;
float maxAge = 10;

/*
 * Flow type:
 *	horizontal: horizontal flow
 *	vertical: vertical flow
 *	other: quantities like temperature, density, etc.
 */
enum FlowType { horizontal, vertical, other };
void setBnd(double a[rows + 2][cols + 2], FlowType flowType);

Vertex gridVertices[rows + 1][cols + 1];
GLuint gridIndices[6 * rows * cols];
Vertex velVertices[rows][cols][2];

void initGrid() {
	memset(vx, 0, sizeof(vx));
	memset(vy, 0, sizeof(vx));
	memset(tp, 0, sizeof(vx));
	memset(vxSrc, 0, sizeof(vxSrc));
	memset(vySrc, 0, sizeof(vySrc));
	memset(tpSrc, 0, sizeof(tpSrc));
	filaments.clear();
	ages.clear();
	cur_vx = vx[0];
	cur_vy = vy[0];
	cur_tp = tp[0];
	next_vx = vx[1];
	next_vy = vy[1];
	next_tp = tp[1];

	GLuint indices[rows + 1][cols + 1];
	GLuint idx = 0;
	for (size_t i = 0; i <= rows; ++i) {
		for (size_t j = 0; j <= cols; ++j) {
			gridVertices[i][j] = Vertex((float)j * cellSize, (float)i * cellSize, 0.f, 0.f, 0.f, 1.f);
			indices[i][j] = idx++;
		}
	}

	size_t k = 0;
	for (size_t i = 0; i < rows; ++i) {
		for (size_t j = 0; j < cols; ++j) {
			gridIndices[k++] = indices[i][j];
			gridIndices[k++] = indices[i][j + 1];
			gridIndices[k++] = indices[i + 1][j + 1];
			gridIndices[k++] = indices[i][j];
			gridIndices[k++] = indices[i + 1][j + 1];
			gridIndices[k++] = indices[i + 1][j];

			velVertices[i][j][0] = Vertex((float)(j + 0.5) * cellSize, (float)(i + 0.5) * cellSize, 0.f, 0.f, 1.f, 1.f);
			velVertices[i][j][0] = Vertex((float)(j + 0.5) * cellSize, (float)(i + 0.5) * cellSize, 0.f, 0.f, 1.f, 1.f);
		}
	}

	updateGrid();
	simTime = 0;
}

/*
 * Unit square bilinear-interpolation (https://en.wikipedia.org/wiki/Bilinear_interpolation#On_the_unit_square).
 * You need this for interpolating the value in the cell (e.g., in backtracing advect).
 * Parameters:
 *	x, y:
 *	    position of the point whose value you want to interpolate
 *	v00:
 *	    value at (0, 0)
 *	v01:
 *	    value at (0, 1)
 *	v10:
 *	    value at (1, 0)
 *	v11:
 *	    value at (1, 1)
 *  return:
 *	interpolated value
 */
double bilinearInterpolate(double x, double y, double v00, double v01, double v10, double v11) {
	// DONE: Do a unit square bilinear interpolation for a point at (x, y).
	return v00 * (1-x) * (1-y) + v10 * (x) * (1 - y) + v01 * (1 - x) * (y) + v11 * (x) * (y);
}

/*
 * Set up the boundary.
 * Parameters:
 *	a:
 *	    2D array whose boundary need to be set
 *	flowType:
 *	    type of flow: horizontal flow (the interpolated quantity will goes to zero at the vertical boundaries), vertical flow (the interpolated quantity will goes to zero at the horizontal boundaries), other (e.g., temperature, density, which are not really flow... only continuity need to be guaranteed.)
 */
void setBnd(double a[rows + 2][cols + 2], FlowType flowType) {
	// DONE: Set up the boundary according to the flow type.
	for (int i = 1; i < rows + 1; i++)
	{
		a[i][0] = flowType == horizontal ? -a[i][1] : a[i][1];
		a[i][cols + 1] = flowType == horizontal ? -a[i][cols] : a[i][cols];
	}
	for (int j = 1; j < cols + 1; j++)
	{
		a[0][j] = flowType == vertical ? -a[1][j] : a[1][j];
		a[rows + 1][j] = flowType == vertical ? -a[rows][j] : a[rows][j];
	}
	a[0][0] = 0.5 * (a[1][0] + a[0][1]);
	a[0][cols + 1] = 0.5 * (a[1][cols + 1] + a[0][cols]);
	a[rows + 1][0] = 0.5 * (a[rows][0] + a[rows + 1][1]);
	a[rows + 1][cols + 1] = 0.5 * (a[rows][cols + 1] + a[rows + 1][cols]);
}

/*
 * Add source to the field.
 * Parameters:
 *	src:
 *	    2D array containing the source
 *	a:
 *	    target 2D array storing the quantity to modify
 */
void addSource(double src[rows + 2][cols + 2], double a[rows + 2][cols + 2]) {
	// DONE: Add the source from array *src* to the target array *a*: e.g., a = a + dt * src
	for (int x = 0; x < rows + 2; x++)
	{
		for (int y = 0; y < cols + 2; y++)
		{
			a[x][y] = a[x][y] + dt * src[x][y];
		}
	}
}

/*
 * Compute the diffusion part.
 * Parameters:
 *	a0:
 *	    2D array storing the quantities in current state
 *	a1:
 *	    2D array storing the quantities in next state
 *	nv:
 *	    diffusion rate (nv/kappa in the equations)
 *	flowType:
 *	    flow type
 */
void diffuse(double a0[rows + 2][cols + 2], double a1[rows + 2][cols + 2], double nv, FlowType flowType) {
	// DONE: diffusion
	// Compute the diffusion part and update array *a1*.
	// Use dt and dx for step size and cell size.
	// Do a implicit solve for stability.
	// Use a Gauss Seidel solve (or better, Conjugate Gradients).
	// Use *iterations* as the number of iterations.
	// Call setBnd to fix the boundary.
	double a = dt * nv / dx;
	for (int k = 0; k < iterations; k++)
	{
		for (int i = 1; i <= rows; i++)
		{
			for (int j = 1; j <= cols; j++)
			{
				a1[i][j] = (a0[i][j] + a * ((a1[i - 1][j] + a1[i + 1][j] + a1[i][j - 1] + a1[i][j + 1]) / 4.0)) /
					(1 + a);
			}
		}
		setBnd(a1, flowType);
	}
}

/*
 * Compute the advection part.
 * Parameters:
 *	a0:
 *	    2D array storing the quantities in current state
 *	a1:
 *	    2D array storing the quantities in next state
 *	vx, vy:
 *	    2D arrays storing the velocity field
 *	flowType:
 *	    flow type
 */
void advect(double a0[rows + 2][cols + 2], double a1[rows + 2][cols + 2], double vx[rows + 2][cols + 2], double vy[rows + 2][cols + 2], FlowType flowType) {
	// DONE: advection
	// Compute the advection part and update array *a1*.
	// Use dt and dx for step size and cell size.
	// Do a linear (or better, higher order or adaptive) backtrace for each center of the cells.
	// Compute the quantity in the previous state using the bilinear interpolation.
	// Call setBnd to fix the boundary.r
	for (int i = 1; i <= rows; i++) {
		for (int j = 1; j <= cols; j++) {
			double x = i - dt * (1 / (dx * dx)) * vy[i][j];
			double y = j - dt * (1 / (dx * dx)) * vx[i][j];
			if (x < 0.5) x = 0.5; if (x > rows + 0.5) x = rows + 0.5;
			if (y < 0.5) y = 0.5; if (y > cols + 0.5) y = cols + 0.5;
			int i0 = (int)x;
			int j0 = (int)y;
			int i1 = i0 + 1;
			int j1 = j0 + 1;
			a1[i][j] = bilinearInterpolate(x-i0, y-j0, a0[i0][j0], a0[i0][j1], a0[i1][j0], a0[i1][j1]);
			//std::cout << "Before: " << a0[i][j] << " After: " << a1[i][j] << std::endl;
		}
	}
	setBnd(a1, flowType);

}

double s_p[rows + 2][cols + 2], s_div[rows + 2][cols + 2];

/*
 * Projection for the mass conservation.
 * Parameter:
 *	vx, vy:
 *	    the velocity field to be fixed
 */
void project(double vx[rows + 2][cols + 2], double vy[rows + 2][cols + 2]) {
	// DONE: projection
	// Do a Poisson Solve to get a divergence free velocity field.
	// Use a Gauss Seidel solve (or better, Conjugate Gradients).
	// Use *iterations* as the number of iterations.
	// Call setBnd to fix the boundary.
	double div[rows + 2][cols + 2];
	double p[rows + 2][cols + 2];
	for (int i = 1; i <= rows; i++) {
		for (int j = 1; j <= cols; j++) {
			div[i][j] = 0.5 * dx * (vx[i + 1][j] - vx[i - 1][j] +
				vy[i][j + 1] - vy[i][j - 1]);
			p[i][j] = 0;
		}
	}
	setBnd(div, other); setBnd(p, other);
	for (int k = 0; k < iterations; k++) {
		for (int i = 1; i <= rows; i++) {
			for (int j = 1; j <= cols; j++) {
				p[i][j] = (-div[i][j] + p[i - 1][j] + p[i + 1][j] +
					p[i][j - 1] + p[i][j + 1]) / 4;
			}
		}
		setBnd(p, other);
	}
	for (int i = 1; i <= rows; i++) {
		for (int j = 1; j <= cols; j++) {
			vx[i][j] -= 0.5 * (p[i + 1][j] - p[i - 1][j])/dx;
			vy[i][j] -= 0.5 * (p[i][j + 1] - p[i][j - 1])/dx;
		}
	}
	setBnd(vy, vertical); setBnd(vx, horizontal);
}

/*
 * Get the reference (average) temperature of the grid.
 * Parameters:
 *	tp:
 *	    2D array storing the temperatures
 */
double getReferenceTemperature(double tp[rows + 2][cols + 2]) {
	// TODO: Sum up array *tp* and compute the average.
	double average = 0;
	for (int i = 0; i < rows + 2; i++)
	{
		for (int j = 0; j < cols + 2; j++)
		{
			average += tp[i][j];
		}
	}
	return average/((rows + 2)*(cols + 2));
}

/*
 * Apply the buoyancy force due to the temperature difference.
 * Parameters:
 *	a:
 *	    2D array storing the velocity
 *	tp:
 *	    2D array storing the temperature
 *	beta:
 *	    buoyancy force coefficient
 *	flowType:
 *	    flow type
 */
void applyTemperatureForce(double a[rows + 2][cols + 2], double tp[rows + 2][cols + 2], double beta, FlowType flowType) {
	// TODO: buoyancy forces
	// Apply the buoyancy force and update array *a*.
	// For more details, see Foster and Metaxas [1997] Equation 2.
	double t0 = getReferenceTemperature(tp);
	for (int i = 0; i < cols + 2; i++)
	{
		for (int j = 0; j < rows + 1; j++)
		{
			a[i][j] += -beta * (t0 - ((tp[i][j] + tp[i][j + 1]) / 2));
		}
	}
	
	setBnd(a, flowType);
}

/*
 * One stimulation step.
 */
void step() {
	// TODO: step the simulation
	// Compute and update the velocity and temperature (in cur_vx,cur_vy, cur_tp) in the next state based on the current state.
	// You need to apply source, diffuse, advect forces for temperature and velocity fields.
	// For velocity field, you need to do the projection to get a divergence free field.
	// You also need to apply the buoyancy force to the velocity field.
	// Don't forget to swap pointers (e.g., cur_vx and next_vx, etc.)!
	// Have a look at GDC03.
	// Change the paramters (dt, dx, etc.) in the setting panel to check whether your solver is stable!

	addSource(vxSrc, next_vx);
	addSource(vySrc, next_vy);
	addSource(vxSrc, cur_vx);
	addSource(vySrc, cur_vy);
	diffuse(next_vx, cur_vx, nv, horizontal);
	diffuse(next_vy, cur_vy, nv, vertical);
	project(cur_vx, cur_vy);
	diffuse(cur_vx, next_vx, nv, horizontal);
	diffuse(cur_vy, next_vy, nv, vertical);
	project(next_vx, next_vy);

	addSource(tpSrc, next_tp);
	diffuse(next_tp, cur_tp, kappa, other);
	advect(cur_tp, next_tp, next_vx, next_vy, other);

	applyTemperatureForce(next_vy, cur_tp, dt*buoyancy, vertical);
	// Please DO NOT change the following
	simTime += dt;
	//updateGrid();
	updateFilament();
	memset(vxSrc, 0, sizeof(vxSrc));
	memset(vySrc, 0, sizeof(vySrc));
}

void updateGrid() {
	for (size_t i = 0; i <= rows; ++i)
		for (size_t j = 0; j <= cols; ++j) {
			float v00 = cur_tp[i][j];
			float v01 = cur_tp[i + 1][j];
			float v10 = cur_tp[i][j + 1];
			float v11 = cur_tp[i + 1][j + 1];
			// Weird interpolation because the quad is rendered as 2 triangles?
			// Maybe this can help: https://jcgt.org/published/0011/03/04/paper.pdf
			double t = bilinearInterpolate(0.5, 0.5, v00, v01, v10, v11);
			gridVertices[i][j].r = 0;
			gridVertices[i][j].g = 0;
			gridVertices[i][j].b = 0;
			if (t > 0) {
				gridVertices[i][j].r = (GLfloat)std::min(1.0, t);
			} else if (t < 0) {
				gridVertices[i][j].b = (GLfloat)std::min(1.0, -t);
			}
		}

	for (size_t i = 1; i <= rows; ++i)
		for (size_t j = 1; j <= cols; ++j) {
			velVertices[i - 1][j - 1][1].x = velVertices[i - 1][j - 1][0].x + cur_vx[i][j] * velScale;
			velVertices[i - 1][j - 1][1].y = velVertices[i - 1][j - 1][0].y + cur_vy[i][j] * velScale;
		}
}

double dist(const Vertex &a, const Vertex &b) {
    double x = a.x - b.x;
    double y = a.y - b.y;
    return std::sqrt(x*x + y*y);
}

void addFilament(double x, double y) {
	filaments.push_back(make_shared<vector<Vertex>>());
	for (size_t i = 0; i <= rows; ++i) {
		filaments.back()->push_back(Vertex(x, (float)i * cellSize, 1.0f, 1.0f, 1.0f, 1.0f));
	}
	ages.push_back(0);

	filaments.push_back(make_shared<vector<Vertex>>());
	for (size_t i = 0; i <= cols; ++i) {
		filaments.back()->push_back(Vertex((float)i * cellSize, y, 1.0f, 1.0f, 1.0f, 1.0f));
	}
	ages.push_back(0);
}

void updateFilament() {
	while (!ages.empty() && ages.front() > maxAge) {
		ages.pop_front();
		filaments.pop_front();
	}
	for (size_t i = 0; i < filaments.size(); ++i) {
		ages[i] += dt;
		shared_ptr<vector<Vertex>> tmp = make_shared<vector<Vertex>>();
		tmp->push_back(filaments[i]->front());
		for (size_t j = 1; j < filaments[i]->size(); ++j) {
			const Vertex &a = tmp->back();
			const Vertex &b = filaments[i]->at(j);
			if (dist(a, b) > refinementThreshold) {
				tmp->push_back(Vertex((a.x + b.x)/2, (a.y + b.y)/2, b.r, b.g, b.b, b.a));
			}
			tmp->push_back(b);
		}
		filaments[i] = tmp;
		for (Vertex &v: *filaments[i]) {
			double x = v.x/cellSize + 0.5;
			double y = v.y/cellSize + 0.5;
			size_t j0 = (size_t)x;
			size_t i0 = (size_t)y;
			size_t i1 = i0 + 1, j1 = j0 + 1;
			x -= j0;
			y -= i0;
			double v00, v01, v10, v11;
			double vx, vy;
			v00 = cur_vx[i0][j0];
			v01 = cur_vx[i1][j0];
			v10 = cur_vx[i0][j1];
			v11 = cur_vx[i1][j1];
			vx = bilinearInterpolate(x, y, v00, v01, v10, v11);
			v00 = cur_vy[i0][j0];
			v01 = cur_vy[i1][j0];
			v10 = cur_vy[i0][j1];
			v11 = cur_vy[i1][j1];
			vy = bilinearInterpolate(x, y, v00, v01, v10, v11);
			v.x += (float)vx * dt * cellSize;
			v.y += (float)vy * dt * cellSize;
			v.x = std::max(0.f, v.x);
			v.x = std::min((float)frameWidth, v.x);
			v.y = std::max(0.f, v.y);
			v.y = std::min((float)frameHeight, v.y);
			if (ages[i] > maxAge / 2)
			v.a = 2 - (float)ages[i] / (maxAge / 2);
		}
	}
}
