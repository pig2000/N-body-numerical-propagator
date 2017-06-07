#include <iostream>
#include <boost/numeric/odeint.hpp>
#include <math.h>
#include <cspice/include/SpiceUsr.h>
#include <fstream>


using namespace std;
using namespace boost::numeric::odeint;




double Tx, Ty, Tz, ueq, GM, t0, dt, ax, ay, az;
vector<double> r_x, r_y, r_z, v_x, v_y, v_z, a_x, a_y, a_z, final_state, Time, G_M;
vector<int> body_list;
int list_size, frame;

typedef boost::array < double, 8 > state_type;
typedef runge_kutta_fehlberg78< state_type > error_stepper_type;
typedef controlled_runge_kutta< error_stepper_type > controlled_stepper_type;
controlled_stepper_type controlled_stepper;

void input_file(state_type &X0)
{
	double a[13];
	ifstream fin("input.txt");
	for (int i = 0; i<14; i++)
	{
		fin >> a[i];
	}
	fin.close();
	//X0 = [x (km), y (km), z (km), v_x (km), v_y (km), v_z(km), m (kg), s (km)]
	X0[0] = a[0];
	X0[1] = a[1];
	X0[2] = a[2];
	X0[3] = a[3];
	X0[4] = a[4];
	X0[5] = a[5];
	X0[6] = a[6];
	X0[7] = 0;
	t0=a[7]; //start time(s)
	Tx = a[8]; //thrust, x component(N)
	Ty = a[9]; //thrust, y component(N)
	Tz = a[10]; //thrust, z component(N)
	ueq = a[11]; //equivalent engine exhaust velocity,(m/s)
	dt = a[12]; //duration(s)
	frame=a[13];
}

void update_bodies()
{		
	int id,i;	
	ifstream body_file("body_list.txt");
	while (!body_file.eof())
	{
		body_file >> id;
		body_list.push_back(id);
	}
	body_file.close();
	list_size = body_list.size()-1;
	SpiceDouble GM;
	char k[10];
	for (int i = 0; i < list_size; i++)
	{
		sprintf(k,"%d",body_list[i]);
		SpiceInt dim = 1;
		bodvrd_c(k, "GM", 1, &dim, &GM);
		G_M.push_back(GM);
	}
}

void N_body(const state_type &X0, state_type &Xdot, const double t)
{		
	Xdot[0] = X0[3];
	Xdot[1] = X0[4];
	Xdot[2] = X0[5];
	Xdot[3] = ax;
	Xdot[4] = ay;
	Xdot[5] = az;
	double X,Y,Z;
	char k[10];
	int i;
	SpiceDouble state[6];	//state of the celestial body
	SpiceDouble et = t;	//current time
	SpiceDouble lt = 0;
	for (int j = 0; j < list_size; j++)
	{
		sprintf(k,"%d",body_list[j]);	//get the ID of the celestial body
		spkezr_c(k, et, "J2000", "NONE", "SSB", state, &lt);	//get the state of the celestial body
		X = X0[0] - state[0];	//position difference between the spacecraft and the celestial body
		Y = X0[1] - state[1];
		Z = X0[2] - state[2];
		Xdot[3] = Xdot[3] - X * G_M[j] / pow((X * X + Y * Y + Z * Z), 1.5);	//calculate the gravitational acceleration
		Xdot[4] = Xdot[4] - Y * G_M[j] / pow((X * X + Y * Y + Z * Z), 1.5);
		Xdot[5] = Xdot[5] - Z * G_M[j] / pow((X * X + Y * Y + Z * Z), 1.5);
	}
	Xdot[6] = -1.0 * pow(Tx*Tx + Ty*Ty + Tz*Tz, 0.5) / ueq;
	Xdot[7] = 1.0 / pow(X0[0] * X0[0] + X0[1] * X0[1] + X0[2] * X0[2], 0.5);	
}


void convert_to_inertial(const state_type &X0)
{
	double R[3] = { X0[0], X0[1], X0[2] };
	double normR = pow(R[0] * R[0] + R[1] * R[1] + R[2] * R[2], 0.5);
	double Rhat[3] = { R[0] / normR, R[1] / normR, R[2] / normR };
	double V[3] = { X0[3], X0[4], X0[5] };
	double normV = pow(V[0] * V[0] + V[1] * V[1] + V[2] * V[2], 0.5);
	double Vhat[3] = { V[0] / normV, V[1] / normV, V[2] / normV };
	double Nhat[3];
	if ((1 - abs(R[0] * V[0] + R[1] * V[1] + R[2] * V[2])) / (normR * normV) < 1E-08)
	{
		double Nhat[3] = { 0, 0, 1 };
	}
	else
	{
		double N[3] = { R[1] * V[2] - R[2] * V[1], R[2] * V[0] - R[0] * V[2], R[0] * V[1] - R[1] * V[0] };
		double normN = pow(N[0] * N[0] + N[1] * N[1] + N[2] * N[2], 0.5);
		double Nhat[3] = { N[0] / normN, N[1] / normN, N[2] / normN };
	}
	double Chat[3] = { Nhat[1] * Vhat[2] - Nhat[2] * Vhat[1], Nhat[2] * Vhat[0] - Nhat[0] * Vhat[2], Nhat[0] * Vhat[1] - Nhat[1] * Vhat[0]};
	double h[3] = { R[1] * V[2] - R[2] * V[1], R[2] * V[0] - R[0] * V[2], R[0] * V[1] - R[1] * V[0] };
	double normh = pow(h[0] * h[0] + h[1] * h[1] + h[2] * h[2], 0.5);
	double hhat[3] = { h[0] / normh, h[1] / normh, h[2] / normh };
	double theta_hat[3] = { hhat[1] * Rhat[2] - hhat[2] - Rhat[1], hhat[2] * Rhat[0] - hhat[0] * Rhat[2], hhat[0] * Rhat[1] - hhat[1] * Rhat[0] };
	switch (frame)
	{
	case 1:
		ax = Tx / (X0[6]*1000);
		ay = Ty / (X0[6]*1000);
		az = Tz / (X0[6]*1000);
		break;
	case 2:
		ax = (Vhat[0] * Tx + Chat[0] * Ty + Nhat[0] * Tz) / (X0[6]*1000);
		ay = (Vhat[1] * Tx + Chat[1] * Ty + Nhat[1] * Tz) / (X0[6]*1000);
		az = (Vhat[2] * Tx + Chat[2] * Ty + Nhat[2] * Tz) / (X0[6]*1000);
		break;
	case 3:
		ax = (Rhat[0] * Tx + theta_hat[0] * Ty + hhat[0] * Tz)/(X0[6]*1000);
		ay = (Rhat[1] * Tx + theta_hat[1] * Ty + hhat[1] * Tz)/(X0[6]*1000);
		az = (Rhat[2] * Tx + theta_hat[2] * Ty + hhat[2] * Tz)/(X0[6]*1000);
		break;
	}
}

void write_cout(const state_type &X0, const double t)
{
	r_x.push_back(X0[0]);
	r_y.push_back(X0[1]);
	r_z.push_back(X0[2]);
	v_x.push_back(X0[3]);
	v_y.push_back(X0[4]);
	v_z.push_back(X0[5]);
	Time.push_back(t);
	if(t==dt + t0)
	{
		final_state.push_back(X0[0]);
		final_state.push_back(X0[1]);
		final_state.push_back(X0[2]);
		final_state.push_back(X0[3]);
		final_state.push_back(X0[4]);
		final_state.push_back(X0[5]);
		final_state.push_back(X0[6]);
		final_state.push_back(t);
	}	
}

void output_file()
{	
	const int r_size = r_x.size();
	int i;
	ofstream R_x, R_y, R_z, V_x, V_y, V_z,T_file,finalstate;
	R_x.open("r_x.txt");
	for (i = 0; i < r_size; i++)
	{
		R_x << std::setprecision(20) << r_x[i] << endl;
	}
	R_x.close();
	R_y.open("r_y.txt");
	for (i = 0; i < r_size; i++)
	{
		R_y << std::setprecision(20) << r_y[i] << endl;
	}
	R_y.close();
	R_z.open("r_z.txt");
	for (i = 0; i < r_size; i++)
	{
		R_z << std::setprecision(20) << r_z[i] << endl;
	}
	R_z.close();
	V_x.open("v_x.txt");
	for (i = 0; i < r_size; i++)
	{
		V_x << std::setprecision(20) << v_x[i] << endl;
	}
	V_x.close();
	V_y.open("v_y.txt");
	for (i = 0; i < r_size; i++)
	{
		V_y << std::setprecision(20) << v_y[i] << endl;
	}
	V_y.close();
	V_z.open("v_z.txt");
	for (i = 0; i < r_size; i++)
	{
		V_z << std::setprecision(20) << v_z[i] << endl;
	}
	V_z.close();
	T_file.open("t.txt");
	for (i = 0; i < r_size; i++)
	{
		T_file << std::setprecision(20) << Time[i] << endl;
	}
	T_file.close();
	finalstate.open("finalstate.txt");
	for (i =0; i < 8; i++)
	{
		finalstate << std::setprecision(20) << final_state[i] << endl;
	}
	finalstate.close();
}

int main()
{
	state_type X0;
	input_file(X0);
	convert_to_inertial(X0);
	furnsh_c ("/usr/include/cspice/de423.bsp");
	furnsh_c ("/usr/include/cspice/de430.bsp");
	furnsh_c ("/usr/include/cspice/naif0010.tls.pc");
	furnsh_c ("/usr/include/cspice/pck00010.tpc");
	furnsh_c ("/usr/include/cspice/de-403-masses.tpc");
	update_bodies();
	integrate_adaptive(make_controlled< error_stepper_type >( 1.0e-14 , 1.0e-14 ), N_body, X0, t0, dt + t0, 1.0e-10, write_cout);
	output_file();
}
