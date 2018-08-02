#include <functional>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <set>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <stdexcept>
#include "bending.h"
#include<../iFluid/ifluid_state.h>

static void DebugShow(const double&);
bool findAndLocate(std::ifstream&, const char*);
bool getString(std::ifstream&, const char*);
static double divEx(double, double); 

BendingForce::BendingForce(INTERFACE* _intfc, double s, double d) :
        intfc(_intfc) {
    bends = s;
    bendd = d;
    index = 2;
    method[0] = &BendingForce::calculateBendingForce3d2003;
    method[1] = &BendingForce::calculateBendingForce3d2006;
    method[2] = &BendingForce::calculateBendingForce3dparti;
}

void BendingForce::getParaFromFile(const char* inname)
{
    std::ifstream fin(inname);

        if (!fin.is_open())
        {
            std::cerr << "Can't open file!\n";
            clean_up(ERROR);
        }
        if (!findAndLocate(fin, "Enter fabric bend stiffness constant:"))
            bends = 0.0;
        else
            fin >> bends;
        std::cout << bends << std::endl;
        if (!findAndLocate(fin, "Enter fabric bend damping constant:"))
            bendd = 0.0;
        else
            fin >> bendd;
        std::cout << bendd << std::endl;
        if (!findAndLocate(fin, "Enter bending method index:"))
            index = 2;
        else
            fin >> index;
        std::cout << index << std::endl;
        fin.close();
}
double* BendingForce::getExternalForce(SpringVertex* sv) 
{
    	return ((POINT*)(sv->org_vtx))->force;
}

void BendingForce::computeExternalForce()
{
        SURFACE **surf;
        TRI *tri;

        intfc_surface_loop(intfc, surf)
        {
            if (wave_type(*surf) != ELASTIC_BOUNDARY) continue;
            if (is_bdry(*surf)) continue;
            clear_surf_point_force(*surf);
            surf_tri_loop(*surf, tri)
            {
                if (wave_type(*surf) != ELASTIC_BOUNDARY) continue;
                if (is_bdry(*surf)) continue;
                for (int i = 0; i < 3; ++i)
                {
                    POINT *p1 = Point_of_tri(tri)[i];

                    sorted(p1) = NO;
                    TRI* n_tri = Tri_on_side(tri, (i+1)%3);
                    if (is_side_bdry(tri, (i+1)%3)) continue;
                    (this->*method[methodIndex()])(p1, tri, n_tri);
                }
            }
        }
}       /* setBendingForce3d */

static void DebugShow(const double & sva)
{
	std::cout << std::setw(20) << sva << " ";
}

static double divEx(double numerator, double denominator) {
    if (fabs(denominator) < 1.0e-10) 
	throw std::overflow_error("Divide by zero exception!\n"); 
    return numerator / denominator; 
}

double BendingForce::calOriLeng(int index1, int index2, TRI* tri, TRI* n_tri) {
    double c = tri->side_length0[(index1+1)%3];
    double b1 = tri->side_length0[index1];
    double a1 = tri->side_length0[3-index1-(index1+1)%3];
    double a2, b2;
    POINT* p1 = Point_of_tri(tri)[(index1+1)%3];
    int pInN;

    for (int i = 0; i < 3; i++)
         if (Point_of_tri(n_tri)[i] == p1) {
             pInN = i;
             break;
         }
    if (pInN < index2 || pInN > index2 && pInN == 2) {
        b2 = n_tri->side_length0[pInN];
        a2 = n_tri->side_length0[index2];
    }
    else {
        b2 = n_tri->side_length0[index2];
        a2 = n_tri->side_length0[3-index2-pInN];
    }
    
    double cangle1, cangle2; 

    try {
        cangle1 = divEx(sqr(a1) + sqr(c) - sqr(b1), (2 * a1 * c));
    } catch (std::overflow_error e) {
	std::cout << e.what() << " -> "; 
	return a2; 
    }
    try {
        cangle2 = divEx(sqr(a2) + sqr(c) - sqr(b2), (2 * a2 * c));
    } catch (std::overflow_error e) {
	std::cout << e.what() << " -> "; 
	return a1; 
    }

    double sangle1 = sqrt(1 - std::min(sqr(cangle1), 1.0));
    double sangle2 = sqrt(1 - std::min(sqr(cangle2), 1.0));
    double cangle = cangle1 * cangle2 - sangle1 * sangle2;

    return sqrt(sqr(a1) + sqr(a2) - 2 * cangle * a1 * a2);
}

void BendingForce::calculateBendingForce3dparti(POINT* p1,
                TRI* tri, TRI* n_tri) {
        int index1, index2;
        std::set<POINT*> pointset;

        for (int i = 0; i < 3; i++) {
             if (Point_of_tri(tri)[i] == p1)
                 index1 = i;
             pointset.insert(Point_of_tri(tri)[i]);
        }

        POINT* p2;

        for (int i = 0; i < 3; i++)
             if (pointset.find(Point_of_tri(n_tri)[i]) == pointset.end()) {
                 p2 = Point_of_tri(n_tri)[i];
                 index2 = i;
                 break;
             }

        double length0 = calOriLeng(index1, index2, tri, n_tri);
        double length = separation(p1, p2, 3);
        STATE* state = static_cast<STATE*>(left_state(p1));
        double *vel1 = state->vel;
        state = static_cast<STATE*>(left_state(p2));
        double *vel2 = state->vel;
        double velr[3] = {0.0};
        double dir[3] = {0.0};

        for (int i = 0; i < 3; i++) {
             velr[i] = vel2[i] - vel1[i];
             dir[i] = Coords(p2)[i] - Coords(p1)[i];
        }
        if (length > 1.0e-10)
            for (int i = 0; i < 3; i++)
                 dir[i] /= length;
        for (int i = 0; i < 3; i++) {
             p1->force[i] += bends * (length - length0) * dir[i];
             p1->force[i] += bendd * velr[i];
        }
}

void BendingForce::calculateBendingForce3d2003(
        POINT* p1,
        TRI* t1,
        TRI* t2)
{
        if (Mag3d(Tri_normal(t1)) < MACH_EPS ||
            Mag3d(Tri_normal(t2)) < MACH_EPS) return;

        int index = Vertex_of_point(t1, p1);
        POINT *p3 = Point_of_tri(t1)[(index+1)%3];
        POINT *p4 = Point_of_tri(t1)[(index+2)%3];
        index = 3 - Vertex_of_point(t2, p3) - Vertex_of_point(t2, p4);
        POINT *p2 = Point_of_tri(t2)[index];
        double x13[3], x14[3], x23[3], x24[3], E[3];

	std::transform(Coords(p1), Coords(p1) + 3, Coords(p3), 
		x13, std::minus<double>());
	std::transform(Coords(p1), Coords(p1) + 3, Coords(p4), 
		x14, std::minus<double>());
	std::transform(Coords(p2), Coords(p2) + 3, Coords(p3), 
		x23, std::minus<double>());
	std::transform(Coords(p2), Coords(p2) + 3, Coords(p4), 
		x24, std::minus<double>());
	std::transform(Coords(p4), Coords(p4) + 3, Coords(p3), 
		E, std::minus<double>());

        double N1[3], N2[3], E_mag, N1_mag, N2_mag;
        double n1[3], n2[3];

        Cross3d(x13, x14, N1);
        Cross3d(x24, x23, N2);
        E_mag = Mag3d(E);
        N1_mag = Mag3d(N1);
        N2_mag = Mag3d(N2);
        std::transform(E, E + 3, E, 
		std::bind2nd(std::divides<double>(), E_mag));
	std::transform(N1, N1 + 3, n1, 
		std::bind2nd(std::divides<double>(), N1_mag));
	std::transform(N2, N2 + 3, n2, 
		std::bind2nd(std::divides<double>(), N2_mag));
        std::transform(n1, n1 + 3, N1, 
		std::bind2nd(std::divides<double>(), N1_mag));
	std::transform(n2, n2 + 3, N2, 
		std::bind2nd(std::divides<double>(), N2_mag));

        double u1[3], u2[3], u3[3], u4[3];

	std::transform(N1, N1 + 3, u1,
                std::bind1st(std::multiplies<double>(), Dot3d(x14, E)));
	std::transform(N2, N2 + 3, u2,
                std::bind1st(std::multiplies<double>(), Dot3d(x24, E)));
	std::transform(u1, u1 + 3, u2, u3, std::plus<double>()); 
        std::transform(N1, N1 + 3, u1,
                std::bind1st(std::multiplies<double>(), -Dot3d(x13, E)));
        std::transform(N2, N2 + 3, u2,
                std::bind1st(std::multiplies<double>(), -Dot3d(x23, E)));
        std::transform(u1, u1 + 3, u2, u4, std::plus<double>());
	std::transform(N1, N1 + 3, u1, 
	        std::bind1st(std::multiplies<double>(), E_mag));
	std::transform(N2, N2 + 3, u2, 
	        std::bind1st(std::multiplies<double>(), E_mag));

	double bend_stiff = getBendStiff();
        double coeff = bend_stiff * sqr(E_mag) / (N1_mag + N2_mag);

        if (Dot3d(n1, n2) > 1.0 + 1.0e-10)
        {
	    std::cout << std::fixed << std::setprecision(14); 
	    std::cout << "t1 = "; 
	    std::for_each(Tri_normal(t1),Tri_normal(t1) + 3, DebugShow); 
	    std::cout << std::endl; 
	    std::cout << "n1 = ";
            std::for_each(n1,n1 + 3, DebugShow);
            std::cout << std::endl;
	    std::cout << "t2 = ";
            std::for_each(Tri_normal(t2),Tri_normal(t2) + 3, DebugShow);
            std::cout << std::endl;
	    std::cout << "n2 = ";
            std::for_each(n2,n2 + 3, DebugShow);
            std::cout << std::endl;
	    std::cout << "Dot3d(n1, n2) = "; 
            DebugShow(Dot3d(n1, n2)); 
            std::cout << std::endl; 
	    std::cout << "u1 = ";
            std::for_each(u1,u1 + 3, DebugShow);
            std::cout << std::endl;
            std::cout << "u2 = ";
            std::for_each(u2,u2 + 3, DebugShow);
            std::cout << std::endl;
            std::cout << "u3 = ";
            std::for_each(u3,u3 + 3, DebugShow);
            std::cout << std::endl;
            std::cout << "u4 = ";
            std::for_each(u4,u4 + 3, DebugShow);
            std::cout << std::endl;
            std::cout << "N1 = ";
            std::for_each(N1,N1 + 3, DebugShow);
            std::cout << std::endl;
            std::cout << "N2 = ";
            std::for_each(N2,N2 + 3, DebugShow);
            std::cout << std::endl;
            clean_up(0);
        }

        double sine_half_theta = sqrt(0.5 * std::max(0.0, 1.0 - Dot3d(n1, n2)));
        double tmp[3];

        Cross3d(n1, n2, tmp);
        if (Dot3d(tmp, E) < 0)
            sine_half_theta *= -1.0;
        coeff *= sine_half_theta;
	
	double bend_damp = getBendDamp(); 
        double dtheta = 0.0;

        dtheta = Dot3d(u1, p1->vel) + Dot3d(u2, p2->vel) +
                 Dot3d(u3, p3->vel) + Dot3d(u4, p4->vel);
        if (fabs(dtheta) < 1.0e-10) dtheta = 0.0;
// used for debugging        
//	double sum[3] = {0.0};

//	for (int i = 0; i < 3; i++) 
//	     sum[i] = u1[i] + u2[i] + u3[i] + u4[i]; 
//	int summ = Mag3d(sum); 
//	if (summ > 1.0e-6)
//	{
//	    std::cout << "u vector is too large\n"; 
//	    std::cout << "sum is: " << sum[0] << ' ' << sum[1] 
//		<< ' ' << sum[2];
//	    std::cout << "p1 is: " << Coords(p1)[0] << ' ' 
//		<< Coords(p1)[1] << ' ' << Coords(p1)[2] << std::endl;    
//	}
        

// used for debugging
//      if (dtheta > 0.0)
//        {
//            double sum[3];
//            for (int i = 0; i < 3; ++i)
//                sum[i] = u1[i] + u2[i] + u3[i] + u4[i];
//           printf("sum = %20.14f %20.14f %20.14f\n",sum[0],sum[1],sum[2]);
//            printf("p1->vel = %f %f %f\n",p1->vel[0],p1->vel[1],p1->vel[2]);
//            printf("dtheta = %20.14f\n", dtheta);
//            clean_up(0);
//        }     
        coeff += -bend_damp * E_mag * dtheta;

        // each tri_pair will be calculated twice
        std::transform(u1, u1 + 3, tmp, 
		std::bind1st(std::multiplies<double>(), coeff)); 
	std::transform(p1->force, p1->force + 3, tmp, 
		p1->force, std::plus<double>());
	std::transform(u2, u2 + 3, tmp, 
                std::bind1st(std::multiplies<double>(), coeff));
        std::transform(p2->force, p2->force + 3, tmp, 
                p2->force, std::plus<double>());
        std::transform(u3, u3 + 3, tmp, 
                std::bind1st(std::multiplies<double>(), coeff));
        std::transform(p3->force, p3->force + 3, tmp, 
                p3->force, std::plus<double>());
	std::transform(u4, u4 + 3, tmp, 
                std::bind1st(std::multiplies<double>(), coeff));
        std::transform(p4->force, p4->force + 3, tmp, 
                p4->force, std::plus<double>());
}       /* calculateBendingForce3d */

void BendingForce::calculateBendingForce3d2006(
        POINT* p1,
        TRI* t1,
        TRI* t2)
{
        if (Mag3d(Tri_normal(t1)) < MACH_EPS ||
            Mag3d(Tri_normal(t2)) < MACH_EPS) return;

        int index = Vertex_of_point(t1, p1);
        POINT *p3 = Point_of_tri(t1)[(index+1)%3];
        POINT *p4 = Point_of_tri(t1)[(index+2)%3];
        index = 3 - Vertex_of_point(t2, p3) - Vertex_of_point(t2, p4);
        POINT *p2 = Point_of_tri(t2)[index];

        double x13[3], x14[3], x23[3], x24[3], E[3];
        for (int i = 0; i < 3; ++i)
        {
            x13[i] = Coords(p1)[i] - Coords(p3)[i];
            x14[i] = Coords(p1)[i] - Coords(p4)[i];
            x23[i] = Coords(p2)[i] - Coords(p3)[i];
            x24[i] = Coords(p2)[i] - Coords(p4)[i];
            E[i] = Coords(p4)[i] - Coords(p3)[i];
        }
        double N1[3], N2[3], N3[3], N4[3];
        Cross3d(x13, x14, N1);
        Cross3d(x24, x23, N2);
        Cross3d(x23, x13, N3);
        Cross3d(x14, x24, N4);
        double N1_mag, N2_mag, N3_mag, N4_mag;
        N1_mag = Mag3d(N1);
        N2_mag = Mag3d(N2);
        N3_mag = Mag3d(N3);
        N4_mag = Mag3d(N4);
        double a1, a2, a3, a4;
        a1 = N2_mag / (N1_mag + N2_mag);
        a2 = 1 - a1;
        a3 = - N4_mag / (N3_mag + N4_mag);
        a4 = -1 - a3;
        double R[3];
        for (int i = 0; i < 3; ++i)
            R[i] = a1 * Coords(p1)[i] + a2 * Coords(p2)[i] +
                   a3 * Coords(p3)[i] + a4 * Coords(p4)[i];
        double E_mag, h1, h2;
        E_mag = Mag3d(E);
        h1 = fabs(Dot3d(x13, E) / E_mag);
        h1 = sqrt(Dot3d(x13, x13) - h1 * h1);
        h2 = fabs(Dot3d(x23, E) / E_mag);
        h2 = sqrt(Dot3d(x23, x23) - h2 * h2);

	double bend_stiff = getBendStiff(); 
        double lambda = 2.0*(h1 + h2)*E_mag*bend_stiff / (3.0*h1*h2*h1*h2);

        for (int i = 0; i < 3; ++i)
        {
            p1->force[i] += -lambda * a1 * R[i] * 0.5;
            p2->force[i] += -lambda * a2 * R[i] * 0.5;
            p3->force[i] += -lambda * a3 * R[i] * 0.5;
            p4->force[i] += -lambda * a4 * R[i] * 0.5;
        }
	if (fabs(Coords(p1)[0] -0.75) < 0.1)
	{
	    double R_mag = Mag3d(R);
	    int count = 0;
	    if (Coords(p1)[2] > 0.8) count++;
	    if (Coords(p2)[2] > 0.8) count++;
	    if (Coords(p3)[2] > 0.8) count++;
	    if (Coords(p4)[2] > 0.8) count++;
	    if (R_mag > 0 && count == 3)
	    {
	    printf("R = [%f %f %f], R_mag = %e, lambda = %e\n", R[0],R[1],R[2],R_mag,lambda);
	    printf("nor R = [%f %f %f]\n", R[0]/R_mag, R[1]/R_mag, R[2]/R_mag);
	    printf("h1 = %f, h2 = %f, E_mag = %f\n", h1, h2, E_mag);
	    printf("p1 = [%f %f %f]\n", Coords(p1)[0], Coords(p1)[1], Coords(p1)[2]); 
	    printf("p2 = [%f %f %f]\n", Coords(p2)[0], Coords(p2)[1], Coords(p2)[2]); 
	    printf("p3 = [%f %f %f]\n", Coords(p3)[0], Coords(p3)[1], Coords(p3)[2]); 
	    printf("p4 = [%f %f %f]\n", Coords(p4)[0], Coords(p4)[1], Coords(p4)[2]); 
	    }
	}
        if (0) //Mag3d(p1->force) == 0)
        {
        printf("f1 = [%f %f %f]\n", p1->force[0], p1->force[1], p1->force[2]);
        printf("f2 = [%f %f %f]\n", p2->force[0], p2->force[1], p2->force[2]);
        printf("f3 = [%f %f %f]\n", p3->force[0], p3->force[1], p3->force[2]);
        printf("f4 = [%f %f %f]\n", p4->force[0], p4->force[1], p4->force[2]);
        printf("a = [%f %f %f %f]\n\n", a1, a2, a3, a4);
        printf("h1 = %e, h2 = %e\n", h1, h2);
        }
}       /* end calculateBendingForce3d */

void BendingForce::clear_surf_point_force(SURFACE* surf)
{
        TRI* tri;
        surf_tri_loop(surf, tri)
        {
            for (int i = 0; i < 3; ++i)
            {
                POINT *p = Point_of_tri(tri)[i];
                for (int j = 0; j < 3; ++j)
                    p->force[j] = 0.0;
            }
        }
}       /* clear_surf_point_force */

