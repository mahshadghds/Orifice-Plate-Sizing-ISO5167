#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>
#include <cstdlib>

using namespace std;

double calc_c(double Re, double beta, double D, string tap_type) {
    double L1 = 0.0;
    double L2_prime = 0.0;

    if (tap_type == "Flange" || tap_type == "flange") {
        L1 = 25.4 / D;
        L2_prime = 25.4 / D;
    } else if (tap_type == "D-D2" || tap_type == "d-d2") {
        L1 = 1.0;
        L2_prime = 0.47;
    } else { //for corner tapping 
        L1 = 0.0;
        L2_prime = 0.0;
    }

    double M2_prime = (2.0 * L2_prime) / (1.0 - beta);
    double A = pow((19000.0 * beta) / Re, 0.8);

    double term1 = 0.5961;
    double term2 = 0.0261 * pow(beta, 2);
    double term3 = -0.216 * pow(beta, 8);
    double term4 = 0.000521 * pow((1000000.0 * beta) / Re, 0.7);
    double term5 = (0.0188 + 0.0063 * A) * pow(beta, 3.5) * pow(1000000.0 / Re, 0.3);
    double term6 = (0.043 + 0.080 * exp(-10.0 * L1) - 0.123 * exp(-7.0 * L1)) * (1.0 - 0.11 * A) * (pow(beta, 4) / (1.0 - pow(beta, 4)));
    double term7 = -0.031 * (M2_prime - 0.8 * pow(M2_prime, 1.1)) * pow(beta, 1.3);

    double C = term1 + term2 + term3 + term4 + term5 + term6 + term7;

    if (D < 71.12) {
        C += 0.011 * (0.75 - beta) * (2.8 - (D / 25.4));
    }

    return C;
}

double calc_epsilon(string phase, double beta, double p1_bar_g, double deltaP_mbar, double kappa) {
    if (phase == "Liquid" || phase == "liquid") {
        return 1.0;
    }
    
    double p1_abs_mbar = (p1_bar_g * 1000.0) + 1013.25;
    double p2_abs_mbar = p1_abs_mbar - deltaP_mbar;
    double p_ratio = p2_abs_mbar / p1_abs_mbar;
    
    if (p_ratio < 0.75) {
        cerr << "ERROR: p2/p1 limit\n";
        exit(1);
    }
    
    double term1 = 0.351 + 0.256 * pow(beta, 4) + 0.93 * pow(beta, 8);
    double term2 = 1.0 - pow(p_ratio, 1.0 / kappa);
    
    return 1.0 - (term1 * term2);
}

void check_limits(double D, double qv, double rho, double mu_cp, double d, double beta, double Re, double C) {
    if (D < 50 || D > 1000) {
        cerr << "ERROR:D limit\n";
        exit(1);
    }
    if (qv <= 0 || rho <= 0 || mu_cp <= 0) {
        cerr << "ERROR:Negative or zero inputs\n";
        exit(1);
    }
    if (d < 12.5 || d > 1000) {
        cerr << "ERROR:d MUST BE bigger THAN 12.5mm\n";
        exit(1);
    }
    if (beta < 0.1 || beta > 0.75) {
        cerr << "ERROR:Beta (d/D) limit\n";
        exit(1);
    }
    
    double Re_min;
    if (beta <= 0.56) {
        Re_min = 5000;
    } else {
        Re_min = 16000;
    }
    
    if (Re < Re_min) {
        cerr << "ERROR:Re_d limit\n";
        exit(1);
    }
    if (C < 0.55 || C > 0.65) {
        cerr << "ERROR:C limit\n";
    }
}

int main() {
    string line_no, fluid_name, phase, tap_type, material;
    double D, qv, T_normal, T_min, T_max, P_normal, P_min, P_max, rho, mu_cp, kappa,deltaP_mbar;

    cout << "ENTER D: \n";
    cin >> D;
    cout << "ENTER phase (Liquid/Gas): \n";
    cin >> phase;
    cout << "ENTER qv: \n";
    cin >> qv;
    cout << "ENTER P normal: \n";
    cin >> P_normal;
    cout << "ENTER rho: \n";
    cin >> rho;
    cout << "ENTER mu: \n";
    cin >> mu_cp;
    cout << "ENTER kappa (Cp/Cv): \n";
    cin >> kappa;
    cout << "ENTER tap_type: \n";
    cin >> tap_type;
        cout << "ENTER delta p (mbar): \n";
    cin >> deltaP_mbar;
   
    double D_m = D / 1000.0;
    double mu_Pa_s = mu_cp / 1000.0;
    
    double qm_kg_h=qv;
    double qm = qm_kg_h/3600.0;

 
    
    double Re = (4.0 * qm) / (M_PI * D_m * mu_Pa_s);
    
   
    double deltaP = deltaP_mbar * 100.0;
    
    double d = D * 0.5;
    double beta = d / D;
    double C = 0.6;
    double epsilon = 1.0;
    double tolerance = 1e-6;
    double diff = 1.0;
    int max_iter = 100;
    int iter_count = 0;

    while (diff > tolerance && iter_count < max_iter) {
        C = calc_c(Re, beta, D, tap_type);
        epsilon = calc_epsilon(phase, beta, P_normal, deltaP_mbar, kappa);
        
        double term_sqrt = sqrt(1.0 - pow(beta, 4));
        double den = M_PI * C * epsilon * sqrt(2.0 * rho * deltaP);
        double d_new = 1000.0 * sqrt((4.0 * qm * term_sqrt) / den);
        
        diff = abs(d_new - d);
        d = d_new;
        beta = d / D;
        iter_count++;
    }
    
    check_limits(D, qv, rho, mu_cp, d, beta, Re, C);
    
    double beta_4 = pow(beta, 4);
    double C_sq = pow(C, 2);
    
    double num = sqrt(1.0 - beta_4 * (1.0 - C_sq)) - C * pow(beta, 2);
    double den = sqrt(1.0 - beta_4 * (1.0 - C_sq)) + C * pow(beta, 2);
    double ratio = num / den;
    double Dsigma_Pa = ratio * deltaP;
    double pressure_loss = Dsigma_Pa / 100.0;
    
  
    
    double UC;
    if (beta < 0.2) {
        UC = (0.7 - beta);
    } else if (beta <= 0.6) {
        UC = 0.005;
    } else {
        UC = (1.667 * beta - 0.5) / 100.0;
    }
    
    if (D < 71.12) {
        UC += 0.9 * (0.75 - beta) * (2.8 - D / 25.4) / 100.0;
    }
    if (beta > 0.5 && Re < 10000) {
        UC += 0.005;
    }
    double uncertainty = UC * 100.0;

    cout << "\nqm_kg_h: " << fixed << setprecision(3) << qm_kg_h << "\n";
    cout << "d(mm): " << fixed << setprecision(2) << d << "\n";
    cout << "d(inch): " << fixed << setprecision(3) << d/25.4 << "\n";
    cout << "beta: " << fixed << setprecision(4) << beta << "\n";
    cout << "Re: " << fixed << setprecision(0) << Re << "\n";
    cout << "C: " << fixed << setprecision(4) << C << "\n";
    cout << "epsilon: " << fixed << setprecision(4) << epsilon << "\n";
    cout << "deltaP(mbar): " << fixed << setprecision(0) << deltaP_mbar << "\n";
    cout << "pressure_loss(mbar): " << fixed << setprecision(2) << pressure_loss << "\n";
    cout << "uncertainty: " << fixed << setprecision(2) << uncertainty << "\n";
    cout << "Total Iterations: " << iter_count << setprecision(2)  << "\n";

    return 0;
}
