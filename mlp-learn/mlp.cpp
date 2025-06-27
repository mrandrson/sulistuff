#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <functional>

void load_csv(const std::string& filename, std::vector<double>& x, std::vector<double>& y) {
    std::ifstream file(filename);
    std::string line;
    std::getline(file, line); // skip header

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string xs, ys;
        std::getline(ss, xs, ',');
        std::getline(ss, ys, ',');
        x.push_back(std::stod(xs));
        y.push_back(std::stod(ys));
    }
}

#define MAX_ITER 40000
#define ALPHA 0.005
#define BETA1 0.9
#define BETA2 0.999
#define EPSILON 1e-8
#define BATCH_SIZE 25
#define DIFF_H 1e-6

const int H1 = 30;  
const int H2 = 60;  
const int H3 = 30;  

double mlp_model(double x, const std::vector<double>& params) {
    std::vector<double> h1(H1), h2(H2), h3(H3);
    int offset = 0;

    // Layer 1: Input -> H1
    for (int i = 0; i < H1; ++i)
        h1[i] = std::tanh(params[offset + i] * x + params[offset + H1 + i]); // ReLU(W1 * x + b1)
    offset += 2 * H1;

    // Layer 2: H1 -> H2
    for (int i = 0; i < H2; ++i) {
        double sum = 0.0;
        for (int j = 0; j < H1; ++j)
            sum += params[offset + i * H1 + j] * h1[j];
        h2[i] = std::max(0.0, sum + params[offset + H1 * H2 + i]); // tanh(W2 * h1 + b2)
    }
    offset += H2 * H1 + H2;

    // Layer 3: H2 -> H3
    for (int i = 0; i < H3; ++i) {
        double sum = 0.0;
        for (int j = 0; j < H2; ++j)
            sum += params[offset + i * H2 + j] * h2[j];
        h3[i] = std::tanh(sum + params[offset + H2 * H3 + i]); // ReLU(W3 * h2 + b3)
    }
    offset += H3 * H2 + H3;

    // Output layer: H3 -> Output
    double out = 0.0;
    for (int i = 0; i < H3; ++i)
        out += params[offset + i] * h3[i]; // W4 * h3
    out += params[offset + H3]; // b4

    return out;
}


double compute_error(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& params, std::function<double(double, const std::vector<double>&)> model) {
    double error = 0.0;
    for (int i = 0; i < x.size(); ++i) {
        double f = model(x[i], params);
        error += std::pow(f - y[i], 2);
    }

    std::cout << "[compute_error] Error: " << error << std::endl;
    return error;
}

/*
void compute_gradient(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& params, std::vector<double>& grad, int batch_size,std::function<double(double, const std::vector<double>&)> model) {
    grad.resize(params.size());
    std::vector<double> params_perturbed = params;

    double base_error = 0.0;
    if (batch_size == x.size()) base_error = compute_error(x, y, params, model);

    for (int i = 0; i < params.size(); ++i) {
        params_perturbed[i] += DIFF_H;

        double perturbed_error;
        if (batch_size == x.size()) {
            perturbed_error = compute_error(x, y, params_perturbed, model);
        } else {
            perturbed_error = 0.0;
            base_error = 0.0;
            for (int j = 0; j < batch_size; ++j) {
                int idx = rand() % x.size();
                double fxp = model(x[idx], params_perturbed);
                double fx = model(x[idx], params);
                perturbed_error += std::pow(fxp - y[idx], 2);
                base_error += std::pow(fx - y[idx], 2);
            }
            perturbed_error /= batch_size;
            base_error /= batch_size;
        }

        grad[i] = (perturbed_error - base_error) / DIFF_H;
        params_perturbed[i] = params[i];
    }
}
*/

inline double d_tanh(double x) {
    double t = std::tanh(x);
    return 1.0 - t * t;
}

inline double d_relu(double x) {
    return x > 0.0 ? 1.0 : 0.0;
}

void compute_gradient(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& params, std::vector<double>& grad, int batch_size, std::function<double(double, const std::vector<double>&)> model) {
    grad.assign(params.size(), 0.0);
    std::vector<double> local_grad(params.size(), 0.0);

    // Helper functions
    auto d_tanh = [](double z) { double t = std::tanh(z); return 1.0 - t * t; };
    auto d_relu = [](double z) { return z > 0.0 ? 1.0 : 0.0; };

    for (int b = 0; b < batch_size; ++b) {
        int idx = rand() % x.size();
        double input = x[idx];
        double target = y[idx];

        // === FORWARD PASS ===
        int offset = 0;
        std::vector<double> z1(H1), a1(H1);
        std::vector<double> z2(H2), a2(H2);
        std::vector<double> z3(H3), a3(H3);

        for (int i = 0; i < H1; ++i) {
            z1[i] = params[offset + i] * input + params[offset + H1 + i];
            a1[i] = std::tanh(z1[i]);
        }
        offset += 2 * H1;

        for (int i = 0; i < H2; ++i) {
            z2[i] = 0.0;
            for (int j = 0; j < H1; ++j)
                z2[i] += params[offset + i * H1 + j] * a1[j];
            z2[i] += params[offset + H1 * H2 + i];
            a2[i] = std::max(0.0, z2[i]);
        }
        offset += H2 * H1 + H2;

        for (int i = 0; i < H3; ++i) {
            z3[i] = 0.0;
            for (int j = 0; j < H2; ++j)
                z3[i] += params[offset + i * H2 + j] * a2[j];
            z3[i] += params[offset + H2 * H3 + i];
            a3[i] = std::tanh(z3[i]);
        }
        offset += H3 * H2 + H3;

        double out = 0.0;
        for (int i = 0; i < H3; ++i)
            out += params[offset + i] * a3[i];
        out += params[offset + H3];

        double loss_grad = 2.0 * (out - target);  // dL/dyhat

        // === BACKWARD PASS ===
        std::vector<double> d3(H3), d2(H2), d1(H1);
        int out_offset = offset;

        // Output layer
        for (int i = 0; i < H3; ++i) {
            grad[out_offset + i] += loss_grad * a3[i];
            d3[i] = loss_grad * params[out_offset + i] * d_tanh(z3[i]);
        }
        grad[out_offset + H3] += loss_grad;

        offset -= (H3 * H2 + H3);

        // Layer 3
        for (int i = 0; i < H3; ++i) {
            for (int j = 0; j < H2; ++j)
                grad[offset + i * H2 + j] += d3[i] * a2[j];
            grad[offset + H3 * H2 + i] += d3[i];
        }

        for (int j = 0; j < H2; ++j)
            for (int i = 0; i < H3; ++i)
                d2[j] += d3[i] * params[offset + i * H2 + j] * d_relu(z2[j]);

        offset -= (H2 * H1 + H2);

        // Layer 2
        for (int i = 0; i < H2; ++i) {
            for (int j = 0; j < H1; ++j)
                grad[offset + i * H1 + j] += d2[i] * a1[j];
            grad[offset + H1 * H2 + i] += d2[i];
        }

        for (int j = 0; j < H1; ++j)
            for (int i = 0; i < H2; ++i)
                d1[j] += d2[i] * params[offset + i * H1 + j] * d_tanh(z1[j]);

        offset -= 2 * H1;

        // Layer 1
        for (int i = 0; i < H1; ++i) {
            grad[offset + i] += d1[i] * input;
            grad[offset + H1 + i] += d1[i];
        }
    }

    // Average over batch
    for (auto& g : grad) g /= batch_size;
}

void fit(const std::vector<double>& x, const std::vector<double>& y, std::vector<double>& params,
         std::function<double(double, const std::vector<double>&)> model) {
    std::ofstream file("results.csv");
    file << "Iter,Error";
    for (int i = 0; i < params.size(); ++i) file << ",p" << i;
    file << "\n";

    std::vector<double> m(params.size(), 0.0), v(params.size(), 0.0), grad(params.size(), 0.0);

    for (int t = 1; t <= MAX_ITER; ++t) {
        compute_gradient(x, y, params, grad, BATCH_SIZE, model);

        for (int i = 0; i < params.size(); ++i) {
            m[i] = BETA1 * m[i] + (1 - BETA1) * grad[i];
            v[i] = BETA2 * v[i] + (1 - BETA2) * grad[i] * grad[i];

            double m_hat = m[i] / (1 - std::pow(BETA1, t));
            double v_hat = v[i] / (1 - std::pow(BETA2, t));

            params[i] -= ALPHA * m_hat / (std::sqrt(v_hat) + EPSILON);
        }

        double err = compute_error(x, y, params, model);
        file << t << "," << err;
        for (const auto& p : params) file << "," << p;
        file << "\n";

        if (err < 1e-6) break;
    }
    file.close();
}

/*
int main() {
    // srand(42);

    std::vector<double> x, y;
    load_csv("synthetic_data.csv", x, y);

    // std::vector<double> params(2 * H1 + H2 * H1 + H2 + H2 + 1);  // For 2 Layer MLP
    std::vector<double> params(2*H1 + H2*H1 + H2 + H3*H2 + H3 + H3 + 1); // For 3 Layer MLP
    for (auto& p : params) p = ((rand() % 2000) / 1000.0 - 1.0);  // [-1, 1] 

    fit(x, y, params, mlp_model);

    return 0;
}
*/

int main() {
    std::vector<double> x, y;
    load_csv("synthetic_data.csv", x, y);

    // Normalize x to [0, 1]
    double x_min = *std::min_element(x.begin(), x.end());
    double x_max = *std::max_element(x.begin(), x.end());
    for (auto& xi : x)
        xi = (xi - x_min) / (x_max - x_min);

    // Normalize y to [0, 1]
    double y_min = *std::min_element(y.begin(), y.end());
    double y_max = *std::max_element(y.begin(), y.end());
    for (auto& yi : y)
        yi = (yi - y_min) / (y_max - y_min);

    // Truncate to first quarter of entries
    size_t quarter_size = x.size();
    x.resize(quarter_size);
    y.resize(quarter_size);

    // Initialize parameters
    std::vector<double> params(2*H1 + H2*H1 + H2 + H3*H2 + H3 + H3 + 1); // For 3 Layer MLP
    for (auto& p : params)
        p = ((rand() % 2000) / 1000.0 - 1.0);  // Uniform [-1, 1]

    fit(x, y, params, mlp_model);

    return 0;
}



