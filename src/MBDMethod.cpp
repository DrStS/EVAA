
#include "MBDMethod.h"

namespace EVAA {

#if MIGHT_BE_USEFUL

void create_basis_car(floatEVAA* qc, floatEVAA* basis_c)
{
    const MKL_INT n = 4, incx = 1;
    floatEVAA q_norm_inv = 2 / (floatEVAA)cblas_ddot(n, qc, incx, qc, incx);  // 2 / (||q||_2)^2
    basis_c[0] = 1 - q_norm_inv * (qc[1] * qc[1] + qc[2] * qc[2]);
    basis_c[1] = q_norm_inv * (qc[0] * qc[1] + qc[2] * qc[3]);
    basis_c[2] = q_norm_inv * (qc[0] * qc[2] - qc[1] * qc[3]);
    basis_c[3] = q_norm_inv * (qc[0] * qc[1] - qc[2] * qc[3]);
    basis_c[4] = 1 - q_norm_inv * (qc[0] * qc[0] + qc[2] * qc[2]);
    basis_c[5] = q_norm_inv * (qc[1] * qc[2] + qc[0] * qc[3]);
    basis_c[6] = q_norm_inv * (qc[0] * qc[2] + qc[1] * qc[3]);
    basis_c[7] = q_norm_inv * (qc[1] * qc[2] - qc[0] * qc[3]);
    basis_c[8] = 1 - q_norm_inv * (qc[0] * qc[0] + qc[1] * qc[1]);
}

void get_tilda(floatEVAA* x, floatEVAA* x_tilda)
{
    // x - vector 3 elems (floatEVAA* x); x - tilda matrix 3x3 (floatEVAA* x_tilda)
    // given y: x_tilda*y=cross(x,y) (Stoneking, page 3 bottom)
    x_tilda[0] = 0;
    x_tilda[1] = -x[2];
    x_tilda[2] = x[1];
    x_tilda[3] = x[2];
    x_tilda[4] = 0;
    x_tilda[5] = -x[0];
    x_tilda[6] = -x[1];
    x_tilda[7] = x[0];
    x_tilda[8] = 0;
}

void get_tilda_r(floatEVAA* r, const int Constants::DIM, const int num_wheels, floatEVAA* r_tilda)
{
    // r - matrix 3x4 elems (floatEVAA* r); r - tilda matrix 3x12 (floatEVAA* x_tilda)
    // given y: x_tilda*y=cross(x,y) (Stoneking, page 3 bottom)
    for (int i = 0; i < num_wheels; i++) {
        r_tilda[0 + i * Constants::DIM] = 0;
        r_tilda[1 + i * Constants::DIM] = -r[2 * num_wheels + i];
        r_tilda[2 + i * Constants::DIM] = r[num_wheels + i];
        r_tilda[Constants::DIM * num_wheels + i * Constants::DIM] = r[2 * num_wheels + i];
        r_tilda[Constants::DIM * num_wheels + 1 + i * Constants::DIM] = 0;
        r_tilda[Constants::DIM * num_wheels + 2 + i * Constants::DIM] = -r[i];
        r_tilda[2 * Constants::DIM * num_wheels + i * Constants::DIM] = -r[num_wheels + i];
        r_tilda[2 * Constants::DIM * num_wheels + 1 + i * Constants::DIM] = r[i];
        r_tilda[2 * Constants::DIM * num_wheels + 2 + i * Constants::DIM] = 0;
    }
}

void C_cos_transf(floatEVAA* Y, floatEVAA* C_transf)
{
    // This function defines the cosine transformation matrix, used in changing coordinates from
    // a frame to another.The X and Y matrices represent basis in the X and, respectively,
    //    Y frame.x = C * y
    // Cosinus transformation to N base (eye(3)). X basis is the identity implicitly; Y basis is
    // the basis_car
    //    C = [ Y_1 /||Y_1||, Y_2 /||Y_2||, Y_3 /||Y_3||]; Y_i - columns of Y

    floatEVAA Y_norms_on_columns;
    for (auto i = 0; i < 3; ++i) {
        Y_norms_on_columns = 1. / (floatEVAA)(cblas_dnrm2(3, Y + i, 3));
        Y[i] /= Y_norms_on_columns;
        Y[i + 3] /= Y_norms_on_columns;
        Y[i + 6] /= Y_norms_on_columns;
    }
}

void get_quaternion_derivative(floatEVAA* qc, floatEVAA* Qc)
{
    // input: qc - the quaternion, 4x1 vector
    // output: Qc - the matrix of quaternion's derivative; 4x3 matrix in vector representation
    with row major Qc[0] = 0.5 * qc[3];
    Qc[1] = -0.5 * qc[2];
    Qc[2] = 0.5 * qc[1];

    Qc[3] = 0.5 * qc[2];
    Qc[4] = 0.5 * qc[3];
    Qc[5] = -0.5 * qc[0];

    Qc[6] = -0.5 * qc[1];
    Qc[7] = 0.5 * qc[0];
    Qc[8] = 0.5 * qc[3];

    Qc[9] = -0.5 * qc[0];
    Qc[10] = -0.5 * qc[1];
    Qc[11] = -0.5 * qc[2];
}
#endif

}  // namespace EVAA
