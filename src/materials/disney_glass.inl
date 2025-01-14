#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyGlass& bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
        dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!

    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    assert(eta > 0);
    Vector3 h;
    if (reflect) {
        h = normalize(dir_in + dir_out);
    }
    else {
        h = normalize(dir_in + dir_out * eta);
    }
    if (dot(h, frame.n) < 0) {
        h = -h;
    }

    Vector3 hl = to_local(frame, h);
    Real n_dot_in = dot(frame.n, dir_in);
    Real n_dot_out = dot(frame.n, dir_out);
    Real h_dot_in = dot(h, dir_in);
    Real h_dot_out = dot(h, dir_out);
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real aspect = sqrt(1.0 - 0.9 * eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool));
    Real alpha_x, alpha_y, roughness;
    roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    constexpr Real alpha_min = 0.0001;
    alpha_x = max(alpha_min, roughness * roughness / aspect);
    alpha_y = max(alpha_min, roughness * roughness * aspect);
    Real dg = 1.0 / (c_PI * alpha_x * alpha_y * pow(pow(hl.x / alpha_x, 2) + pow(hl.y / alpha_y, 2) + hl.z * hl.z, 2));
    Vector3 local_in = to_local(frame, dir_in);
    Vector3 local_out = to_local(frame, dir_out);
    Real lambda_in = sqrt(1.0 + (pow(local_in.x * alpha_x, 2) + pow(local_in.y * alpha_y, 2)) / (local_in.z * local_in.z)) / 2.0 - 0.5;
    Real lambda_out = sqrt(1.0 + (pow(local_out.x * alpha_x, 2) + pow(local_out.y * alpha_y, 2)) / (local_out.z * local_out.z)) / 2.0 - 0.5;
    Real g_in = 1.0 / (1.0 + lambda_in);
    Real g_out = 1.0 / (1.0 + lambda_out);
    Real gg = g_in * g_out;
    Real fg = fresnel_dielectric(abs(h_dot_in), abs(h_dot_out), eta);
    if (reflect) {
        return base_color * fg * dg * gg / (4.0 * abs(n_dot_in));
    }
    else {
        return sqrt(base_color) * (1.0 - fg) * dg * gg * abs(h_dot_out * h_dot_in) / (abs(n_dot_in) * pow(h_dot_in + eta * h_dot_out, 2));
    }
}

Real pdf_sample_bsdf_op::operator()(const DisneyGlass& bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
        dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    // If we are going into the surface, then we use normal eta
    // (internal/external), otherwise we use external/internal.
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    assert(eta > 0);

    Vector3 half_vector;
    if (reflect) {
        half_vector = normalize(dir_in + dir_out);
    }
    else {
        half_vector = normalize(dir_in + dir_out * eta);
    }
    if (dot(half_vector, frame.n) < 0) {
        half_vector = -half_vector;
    }

    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    roughness = std::clamp(roughness, Real(0.01), Real(1));

    // We sample the visible normals, also we use F to determine
    // whether to sample reflection or refraction
    // so PDF ~ F * D * G_in for reflection, PDF ~ (1 - F) * D * G_in for refraction.
    Real h_dot_in = dot(half_vector, dir_in);
    Real h_dot_out = dot(half_vector, dir_out);
    Real F = fresnel_dielectric(abs(h_dot_in), abs(h_dot_out), eta);
    constexpr Real alpha_min = 0.0001;
    Real aspect = sqrt(1.0 - 0.9 * eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool));
    Real alpha_x = max(alpha_min, roughness * roughness / aspect);
    Real alpha_y = max(alpha_min, roughness * roughness * aspect);
    Vector3 hl = to_local(frame, half_vector);
    Real D = 1.0 / (c_PI * alpha_x * alpha_y * pow(pow(hl.x / alpha_x, 2) + pow(hl.y / alpha_y, 2) + hl.z * hl.z, 2));
    //Real G_in = smith_masking_gtr2(to_local(frame, dir_in), roughness);
    Vector3 local_in = to_local(frame, dir_in);
    Real lambda_in = sqrt(1.0 + (pow(local_in.x * alpha_x, 2) + pow(local_in.y * alpha_y, 2)) / (local_in.z * local_in.z)) / 2.0 - 0.5;
    Real G_in = 1.0 / (1.0 + lambda_in);
    if (reflect) {
        return (F * D * G_in) / (4 * fabs(dot(frame.n, dir_in)));
    }
    else {
        Real h_dot_out = dot(half_vector, dir_out);
        Real sqrt_denom = h_dot_in + eta * h_dot_out;
        Real dh_dout = eta * eta * h_dot_out / (sqrt_denom * sqrt_denom);
        return (1 - F) * D * G_in * fabs(dh_dout * h_dot_in / dot(frame.n, dir_in));
    }
}

std::optional<BSDFSampleRecord>
sample_bsdf_op::operator()(const DisneyGlass& bsdf) const {
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    assert(eta > 0);
    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    Vector3 local_dir_in = to_local(frame, dir_in);
    constexpr Real alpha_min = 0.0001;
    Real aspect = sqrt(1.0 - 0.9 * eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool));
    Real alpha_x = max(alpha_min, roughness * roughness / aspect);
    Real alpha_y = max(alpha_min, roughness * roughness * aspect);
    Vector3 local_micro_normal =
        sample_visible_normals_anisotropic(local_dir_in, alpha_x, alpha_y, rnd_param_uv);

    Vector3 half_vector = to_world(frame, local_micro_normal);
    if (dot(half_vector, frame.n) < 0) {
        half_vector = -half_vector;
    }
    Real h_dot_in = dot(half_vector, dir_in);
    Real F = fresnel_dielectric(h_dot_in, eta);

    if (rnd_param_w <= F) {
        Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
        return BSDFSampleRecord{ reflected, Real(0) /* eta */, roughness };
    }
    else {
        Real h_dot_out_sq = 1 - (1 - h_dot_in * h_dot_in) / (eta * eta);
        if (h_dot_out_sq <= 0) {
            return {};
        }
        if (h_dot_in < 0) {
            half_vector = -half_vector;
        }
        Real h_dot_out = sqrt(h_dot_out_sq);
        Vector3 refracted = -dir_in / eta + (fabs(h_dot_in) / eta - h_dot_out) * half_vector;
        return BSDFSampleRecord{ refracted, eta, roughness };
    }
}

TextureSpectrum get_texture_op::operator()(const DisneyGlass& bsdf) const {
    return bsdf.base_color;
}