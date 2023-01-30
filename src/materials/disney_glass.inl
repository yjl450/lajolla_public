#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyGlass &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    Vector3 h = dir_in + dir_out;
    h = normalize(h);
    Vector3 hl = to_local(frame, h);
    Real n_dot_in = dot(frame.n, dir_in);
    Real n_dot_out = dot(frame.n, dir_out);
    Real h_dot_in = dot(h, dir_in);
    Real h_dot_out = dot(h, dir_out);
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real aspect = sqrt(1.0 - 0.9 * eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool));
    Real alpha_x, alpha_y, roughness;
    roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    constexpr Real alpha_min = 0.0001;
    alpha_x = max(alpha_min, roughness * roughness / aspect);
    alpha_y = max(alpha_min, roughness * roughness * aspect);
    Real dg = c_INVPI * 1.0 / (alpha_x * alpha_y * pow(pow(hl.x / alpha_x, 2) + pow(hl.y / alpha_y, 2) + hl.z * hl.z, 2));
    Real lambda_in = sqrt(1.0 + (pow(dir_in.x * alpha_x, 2) + pow(dir_in.y * alpha_y, 2)) / (dir_in.z * dir_in.z)) / 2.0 - 0.5;
    Real lambda_out = sqrt(1.0 + (pow(dir_out.x * alpha_x, 2) + pow(dir_out.y * alpha_y, 2)) / (dir_out.z * dir_out.z)) / 2.0 - 0.5;
    Real g_in = 1.0 / (1.0 + lambda_in);
    Real g_out = 1.0 / (1.0 + lambda_out);
    Real gg = g_in * g_out;
    Real rs = (h_dot_in - bsdf.eta * h_dot_out) / (h_dot_in + bsdf.eta * h_dot_out);
    Real rp = (bsdf.eta * h_dot_in - h_dot_out) / (bsdf.eta * h_dot_in + h_dot_out);
    Real fg = 0.5 * (rs * rs + rp * rp);
    if (n_dot_in * n_dot_out > 0) {
        return base_color * fg * dg * gg / (4.0*abs(n_dot_in));
    }
    else {
        return sqrt(base_color) * (1.0 - fg) * dg * gg * abs(h_dot_out * h_dot_in) / (abs(n_dot_in) * pow(h_dot_in + bsdf.eta * h_dot_out, 2));
    }
}

Real pdf_sample_bsdf_op::operator()(const DisneyGlass &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    return fmax(dot(frame.n, dir_out), Real(0)) / c_PI;
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyGlass &bsdf) const {
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!

    return BSDFSampleRecord{
        to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
        bsdf.eta /* eta */, eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool) };
}

TextureSpectrum get_texture_op::operator()(const DisneyGlass &bsdf) const {
    return bsdf.base_color;
}
