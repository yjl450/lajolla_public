#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyMetal& bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
        dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        return make_zero_spectrum();
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    Vector3 h = normalize(dir_in + dir_out);
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Spectrum fm = base_color + (1.0 - base_color) * pow(1.0 - abs(dot(h, dir_out)), 5);
    Real aspect = sqrt(1.0 - 0.9 * eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool));
    Real alpha_x, alpha_y, roughness;
    roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    constexpr Real alpha_min = 0.0001;
    alpha_x = max(alpha_min, roughness * roughness / aspect);
    alpha_y = max(alpha_min, roughness * roughness * aspect);
    Vector3 hl = to_local(frame, h);
    Real dm = 1.0 / (c_PI * alpha_x * alpha_y * pow(pow(hl.x / alpha_x, 2) + pow(hl.y / alpha_y, 2) + hl.z * hl.z, 2));
    Vector3 local_in = to_local(frame, dir_in);
    Vector3 local_out = to_local(frame, dir_out);
    Real lambda_in = lambda(dir_in, frame, alpha_x, alpha_y);
    Real lambda_out = lambda(dir_out, frame, alpha_x, alpha_y);
    Real g_in = 1.0 / (1.0 + lambda_in);
    Real g_out = 1.0 / (1.0 + lambda_out);
    Real gm = g_in * g_out;
    return fm * dm * gm / (4.0 * abs(dot(frame.n, dir_in)));
}

Real pdf_sample_bsdf_op::operator()(const DisneyMetal& bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
        dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        return 0;
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    Vector3 half_vector = normalize(dir_in + dir_out);
    Real n_dot_in = dot(frame.n, dir_in);
    Real n_dot_out = dot(frame.n, dir_out);
    Real n_dot_h = dot(frame.n, half_vector);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    Real aspect = sqrt(1.0 - 0.9 * eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool));
    constexpr Real alpha_min = 0.0001;
    Real alpha_x = max(alpha_min, roughness * roughness / aspect);
    Real alpha_y = max(alpha_min, roughness * roughness * aspect);
    Vector3 hl = to_local(frame, half_vector);
    Vector3 local_in = to_local(frame, dir_in);
    Real lambda_in = lambda(dir_in, frame, alpha_x, alpha_y);
    Real G = 1.0 / (1.0 + lambda_in);
    Real D = 1.0 / (c_PI * alpha_x * alpha_y * pow(pow(hl.x / alpha_x, 2) + pow(hl.y / alpha_y, 2) + hl.z * hl.z, 2));
    return (G * D) / (4.0 * n_dot_in);
}




std::optional<BSDFSampleRecord>
sample_bsdf_op::operator()(const DisneyMetal& bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0) {
        // No light below the surface
        return {};
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    Real aspect = sqrt(1.0 - 0.9 * eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool));
    constexpr Real alpha_min = 0.0001;
    Real alpha_x = max(alpha_min, roughness * roughness / aspect);
    Real alpha_y = max(alpha_min, roughness * roughness * aspect);
    Vector3 local_micro_normal =
        sample_visible_normals_anisotropic(to_local(frame, dir_in), alpha_x, alpha_y, rnd_param_uv);
    Vector3 half_vector = to_world(frame, local_micro_normal);
    Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
    return BSDFSampleRecord{
    reflected,
    Real(0), roughness };
}

TextureSpectrum get_texture_op::operator()(const DisneyMetal& bsdf) const {
    return bsdf.base_color;
}