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
    Real lambda_in = sqrt(1.0 + (pow(local_in.x * alpha_x, 2) + pow(local_in.y * alpha_y, 2)) / (local_in.z * local_in.z)) / 2.0 - 0.5;
    Real lambda_out = sqrt(1.0 + (pow(local_out.x * alpha_x, 2) + pow(local_out.y * alpha_y, 2)) / (local_out.z * local_out.z)) / 2.0 - 0.5;
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
    Real lambda_in = sqrt(1.0 + (pow(local_in.x * alpha_x, 2) + pow(local_in.y * alpha_y, 2)) / (local_in.z * local_in.z)) / 2.0 - 0.5;
    Real G = 1.0 / (1.0 + lambda_in);
    Real D = 1.0 / (c_PI * alpha_x * alpha_y * pow(pow(hl.x / alpha_x, 2) + pow(hl.y / alpha_y, 2) + hl.z * hl.z, 2));
    return (G * D) / (4.0 * n_dot_in);
}


inline Vector3 sample_visible_normals_m(const Vector3& local_dir_in, Real alpha_x, Real alpha_y, const Vector2& rnd_param) {
    if (local_dir_in.z < 0) {
        return -sample_visible_normals_m(-local_dir_in, alpha_x, alpha_y, rnd_param);
    }
    Vector3 hemi_dir_in = normalize(
        Vector3{ alpha_x * local_dir_in.x, alpha_y * local_dir_in.y, local_dir_in.z });

    Real r = sqrt(rnd_param.x);
    Real phi = 2 * c_PI * rnd_param.y;
    Real t1 = r * cos(phi);
    Real t2 = r * sin(phi);
    Real s = (1 + hemi_dir_in.z) / 2;
    t2 = (1 - s) * sqrt(1 - t1 * t1) + s * t2;
    Vector3 disk_N{ t1, t2, sqrt(max(Real(0), 1 - t1 * t1 - t2 * t2)) };
    Frame hemi_frame(hemi_dir_in);
    Vector3 hemi_N = to_world(hemi_frame, disk_N);
    return normalize(Vector3{ alpha_x * hemi_N.x, alpha_y * hemi_N.y, max(Real(0), hemi_N.z) });
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
        sample_visible_normals_m(to_local(frame, dir_in), alpha_x, alpha_y, rnd_param_uv);
    Vector3 half_vector = to_world(frame, local_micro_normal);
    Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
    return BSDFSampleRecord{
    reflected,
    Real(0), roughness };
}

TextureSpectrum get_texture_op::operator()(const DisneyMetal& bsdf) const {
    return bsdf.base_color;
}
