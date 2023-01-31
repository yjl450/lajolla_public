#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyClearcoat& bsdf) const {
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
    constexpr Real eta = 1.5;
    Real r0 = pow(eta - 1.0, 2) / pow(eta + 1.0, 2);
    Real fc = r0 + (1.0 - r0) * pow(1.0 - abs(dot(h, dir_out)), 5);
    Real gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real alpha_g2 = pow((1.0 - gloss) * 0.1 + gloss * 0.001, 2);
    Vector3 hl = to_local(frame, h);
    Real dc = (alpha_g2 - 1.0) * c_INVPI / (log(alpha_g2) * (1.0 + (alpha_g2 - 1) * (hl.z * hl.z)));
    Real gc = (1.0 / (1.0 + lambda(dir_in, frame, 0.25, 0.25))) * (1.0 / (1.0 + lambda(dir_out, frame, 0.25, 0.25)));
    Real f_clearcoat = fc * dc * gc / (4.0 * abs(dot(frame.n, dir_in)));
    return make_const_spectrum(f_clearcoat);
}

Real pdf_sample_bsdf_op::operator()(const DisneyClearcoat& bsdf) const {
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
    // Homework 1: implement this!
    //return fmax(dot(frame.n, dir_out), Real(0)) / c_PI;
    Vector3 h = normalize(dir_in + dir_out);
    Real gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real alpha_g2 = pow((1.0 - gloss) * 0.1 + gloss * 0.001, 2);
    Vector3 hl = to_local(frame, h);
    Real dc = (alpha_g2 - 1.0) * c_INVPI / (log(alpha_g2) * (1.0 + (alpha_g2 - 1) * (hl.z * hl.z)));
    return dc * abs(dot(frame.n, h)) / (4.0 * abs(dot(h, dir_out)));
}

std::optional<BSDFSampleRecord>
sample_bsdf_op::operator()(const DisneyClearcoat& bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0) {
        // No light below the surface
        return {};
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    Real gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real alpha_g2 = pow((1.0 - gloss) * 0.1 + gloss * 0.001, 2);
    Real cos_ele = sqrt((1.0 - pow(alpha_g2, 1 - rnd_param_uv.x)) / (1.0 - alpha_g2));
    Real azi = 2.0 * c_PI * rnd_param_uv.y;
    Vector3 micronormal;
    micronormal.x = sin(acos(cos_ele)) * cos(azi);
    micronormal.y = sin(acos(cos_ele)) * sin(azi);
    micronormal.z = cos_ele;
    Vector3 half_vector = to_world(frame, micronormal);
    Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
    return BSDFSampleRecord{
        reflected,
        1.5 /* eta */, Real(1) };
}

TextureSpectrum get_texture_op::operator()(const DisneyClearcoat& bsdf) const {
    return make_constant_spectrum_texture(make_zero_spectrum());
}