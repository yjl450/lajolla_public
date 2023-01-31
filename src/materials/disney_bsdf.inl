#include "../microfacet.h"

#define DIFFUSE
#define METAL
#define CLEARCOAT
#define GLASS
#define SHEEN

Spectrum eval_op::operator()(const DisneyBSDF& bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
        dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    // Material parameters
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real subsurface = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular = eval(bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    Real specular_tint = eval(bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen = eval(bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen_tint = eval(bsdf.sheen_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real eta = bsdf.eta;

    // Derived parameters 
    constexpr Real alpha_min = 0.0001;
    Spectrum c_tint = luminance(base_color) > 0 ? base_color / luminance(base_color) : make_const_spectrum(1.0);
    Real aspect = sqrt(1.0 - 0.9 * anisotropic);
    Real alpha_x, alpha_y;
    alpha_x = max(alpha_min, roughness * roughness / aspect);
    alpha_y = max(alpha_min, roughness * roughness * aspect);
    Real lambda_in = lambda(dir_in, frame, alpha_x, alpha_y);
    Real lambda_out = lambda(dir_out, frame, alpha_x, alpha_y);
    Real g_in = 1.0 / (1.0 + lambda_in);
    Real g_out = 1.0 / (1.0 + lambda_out);

    // Precompute normal product
    Vector3 h = normalize(dir_in + dir_out);
    Real h_dot_in = abs(dot(h, dir_in));
    Real h_dot_out = abs(dot(h, dir_out));
    Real n_dot_in = abs(dot(frame.n, dir_in));
    Real n_dot_out = abs(dot(frame.n, dir_out));
    Spectrum f_disney = make_zero_spectrum();
    Vector3 hl = to_local(frame, h);

#ifdef DIFFUSE
    // Diffuse component
    if (dot(dir_in, vertex.geometric_normal) > 0 && dot(dir_in, vertex.geometric_normal) >= 0)
    {
        Real fss90 = roughness * (h_dot_out * h_dot_out);
        Real fd90 = 0.5 + 2 * fss90;
        Spectrum base_diffuse = base_color * c_INVPI * fss(fd90, dir_in, frame.n) * fss(fd90, dir_out, frame.n) * n_dot_out;
        Spectrum diffuse_subsurface = 1.25 * base_color * c_INVPI * (fss(fss90, dir_in, frame.n) * fss(fss90, dir_out, frame.n) * (1.0 / (n_dot_in + n_dot_out) - 0.5) + 0.5) * n_dot_out;
        Spectrum f_diffuse = (1.0 - subsurface) * base_diffuse + subsurface * diffuse_subsurface;
        f_disney += (1 - specular_transmission) * (1 - metallic) * f_diffuse;
    }
#endif // DIFFUSE

#ifdef SHEEN
    if (dot(dir_in, vertex.geometric_normal) > 0 && dot(dir_in, vertex.geometric_normal) >= 0)
    {
        Spectrum c_sheen = (1.0 - sheen_tint) + sheen_tint * c_tint;
        Spectrum f_sheen = c_sheen * pow(1.0 - h_dot_out, 5) * n_dot_out;
        f_disney += (1 - metallic) * sheen * f_sheen;
    }
#endif // SHEEN

#ifdef METAL
    if (dot(dir_in, vertex.geometric_normal) > 0 && dot(dir_in, vertex.geometric_normal) >= 0)
    {
        Spectrum ks = (1 - specular_tint) + specular_tint * c_tint;
        Real r0_m = pow(eta - 1.0, 2) / pow(eta + 1.0, 2);
        Spectrum c0 = specular * r0_m * (1 - metallic) * ks + metallic * base_color;
        Spectrum fm = c0 + (1.0 - c0) * pow(1.0 - h_dot_out, 5);
        Real dm = 1.0 / (c_PI * alpha_x * alpha_y * pow(pow(hl.x / alpha_x, 2) + pow(hl.y / alpha_y, 2) + hl.z * hl.z, 2));
        Real gm = g_in * g_out;
        Spectrum f_metal = fm * dm * gm / (4.0 * n_dot_in);
        f_disney += (1.0 - specular_transmission * (1 - metallic)) * f_metal;
    }
#endif // METAL

#ifdef CLEARCOAT
    if (dot(dir_in, vertex.geometric_normal) > 0 && dot(dir_in, vertex.geometric_normal) >= 0)
    {
        Real r0_c = pow(1.5 - 1.0, 2) / pow(1.5 + 1.0, 2);
        Real fc = r0_c + (1.0 - r0_c) * pow(1.0 - h_dot_out, 5);
        Real alpha_g2 = pow((1.0 - clearcoat_gloss) * 0.1 + clearcoat_gloss * 0.001, 2);
        Real dc = (alpha_g2 - 1.0) * c_INVPI / (log(alpha_g2) * (1.0 + (alpha_g2 - 1) * (hl.z * hl.z)));
        Real gc = (1.0 / (1.0 + lambda(dir_in, frame, 0.25, 0.25))) * (1.0 / (1.0 + lambda(dir_out, frame, 0.25, 0.25)));
        Real f_clearcoat = fc * dc * gc / (4.0 * n_dot_in);
        f_disney += 0.25 * clearcoat * make_const_spectrum(f_clearcoat);
    }
#endif // CLEARCOAT

#ifdef GLASS
    Spectrum f_glass;
    struct DisneyGlass glass = { bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta };
    f_glass = eval_op{ dir_in, dir_out,vertex,texture_pool,dir }(glass);
    f_disney += (1 - metallic) * specular_transmission * f_glass;
#endif // GLASS

    return f_disney;
}

Real pdf_sample_bsdf_op::operator()(const DisneyBSDF& bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
        dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real w_diffuse = (1 - metallic) * (1 - specular_transmission);
    Real w_metal = (1 - specular_transmission * (1 - metallic));
    Real w_clearcoat = 0.25 * clearcoat;
    Real w_glass = (1 - metallic) * specular_transmission;

    struct DisneyGlass glass { bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta };
    Real pdf_glass = pdf_sample_bsdf(glass, dir_in, dir_out, vertex, texture_pool);
    if (dot(vertex.geometric_normal, dir_in) <= 0) {
        return pdf_glass;
    }

    struct DisneyMetal metal { bsdf.base_color, bsdf.roughness, bsdf.anisotropic };
    struct DisneyDiffuse diffuse { bsdf.base_color, bsdf.roughness, bsdf.subsurface };
    struct DisneyClearcoat m_clearcoat { bsdf.clearcoat_gloss };

    Real pdf_metal = pdf_sample_bsdf(metal, dir_in, dir_out, vertex, texture_pool);
    Real pdf_diffuse = pdf_sample_bsdf(diffuse, dir_in, dir_out, vertex, texture_pool);
    Real pdf_clearcoat = pdf_sample_bsdf(m_clearcoat, dir_in, dir_out, vertex, texture_pool);
    Real pdf = w_metal * pdf_metal + w_diffuse * pdf_diffuse + w_clearcoat * pdf_clearcoat +w_glass * pdf_glass;
    return pdf / (w_metal + w_diffuse + w_clearcoat+w_glass);
}

std::optional<BSDFSampleRecord>
sample_bsdf_op::operator()(const DisneyBSDF& bsdf) const {
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    struct DisneyGlass glass { bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta };
    if (dot(vertex.geometric_normal, dir_in) <= 0) {
        return sample_bsdf(glass, dir_in, vertex, texture_pool, rnd_param_uv, rnd_param_w);
    }

    std::optional<BSDFSampleRecord> record;
    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);

    struct DisneyMetal metal { bsdf.base_color, bsdf.roughness, bsdf.anisotropic };
    struct DisneyDiffuse diffuse { bsdf.base_color, bsdf.roughness, bsdf.subsurface };
    struct DisneyClearcoat m_clearcoat { bsdf.clearcoat_gloss };
    Real w_diffuse = (1 - metallic) * (1 - specular_transmission);
    Real w_metal = (1 - specular_transmission * (1 - metallic));
    Real w_clearcoat = 0.25 * clearcoat;
    Real w_glass = (1 - metallic) * specular_transmission;
    Real w_sum = w_diffuse + w_metal + w_clearcoat + w_glass;
    Real choice = rnd_param_w * w_sum;
    if (choice <= w_diffuse) {
        record = sample_bsdf(diffuse, dir_in, vertex, texture_pool, rnd_param_uv, rnd_param_w);
    }
    else if (choice > w_diffuse && choice <= w_diffuse + w_metal) {
        record = sample_bsdf(metal, dir_in, vertex, texture_pool, rnd_param_uv, rnd_param_w);
    }
    else if (choice > w_diffuse + w_metal && choice <= w_diffuse + w_metal + w_clearcoat) {
        record = sample_bsdf(m_clearcoat, dir_in, vertex, texture_pool, rnd_param_uv, rnd_param_w);
    }
    else if (choice > w_diffuse + w_metal + w_clearcoat ) {
        record = sample_bsdf(glass, dir_in, vertex, texture_pool, rnd_param_uv, rnd_param_w);
    }
    return record;
}

TextureSpectrum get_texture_op::operator()(const DisneyBSDF& bsdf) const {
    return bsdf.base_color;
}
