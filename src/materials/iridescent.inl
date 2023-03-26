#include "../microfacet.h"

Spectrum eval_op::operator()(const Iridescent& bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
        dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        return make_zero_spectrum();
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Implement spectral bsdf
    // External parameters
    // Number thin film layer
    Real N = eval(bsdf.N, vertex.uv, vertex.uv_screen_size, texture_pool);
    // IOR of the thin film
    Real eta = eval(bsdf.eta, vertex.uv, vertex.uv_screen_size, texture_pool);
    // width of a cuticle layer
    Real a = eval(bsdf.a, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Number of cuticle layer
    Real M = eval(bsdf.M, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Interval of cuticle layers
    Real d = eval(bsdf.d, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real thick_film = eval(bsdf.thick_film, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real thick_air = eval(bsdf.thick_air, vertex.uv, vertex.uv_screen_size, texture_pool);

    // Internal parameters
    Real coeff = 0.1;
    Real m_order = 1;

    WaveSpectrum response(Real(0));
    WaveSpectrum diffuse(Real(0));
    WaveSpectrum lamellae(Real(0));

    Real cos_in = dot(dir_in, frame.n);
    Real cos_out = dot(dir_out, frame.n);
    Real u = sin(acos(cos_in)) + sin(acos(cos_out));
    Real v = cos_in + cos_out;
    Real sin_out = sin(acos(cos_out));

    for (int i = 0; i < nSpectralSamples; i++)
    {
        Real wavelength = WaveSpectrum::wavelength_beg + i * WaveSpectrum::interval + 0.5 * WaveSpectrum::interval;
        
        //// multiple thin film model///////////////////////
        Real sigma_b = (4 * c_PI / wavelength) * (eta * thick_film * cos(asin(sin(acos(dot(dir_in, frame.n))) / eta)) + thick_air * dot(dir_in, frame.n));
        if (cos(sigma_b) > 0) {
            diffuse[i] = coeff * pow(abs(cos(sigma_b)), m_order);
        }
        else {
            diffuse[i] = 0;
        }
        //////////////////////////////////////

        //// lamellae model//////////////////
        lamellae[i] = N * 0.5 * cos_in * cos_in;
        lamellae[i] *= pow(sin(c_PI * d * v * M / wavelength), 2) / pow(sin(c_PI * d * v / wavelength), 2);
        lamellae[i] *= pow(sin(c_PI * a * u / wavelength), 2) / pow(c_PI * a* u / wavelength, 2);
        /////////////////////////////////////        
    }
    response += (diffuse * cos_out) / (cos_out + sin_out);
    lamellae *= 0.01;
    response += (lamellae * sin_out) / (cos_out + sin_out);
    response.Clamp();
    Spectrum response_rgb = response.to_RGB();
    return response_rgb;
}

Real pdf_sample_bsdf_op::operator()(const Iridescent& bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
        dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        return Real(0);
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }

    // For Lambertian, we importance sample the cosine hemisphere domain.
    return fmax(dot(frame.n, dir_out), Real(0)) / c_PI;
}

std::optional<BSDFSampleRecord>
sample_bsdf_op::operator()(const Iridescent& bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0) {
        // No light below the surface
        return {};
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Implement
    std::optional<BSDFSampleRecord> record;
    return record;
}

TextureSpectrum get_texture_op::operator()(const Iridescent& bsdf) const {
    // TODO: remove base color
    return make_constant_spectrum_texture(make_zero_spectrum());
}
