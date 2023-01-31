inline Real fss(Real f90, Vector3 omega, Vector3 n) {
    return (1 + (f90 - 1) * pow(1 - abs(dot(n, omega)), 5));
}

Spectrum eval_op::operator()(const DisneyDiffuse &bsdf) const {
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
    Real half_out = dot(h, dir_out);
    Real fss90 = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool) * (half_out * half_out);
    Real fd90 = 0.5 + 2 * fss90;
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Spectrum base_diffuse = base_color * c_INVPI * fss(fd90, dir_in, frame.n) * fss(fd90, dir_out, frame.n) * abs(dot(frame.n, dir_out));

    // handle subsurface scattering
    Spectrum subsurface = 1.25 * base_color * c_INVPI * (fss(fss90, dir_in, frame.n) * fss(fss90, dir_out, frame.n) * (1.0 / (abs(dot(frame.n, dir_in)) + abs(dot(frame.n, dir_out))) - 0.5) + 0.5) * abs(dot(frame.n, dir_out));
    Real ss = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);
    return (1.0 - ss) * base_diffuse + ss * subsurface;
}

Real pdf_sample_bsdf_op::operator()(const DisneyDiffuse &bsdf) const {
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
    // importance sample the cosine hemisphere domain.
    return fmax(dot(frame.n, dir_out), Real(0)) * c_INVPI;
}

std::optional<BSDFSampleRecord> sample_bsdf_op::operator()(const DisneyDiffuse &bsdf) const {
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
    return BSDFSampleRecord{
        to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
        Real(0),
        eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool)
    };
}

TextureSpectrum get_texture_op::operator()(const DisneyDiffuse &bsdf) const {
    return bsdf.base_color;
}
