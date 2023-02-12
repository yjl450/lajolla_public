#pragma once

Spectrum L_s1(const Scene& scene, Vector3 dir_in, Vector3 p, Medium medium, Spectrum sigma_t, pcg32_state& rng) {
    int light_id = sample_light(scene, next_pcg32_real<Real>(rng));
    //int light_id = (int)(next_pcg32_real<Real>(rng) * 2);
    Light light = scene.lights[light_id];
    // Random parameters for light sampling
    Vector2 uv(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng));
    Real w = next_pcg32_real<Real>(rng);
    PointAndNormal point_normal = sample_point_on_light(light, p, uv, w, scene);
    // Evaluate phase function
    Vector3 dir_out = normalize(point_normal.position - p); // p(t) to p'(on light)
    PhaseFunction phase_function = get_phase_function(medium);
    Vector3 phase = eval(phase_function, -dir_in, dir_out);
    Spectrum Le = emission(light, -dir_out, Real(0), point_normal, scene);
    Spectrum transmittance = exp(-sigma_t * length(p - point_normal.position));
    // Test visibility TODO
    Real g = max(Real(0), dot(-dir_out, point_normal.normal)) / pow(length(p - point_normal.position), 2);
    // impotance sampling lights + importance sampling point on light
    Real L_s1_pdf = pdf_point_on_light(light, point_normal, p, scene) * light_pmf(scene, light_id);
    Spectrum L_s1_estimate = phase * Le * transmittance * g;
    return L_s1_estimate / L_s1_pdf;
}

// The simplest volumetric renderer: 
// single absorption only homogeneous volume
// only handle directly visible light sources
Spectrum vol_path_tracing_1(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
        (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos); 
    RayDifferential ray_diff = RayDifferential{ Real(0), Real(0) };
    std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
    if (!vertex_) {
        // Hit background. Account for the environment map if needed.
        if (has_envmap(scene)) {
            const Light& envmap = get_envmap(scene);
            return emission(envmap,
                -ray.dir, // pointing outwards from light
                ray_diff.spread,
                PointAndNormal{}, // dummy parameter for envmap
                scene);
        }
        return make_zero_spectrum();
    }
    PathVertex vertex = *vertex_;
    Medium scene_medium = scene.media[vertex.exterior_medium_id];
    Spectrum sigma_a = get_sigma_a(scene_medium, vertex.position);
    Spectrum transmittance = exp(-sigma_a * distance(vertex.position, ray.org));
    Spectrum Le = make_zero_spectrum();
    if (is_light(scene.shapes[vertex.shape_id])) {
        Le = emission(vertex, -ray.dir, scene);
    }
    return transmittance * Le;
}

// The second simplest volumetric renderer: 
// single monochromatic homogeneous volume with single scattering,
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_2(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
        (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{ Real(0), Real(0) };
    std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
    if (!vertex_) {
        //return make_zero_spectrum();
        Medium medium = scene.media[scene.camera.medium_id];
        Spectrum sigma_a = get_sigma_a(medium, Vector3(2,2,2));
        Spectrum sigma_s = get_sigma_s(medium, Vector3(2,2,2));
        Spectrum sigma_t = sigma_a + sigma_s;
        Spectrum t = -log(1 - next_pcg32_real<Real>(rng)) / sigma_t;
        Spectrum transmittance = exp(-sigma_t * t);
        Spectrum trans_pdf = exp(-sigma_t * t) * sigma_t;
        Vector3 p = ray.org + t * ray.dir;
        Spectrum scattering = L_s1(scene, ray.dir, p, medium, sigma_t, rng);
        //trans_pdf = exp(-sigma_t * t) + exp(-sigma_t * t) * sigma_t;
        return (transmittance / trans_pdf) * sigma_s * scattering;
    }
    PathVertex vertex = *vertex_;
    Medium medium = scene.media[vertex.exterior_medium_id];
    Spectrum sigma_a = get_sigma_a(medium, vertex.position);
    Spectrum sigma_s = get_sigma_s(medium, vertex.position);
    Spectrum sigma_t = sigma_a + sigma_s;
    Spectrum t = -log(1 - next_pcg32_real<Real>(rng)) / sigma_t;
    if (t.x < length(vertex.position - ray.org)) {
        //return make_zero_spectrum();
        Spectrum trans_pdf = exp(-sigma_t * t) * sigma_t;
        Spectrum transmittance = exp(-sigma_t * t);
        Vector3 p = ray.org + t * ray.dir;
        Spectrum scattering = L_s1(scene, ray.dir, p, medium, sigma_t, rng);
        return (transmittance / trans_pdf) * sigma_s * scattering;
    }
    else {
        Spectrum trans_pdf = exp(-sigma_t * length(vertex.position - ray.org));
        Spectrum transmittance = exp(-sigma_t * length(vertex.position - ray.org));
        Spectrum Le = make_zero_spectrum();
        if (is_light(scene.shapes[vertex.shape_id])) {
            Le = emission(vertex, -ray.dir, scene);
        }
        return (transmittance / trans_pdf) * Le;
    }
}

// The third volumetric renderer (not so simple anymore): 
// multiple monochromatic homogeneous volumes with multiple scattering
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_3(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    return make_zero_spectrum();
}

// The fourth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// still no surface lighting
Spectrum vol_path_tracing_4(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    return make_zero_spectrum();
}

// The fifth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing_5(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    return make_zero_spectrum();
}

// The final volumetric renderer: 
// multiple chromatic heterogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing(const Scene &scene,
                          int x, int y, /* pixel coordinates */
                          pcg32_state &rng) {
    // Homework 2: implememt this!
    return make_zero_spectrum();
}
