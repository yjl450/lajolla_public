#pragma once
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

Spectrum L_s1(const Scene& scene, Vector3 dir_in, Vector3 p, PhaseFunction phase_function, Spectrum sigma_t, pcg32_state& rng) {
    int light_id = sample_light(scene, next_pcg32_real<Real>(rng));
    //int light_id = (int)(next_pcg32_real<Real>(rng) * 2);
    Light light = scene.lights[light_id];
    // Random parameters for light sampling
    Vector2 uv(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng));
    Real w = next_pcg32_real<Real>(rng);
    PointAndNormal point_normal = sample_point_on_light(light, p, uv, w, scene);
    // Evaluate phase function
    Vector3 dir_out = normalize(point_normal.position - p); // p(t) to p'(on light)
    Vector3 phase = eval(phase_function, -dir_in, dir_out);
    Spectrum Le = emission(light, -dir_out, Real(0), point_normal, scene);
    Spectrum transmittance = exp(-sigma_t * length(p - point_normal.position));
    // Test visibility TODO
    Real g = max(Real(0), dot(-dir_out, point_normal.normal)) / length_squared(p - point_normal.position);
    // impotance sampling lights + importance sampling point on light
    Real L_s1_pdf = pdf_point_on_light(light, point_normal, p, scene) * light_pmf(scene, light_id);
    Spectrum L_s1_estimate = phase * Le * transmittance * g;
    return L_s1_estimate / L_s1_pdf;
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
        Spectrum sigma_a = get_sigma_a(medium, Vector3(1, 1, 1));
        Spectrum sigma_s = get_sigma_s(medium, Vector3(1, 1, 1));
        Spectrum sigma_t = sigma_a + sigma_s;
        PhaseFunction phase_function = get_phase_function(medium);
        Spectrum t = -log(1 - next_pcg32_real<Real>(rng)) / sigma_t;
        Spectrum transmittance = exp(-sigma_t * t);
        Spectrum trans_pdf = transmittance * sigma_t;
        Vector3 p = ray.org + t * ray.dir;
        Spectrum scattering = L_s1(scene, ray.dir, p, phase_function, sigma_t, rng);
        //trans_pdf = exp(-sigma_t * t) + exp(-sigma_t * t) * sigma_t;
        return (transmittance / trans_pdf) * sigma_s * scattering;
    }
    PathVertex vertex = *vertex_;
    Medium medium = scene.media[vertex.exterior_medium_id];
    Spectrum sigma_a = get_sigma_a(medium, vertex.position);
    Spectrum sigma_s = get_sigma_s(medium, vertex.position);
    Spectrum sigma_t = sigma_a + sigma_s;
    PhaseFunction phase_function = get_phase_function(medium);
    Spectrum t = -log(1 - next_pcg32_real<Real>(rng)) / sigma_t;
    if (t.x < length(vertex.position - ray.org)) {
        //return make_zero_spectrum();
        Spectrum trans_pdf = exp(-sigma_t * t) * sigma_t;
        Spectrum transmittance = exp(-sigma_t * t);
        Vector3 p = ray.org + t * ray.dir;
        Spectrum scattering = L_s1(scene, ray.dir, p, phase_function, sigma_t, rng);
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

int update_medium(Ray ray, PathVertex vertex) {
    int medium_id;
    if (vertex.interior_medium_id != vertex.exterior_medium_id) {
        // If ray exit the object, return exterior id
        if (dot(ray.dir, vertex.geometric_normal) > 0) {
            medium_id = vertex.exterior_medium_id;
        }
        else {
            medium_id = vertex.interior_medium_id;
        }
    }
    return medium_id;
}


// The third volumetric renderer (not so simple anymore): 
// multiple monochromatic homogeneous volumes with multiple scattering
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_3(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
        (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{ Real(0), Real(0) };
    int current_medium_id = scene.camera.medium_id;
    int bounce = 0;
    Spectrum radiance = make_zero_spectrum();
    Spectrum current_path_throughput = make_const_spectrum(1.f);
    while (1) {
        Medium medium = scene.media[current_medium_id];
        bool scatter = false;
        std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
        Spectrum transmittance = make_const_spectrum(1.f);
        Spectrum trans_pdf = make_const_spectrum(1.f);
        Spectrum t = make_zero_spectrum();
        if (current_medium_id != -1) {
            Spectrum sigma_a, sigma_s, sigma_t;
            Real t_hit = -1;
            if (vertex_) {
                PathVertex vertex = *vertex_;
                t_hit = length(vertex.position - ray.org);
                sigma_a = get_sigma_a(medium, vertex.position);
                sigma_s = get_sigma_s(medium, vertex.position);
            }
            else {
                sigma_a = get_sigma_a(medium, Vector3(1, 1, 1));
                sigma_s = get_sigma_s(medium, Vector3(1, 1, 1));
            }
            sigma_t = sigma_a + sigma_s;
            t = -log(1 - next_pcg32_real<Real>(rng)) / sigma_t;
            if ((!vertex_) || (vertex_ && t.x < t_hit))
            {
                scatter = true;
                transmittance = exp(-sigma_t * t);
                trans_pdf = transmittance * sigma_t;
            }
            else {
                transmittance = exp(-sigma_t * t_hit);
                trans_pdf = transmittance;
            }
        }
        current_path_throughput *= (transmittance / trans_pdf);
        // Hit a surface, include the emission
        if (!scatter && vertex_) {
            PathVertex vertex = *vertex_;
            if (is_light(scene.shapes[vertex.shape_id])) {
                Spectrum Le = emission(vertex, -ray.dir, scene);
                radiance += current_path_throughput * Le;
                break;
            }
        }
        // Reach maximum bounce
        if (scene.options.max_depth != -1 && bounce == scene.options.max_depth - 1) {
            break;
        }
        // Hit a surface, update medium
        if (!scatter && vertex_) {
            PathVertex vertex = *vertex_;
            if (vertex.material_id == -1) {
                current_medium_id = update_medium(ray, vertex);
                bounce++;
                ray.org = vertex.position;
                continue;
            }
        }
        //sample next direct & update path throughput
        if (scatter) {
            Vector3 p = ray.org + t * ray.dir;
            PhaseFunction phase_function = get_phase_function(medium);
            Spectrum sigma_s = get_sigma_s(medium, p);
            Vector2 uv(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng));
            std::optional<Vector3> next_dir_ = sample_phase_function(phase_function, -ray.dir, uv);
            if (next_dir_) {
                Vector3 next_dir = *next_dir_;
                Spectrum phase = eval(phase_function, -ray.dir, next_dir);
                current_path_throughput *=
                    (phase / pdf_sample_phase(phase_function, -ray.dir, next_dir)) * sigma_s;
                ray = Ray{ p,next_dir, Real(0), infinity<Real>() };
            }
        }
        else {
            break;
        }
        // Russian roulette
        if (bounce >= scene.options.rr_depth) {
            Real rr_prob = min(current_path_throughput.x, 0.95);
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                break;
            }
            else {
                current_path_throughput /= rr_prob;
            }
        }
        bounce++;
    }
    return radiance;
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
