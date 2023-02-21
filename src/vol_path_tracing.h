#pragma once

inline int clamp_channel(Real u) {
    return min(max((int)(u * 3), 0), 2);
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
    Real g = abs(dot(-dir_out, point_normal.normal)) / length_squared(p - point_normal.position);
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
    if (vertex.interior_medium_id != vertex.exterior_medium_id) {
        // If ray exit the object, return exterior id
        if (dot(ray.dir, vertex.geometric_normal) > 0) {
            return vertex.exterior_medium_id;
        }
        else {
            return vertex.interior_medium_id;
        }
    }
    return vertex.exterior_medium_id;
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
                ray.org = vertex.position + get_intersection_epsilon(scene) * ray.dir;
                continue;
            }
        }
        //sample next direction & update path throughput
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
                ray = Ray{ p + get_intersection_epsilon(scene) * next_dir,next_dir, Real(0), infinity<Real>() };
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

Spectrum next_event(const Scene& scene, Vector3 p, Vector3 dir_in, int medium_id, int bounce, pcg32_state& rng, bool eval_mat = false, Material mat = Material(), PathVertex v = PathVertex()) {
    // Sample light
    int light_id = sample_light(scene, next_pcg32_real<Real>(rng));
    Light light = scene.lights[light_id];
    Vector2 uv = Vector2(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng));
    PointAndNormal pt_normal = sample_point_on_light(light, p, uv, next_pcg32_real<Real>(rng), scene);
    Real pdf_NEE = pdf_point_on_light(light, pt_normal, p, scene) * light_pmf(scene, light_id);
    Vector3 p_prime = pt_normal.position;
    // Init var
    int shadow_medium_id = medium_id;
    Spectrum transmittance_light = make_const_spectrum(1.f);
    Spectrum pdf_trans_dir = make_const_spectrum(1.f);
    int shadow_bounce = 0;
    RayDifferential ray_diff = RayDifferential{ Real(0), Real(0) };
    Vector3 init_p = p;
    while (1) {
        Ray shadow_ray = Ray{ p, normalize(p_prime - p), get_shadow_epsilon(scene), (1 - get_shadow_epsilon(scene)) * length(p_prime - p) };
        std::optional<PathVertex> vertex_ = intersect(scene, shadow_ray, ray_diff);
        Real next_t = length(p_prime - p);
        if (vertex_) {
            PathVertex vertex = *vertex_;
            next_t = length(vertex.position - p);
        }
        if (shadow_medium_id != -1) {
            Medium medium = scene.media[shadow_medium_id];
            Vector3 position = shadow_ray.org + next_t * shadow_ray.dir;
            Spectrum sigma_a = get_sigma_a(medium, position);
            Spectrum sigma_s = get_sigma_s(medium, position);
            Spectrum sigma_t = sigma_s + sigma_a;
            transmittance_light *= exp(-sigma_t * next_t);
            pdf_trans_dir *= exp(-sigma_t * next_t);
        }
        if (!vertex_) {
            break;
        }
        else {
            PathVertex vertex = *vertex_;
            if (vertex.material_id >= 0) {
                return make_zero_spectrum(); // blocked
            }
            shadow_bounce++;
            if (scene.options.max_depth != -1 && bounce + shadow_bounce + 1 >= scene.options.max_depth) {
                return make_zero_spectrum();
            }
            shadow_medium_id = update_medium(shadow_ray, vertex);
            p += next_t * shadow_ray.dir;
        }
    }
    if (avg(transmittance_light) > 0) {
        Vector3 dir_out = normalize(init_p - p_prime);
        Spectrum f;
        Real pdf_dir;
        Spectrum Le = emission(light, dir_out, Real(0), pt_normal, scene);
        Real g = abs(dot(dir_out, pt_normal.normal)) / length_squared(init_p - p_prime);
        if (eval_mat) {
            f = eval(mat, -dir_in, -dir_out, v, scene.texture_pool);
            pdf_dir = pdf_sample_bsdf(mat, -dir_in, -dir_out, v, scene.texture_pool) * g * pdf_trans_dir.x;
        } else {
            PhaseFunction phase_function = get_phase_function(scene.media[medium_id]);
            f = eval(phase_function, -dir_in, -dir_out);
            pdf_dir = pdf_sample_phase(phase_function, -dir_in, -dir_out) * g * pdf_trans_dir.x;
        }
        Spectrum contrib = transmittance_light * g * f * Le / pdf_NEE;
        Real w = pow(pdf_NEE, 2) / (pow(pdf_dir, 2) + pow(pdf_NEE, 2));
        return w * contrib;
    }
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
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
        (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{ Real(0), Real(0) };
    int current_medium_id = scene.camera.medium_id;
    int bounce = 0;
    Spectrum radiance = make_zero_spectrum();
    Spectrum current_path_throughput = make_const_spectrum(1.f);
    Real dir_pdf = 1;
    Real multi_trans_pdf = 1;
    Vector3 p; // cached p
    bool never_scatter = true;
    while (1) {
        bool scatter = false;
        std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
        Spectrum transmittance = make_const_spectrum(1.f);
        Spectrum trans_pdf = make_const_spectrum(1.f);
        Spectrum t = make_zero_spectrum();
        Spectrum phase = make_const_spectrum(1);
        if (current_medium_id != -1) {
            Medium medium = scene.media[current_medium_id];
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
                never_scatter = false;
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
                if (never_scatter) {
                    radiance += current_path_throughput * Le;
                    break;
                }
                else {
                    // account for NEE
                    int light_id = get_area_light_id(scene.shapes[vertex.shape_id]);
                    Light light = scene.lights[light_id];
                    PointAndNormal pt_normal{ vertex.position, vertex.geometric_normal };
                    Real pdf_NEE = pdf_point_on_light(light, pt_normal, p, scene) * light_pmf(scene, light_id);
                    Real g = abs(dot(-ray.dir, vertex.geometric_normal)) / length_squared(p - vertex.position);
                    Real pdf_phase = dir_pdf * multi_trans_pdf * g;
                    Real w = pow(pdf_phase, 2) / (pow(pdf_phase, 2) + pow(pdf_NEE, 2));
                    radiance += current_path_throughput * Le * w;
                    break;
                }
            }
        }
        // Reach maximum bounce
        if (scene.options.max_depth != -1 && bounce == scene.options.max_depth - 1) {
            break;
        }
        // Hit a surface, update medium
        if (!scatter && vertex_) {
            PathVertex vertex = *vertex_;
            multi_trans_pdf *= trans_pdf.x;
            if (vertex.material_id == -1) {
                current_medium_id = update_medium(ray, vertex);
                bounce++;
                ray.org = vertex.position + get_intersection_epsilon(scene) * ray.dir;
                continue;
            }
        }
        //sample next direction & update path throughput
        if (scatter) {
            Medium medium = scene.media[current_medium_id];
            p = ray.org + (t + get_intersection_epsilon(scene)) * ray.dir;
            Spectrum sigma_s = get_sigma_s(medium, p);

            // Next event estimation
            Spectrum NEE = next_event(scene, p, ray.dir, current_medium_id, bounce, rng);
            radiance += current_path_throughput * sigma_s * NEE;

            // Next bounce
            PhaseFunction phase_function = get_phase_function(medium);
            Vector2 uv(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng));
            std::optional<Vector3> next_dir_ = sample_phase_function(phase_function, -ray.dir, uv);
            if (next_dir_) {
                multi_trans_pdf = 1;
                Vector3 next_dir = *next_dir_;
                phase = eval(phase_function, -ray.dir, next_dir);
                dir_pdf = pdf_sample_phase(phase_function, ray.dir, next_dir);
                current_path_throughput *=
                    (phase / dir_pdf) * sigma_s;
                ray = Ray{ p ,next_dir, Real(0), infinity<Real>() };
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

// The fifth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing_5(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
        (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = init_ray_differential(w, h);
    int current_medium_id = scene.camera.medium_id;
    int bounce = 0;
    Spectrum radiance = make_zero_spectrum();
    Spectrum current_path_throughput = make_const_spectrum(1.f);
    Real dir_pdf = 1;
    Real multi_trans_pdf = 1;
    Vector3 p; // cached p
    bool never_scatter = true;
    Real eta_scale = Real(1);
    while (1) {
        bool scatter = false;
        std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
        Spectrum transmittance = make_const_spectrum(1.f);
        Spectrum trans_pdf = make_const_spectrum(1.f);
        Spectrum t = make_zero_spectrum();
        Spectrum phase = make_const_spectrum(1);
        if (current_medium_id != -1) {
            Medium medium = scene.media[current_medium_id];
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
                never_scatter = false;
                transmittance = exp(-sigma_t * t);
                trans_pdf = transmittance * sigma_t;
            }
            else {
                transmittance = exp(-sigma_t * t_hit);
                trans_pdf = transmittance;
            }
        }
        else if (!vertex_) {
            break;
        }
        current_path_throughput *= (transmittance / trans_pdf);
        // Hit a surface, include the emission
        if (!scatter && vertex_) {
            PathVertex vertex = *vertex_;
            if (is_light(scene.shapes[vertex.shape_id])) {
                Spectrum Le = emission(vertex, -ray.dir, scene);
                if (never_scatter) {
                    radiance += current_path_throughput * Le;
                }
                else {
                    // account for NEE
                    int light_id = get_area_light_id(scene.shapes[vertex.shape_id]);
                    Light light = scene.lights[light_id];
                    PointAndNormal pt_normal{ vertex.position, vertex.geometric_normal };
                    Real pdf_NEE = pdf_point_on_light(light, pt_normal, p, scene) * light_pmf(scene, light_id);
                    Real g = abs(dot(-ray.dir, vertex.geometric_normal)) / length_squared(p - vertex.position);
                    Real pdf_phase = dir_pdf * multi_trans_pdf * g;
                    Real w = pow(pdf_phase, 2) / (pow(pdf_phase, 2) + pow(pdf_NEE, 2));
                    radiance += current_path_throughput * Le * w;
                }
            }
        }
        // Reach maximum bounce
        if (scene.options.max_depth != -1 && bounce == scene.options.max_depth - 1) {
            break;
        }
        // Hit a surface, update medium
        if (!scatter && vertex_) {
            PathVertex vertex = *vertex_;
            multi_trans_pdf *= trans_pdf.x;
            if (vertex.material_id == -1) {
                current_medium_id = update_medium(ray, vertex);
                bounce++;
                ray.org = vertex.position + get_intersection_epsilon(scene) * ray.dir;
                continue;
            }
        }
        //sample next direction & update path throughput
        if (scatter) {
            Medium medium = scene.media[current_medium_id];
            p = ray.org + (t + get_intersection_epsilon(scene)) * ray.dir;
            Spectrum sigma_s = get_sigma_s(medium, p);

            // Next event estimation
            Spectrum NEE = next_event(scene, p, ray.dir, current_medium_id, bounce, rng);
            radiance += current_path_throughput * sigma_s * NEE;

            // Next bounce
            PhaseFunction phase_function = get_phase_function(medium);
            Vector2 uv(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng));
            std::optional<Vector3> next_dir_ = sample_phase_function(phase_function, -ray.dir, uv);
            if (next_dir_) {
                multi_trans_pdf = 1;
                Vector3 next_dir = *next_dir_;
                phase = eval(phase_function, -ray.dir, next_dir);
                dir_pdf = pdf_sample_phase(phase_function, ray.dir, next_dir);
                current_path_throughput *=
                    (phase / dir_pdf) * sigma_s;
                ray = Ray{ p ,next_dir, Real(0), infinity<Real>() };
            }
        }
        else {
            if (vertex_ && (*vertex_).material_id != -1) {
                PathVertex vertex = *vertex_;
                p = vertex.position;
                const Material mat = scene.materials[vertex.material_id];
                // Importance Sampling
                Spectrum NEE = next_event(scene, vertex.position, ray.dir, current_medium_id, bounce, rng, true, mat, vertex);
                radiance += current_path_throughput * NEE;
                never_scatter = false;

                //// BSDF Sampling direction
                Vector3 dir_view = -ray.dir;
                Vector2 bsdf_rnd_param_uv{ next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng) };
                Real bsdf_rnd_param_w = next_pcg32_real<Real>(rng);
                std::optional<BSDFSampleRecord> bsdf_sample_ =
                    sample_bsdf(mat,
                        dir_view,
                        vertex,
                        scene.texture_pool,
                        bsdf_rnd_param_uv,
                        bsdf_rnd_param_w);
                if (bsdf_sample_) {
                    multi_trans_pdf = 1;
                    const BSDFSampleRecord bsdf_sample = *bsdf_sample_;
                    Vector3 dir_bsdf = bsdf_sample.dir_out;
                    if (bsdf_sample.eta == 0) {
                        ray_diff.spread = reflect(ray_diff, vertex.mean_curvature, bsdf_sample.roughness);
                    }
                    else {
                        ray_diff.spread = refract(ray_diff, vertex.mean_curvature, bsdf_sample.eta, bsdf_sample.roughness);
                        eta_scale /= (bsdf_sample.eta * bsdf_sample.eta);
                    }
                    ray = Ray{ vertex.position, dir_bsdf, get_intersection_epsilon(scene), infinity<Real>() };
                    current_medium_id = update_medium(ray, vertex);
                    Spectrum f = eval(mat, dir_view, dir_bsdf, vertex, scene.texture_pool);
                    dir_pdf = pdf_sample_bsdf(mat, dir_view, dir_bsdf, vertex, scene.texture_pool);
                    current_path_throughput *= f / dir_pdf;
                }
                else {
                    break;
                }
            }
            else {
                break;
            }
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

// The final volumetric renderer: 
// multiple chromatic heterogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing(const Scene &scene,
    int x, int y, /* pixel coordinates */
                          pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
        (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = init_ray_differential(w, h);
    int current_medium_id = scene.camera.medium_id;
    int bounce = 0;
    Spectrum radiance = make_zero_spectrum();
    Spectrum current_path_throughput = make_const_spectrum(1.f);
    Real dir_pdf = 1;
    Spectrum multi_trans_pdf = make_const_spectrum(1);
    Vector3 p; // cached p
    bool never_scatter = true;
    Real eta_scale = Real(1);
    while (1) {
        bool scatter = false;
        std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
        Spectrum transmittance = make_const_spectrum(1.f);
        Spectrum phase = make_const_spectrum(1);
        Spectrum trans_dir_pdf = make_const_spectrum(1); // PDF for free-flight sampling
        Spectrum trans_nee_pdf = make_const_spectrum(1); // PDF for next event estimation
        Real accum_t = 0;
        if (current_medium_id != -1) { // In medium
            Real u = next_pcg32_real<Real>(rng);
            int channel = clamp_channel(u);
            Medium medium = scene.media[current_medium_id];
            Spectrum majorant = get_majorant(medium, ray);
            Real t_hit = infinity<Real>();
            int iteration = 0;
            while (1) {
                if (majorant[channel] <= 0) {
                    break;
                }
                if (iteration >= scene.options.max_null_collisions) {
                    break;
                }
                Spectrum sigma_a, sigma_s, sigma_t;
                if (vertex_) {
                    PathVertex vertex = *vertex_;
                    t_hit = distance(vertex.position, ray.org);
                }
                Real t = -log(1 - next_pcg32_real<Real>(rng)) / majorant[channel];
                Real dt = t_hit - accum_t;
                accum_t = min(accum_t + t, t_hit);
                if (t < dt) {
                    Vector3 current_position = ray.org + accum_t * ray.dir;
                    Spectrum sigma_a = get_sigma_a(medium, current_position);
                    Spectrum sigma_s = get_sigma_s(medium, current_position);
                    Spectrum sigma_t = sigma_a + sigma_s;
                    Spectrum real_prob = sigma_t / majorant;
                    if (next_pcg32_real<Real>(rng) < real_prob[channel]) {
                        scatter = true;
                        never_scatter = false;
                        transmittance *= exp(-majorant * t) / max(majorant);
                        trans_dir_pdf *= exp(-majorant * t) * majorant * real_prob / max(majorant);
                        break;
                    }
                    else {
                        Spectrum sigma_n = majorant * (1 - sigma_t / majorant);
                        transmittance *= exp(-majorant * t) * sigma_n / max(majorant);
                        trans_dir_pdf *= exp(-majorant * t) * majorant * (1 - real_prob) / max(majorant);
                        trans_nee_pdf *= exp(-majorant * t) * majorant / max(majorant);
                    }
                }
                else {
                    transmittance *= exp(-majorant * dt);
                    trans_dir_pdf *= exp(-majorant * dt);
                    trans_nee_pdf *= exp(-majorant * dt);
                    break;
                }
                iteration++;
            }
        }
        else if (!vertex_) {
            break;
        }
        current_path_throughput *= (transmittance / avg(trans_dir_pdf));
        // Hit a surface, include the emission
        if (!scatter && vertex_) {
            PathVertex vertex = *vertex_;
            if (is_light(scene.shapes[vertex.shape_id])) {
                Spectrum Le = emission(vertex, -ray.dir, scene);
                if (never_scatter) {
                    radiance += current_path_throughput * Le;
                }
                else {
                    // account for NEE
                    int light_id = get_area_light_id(scene.shapes[vertex.shape_id]);
                    Light light = scene.lights[light_id];
                    PointAndNormal pt_normal{ vertex.position, vertex.geometric_normal };
                    Real pdf_NEE = pdf_point_on_light(light, pt_normal, p, scene) * light_pmf(scene, light_id);
                    Real g = abs(dot(-ray.dir, vertex.geometric_normal)) / length_squared(p - vertex.position);
                    Real pdf_phase = dir_pdf * multi_trans_pdf.x * g;
                    Real w = pow(pdf_phase, 2) / (pow(pdf_phase, 2) + pow(pdf_NEE, 2));
                    radiance += current_path_throughput * Le * w;
                }
            }
        }
        // Reach maximum bounce
        if (scene.options.max_depth != -1 && bounce == scene.options.max_depth - 1) {
            break;
        }
        // Hit a surface, update medium
        if (!scatter && vertex_) {
            PathVertex vertex = *vertex_;
            multi_trans_pdf *= trans_dir_pdf;
            if (vertex.material_id == -1) {
                current_medium_id = update_medium(ray, vertex);
                bounce++;
                ray.org = vertex.position + get_intersection_epsilon(scene) * ray.dir;
                continue;
            }
        }
        //sample next direction & update path throughput
        if (scatter) {
            Medium medium = scene.media[current_medium_id];
            p = ray.org + (accum_t + get_intersection_epsilon(scene)) * ray.dir;
            Spectrum sigma_s = get_sigma_s(medium, p);

            // Next event estimation
            Spectrum NEE = next_event(scene, p, ray.dir, current_medium_id, bounce, rng);
            radiance += current_path_throughput * sigma_s * NEE;

            // Next bounce
            PhaseFunction phase_function = get_phase_function(medium);
            Vector2 uv(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng));
            std::optional<Vector3> next_dir_ = sample_phase_function(phase_function, -ray.dir, uv);
            if (next_dir_) {
                multi_trans_pdf = make_const_spectrum(1);
                Vector3 next_dir = *next_dir_;
                phase = eval(phase_function, -ray.dir, next_dir);
                dir_pdf = pdf_sample_phase(phase_function, ray.dir, next_dir);
                current_path_throughput *=
                    (phase / dir_pdf) * sigma_s;
                ray = Ray{ p ,next_dir, Real(0), infinity<Real>() };
            }
        }
        else {
            if (vertex_ && (*vertex_).material_id != -1) {
                PathVertex vertex = *vertex_;
                p = vertex.position;
                const Material mat = scene.materials[vertex.material_id];
                // Importance Sampling
                Spectrum NEE = next_event(scene, vertex.position, ray.dir, current_medium_id, bounce, rng, true, mat, vertex);
                radiance += current_path_throughput * NEE;
                never_scatter = false;

                //// BSDF Sampling direction
                Vector3 dir_view = -ray.dir;
                Vector2 bsdf_rnd_param_uv{ next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng) };
                Real bsdf_rnd_param_w = next_pcg32_real<Real>(rng);
                std::optional<BSDFSampleRecord> bsdf_sample_ =
                    sample_bsdf(mat,
                        dir_view,
                        vertex,
                        scene.texture_pool,
                        bsdf_rnd_param_uv,
                        bsdf_rnd_param_w);
                if (bsdf_sample_) {
                    multi_trans_pdf = make_const_spectrum(1);
                    const BSDFSampleRecord bsdf_sample = *bsdf_sample_;
                    Vector3 dir_bsdf = bsdf_sample.dir_out;
                    if (bsdf_sample.eta == 0) {
                        ray_diff.spread = reflect(ray_diff, vertex.mean_curvature, bsdf_sample.roughness);
                    }
                    else {
                        ray_diff.spread = refract(ray_diff, vertex.mean_curvature, bsdf_sample.eta, bsdf_sample.roughness);
                        eta_scale /= (bsdf_sample.eta * bsdf_sample.eta);
                    }
                    ray = Ray{ vertex.position, dir_bsdf, get_intersection_epsilon(scene), infinity<Real>() };
                    current_medium_id = update_medium(ray, vertex);
                    Spectrum f = eval(mat, dir_view, dir_bsdf, vertex, scene.texture_pool);
                    dir_pdf = pdf_sample_bsdf(mat, dir_view, dir_bsdf, vertex, scene.texture_pool);
                    current_path_throughput *= f / dir_pdf;
                }
                else {
                    break;
                }
            }
            else {
                break;
            }
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
