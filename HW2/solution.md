## 2
1. Change the absorption parameters to zero in scenes/volpath_test/volpath_test1.xml. What do you see? Why?
   
   The light source is shown without any absorption. This is because if $\sigma_a$ is set to zero, $\exp(\sigma_a*t)$ is always zero, meaning there is no absorption at all.

2. In the homework, we assume the volume being not emissive. If you were tasked to modify the pseudocode above to add volume emission, how would you do it? Briefly describe your approach.
   
   Since this is a homogeneous volume, we can assume the volume emission is the integral is constant $Le$ at each $dt$. Thus, if the primary ray hit an object, in addition to include the surface emission, we randomly sample a point $t\in(0, t_{hit})$ and include the emission at that point $e^{-\sigma_a*t}Le$. If the ray doesn't hit any object, we sample $t\in(0,\infty)$.

   In terms of code:

   `t_sample = rng()*t__hit`

   `return transmittance * Le_surface + exp(-sigma_a * t_sample) * Le_volume`

## 3
1. In the derivation above, how did we get from p(t) ∝ exp(−σtt) to p(t) = σt exp(σtt)?

    Take the integral of the pdf on the domain $t\in(0,\infty)$.
    $$
    \int_0^\infty \exp(-\sigma_t s)ds 
    = \lim_{s\rightarrow \infty}-\frac{1}{\sigma_t}\exp(-\sigma_t s) + \frac{1}{\sigma_t} 
    = 0 + \frac{1}{\sigma_t}
    $$
    To make this integral equal to 1, the pdf need to be scaled with $\sigma_t$. Thus, the pdf is $p(t)=\sigma_t\exp(-\sigma_t t)$.

2. How was Equation (13) incorporated into the pseudo code above? Why is it done this way?

    Equation 13 is used to compute the trans_pdf so that when the emission from a surface is included, the transmittance is devided by the trans_pdf. This is because the case includes all the possibilities of $t\geq t_{hit}$, so the overall probability is the integral on $t\in(t_{hit},\infty)$. Since the transmittance and pdf are the same, they cancel out each other and the emission wii be returned directly.

3. Play with the parameters σs and σa, how do they affect the final image? Why? (we assume monochromatic volumes, so don’t set the parameters to be different for each channel.)

    Both parameters will increase the energy loss as they get higher. $\sigma_a$ direct affects the absorption so the image gets darker significantly as $\sigma_a$ increases.

    $\sigma_s$ both contributes the energy loss and increase. So, a higher $\sigma_s$ does not dim the image as much. However, as the scattering increases, the image becomes more blurry as the light bounces more in the medium.

4. Change the phase function from isotropic (the default one) to Henyey-Greenstein by uncommenting
the phase function specification in the scene scenes/volpath_test/volpath_test2.xml. Play with the
parameter g (the valid range is (−1, 1)). What does g mean? How does the g parameter change the
appearance? Why?

    g means the intensity of forward scattering (positive) and back scattering (negative). As the absolute value of g gets bigger, the scattering intensity along the incident axis is stonger and scattering from the sides becomes weaker, which means the purple light becomes less visible. Back scattering also leads to less light entering the camera, making the image darker overall.

## 4
1. Play with the parameters σs, σa of different volumes, and change max_depth, how do they affect the final image? How does increasing/decreasing σs and σa of medium1 and medium2 affect the appearance, respectively? Why? Do different σs and σa values affect how high you should set max_depth?

    1) The higher the $\sigma_s$, the brighter part will become dimmer and the dark part will become brigher, which is kind of a blend effect. The higher the $\sigma_a$, every part of the medium simply becomes darker because more light is absorbed.
    2) medium1 controls the fog-like substance filling the entire space. medium2 controls the semi-transparent ball. They behave following the rule described above.
    3) A high $\sigma_s$ needs bigger max_depth because more iterations are needed to get accurate scattering effect. A high $\sigma_a$ is not affected by max_depth as much.

2. Switch to the Henyey-Greenstein phase function again. How does changing the g parameter affect the appearance? Why?
   
    A positive g makes the light scatters more towards the positive axis. So the light source appears brighter and has a bigger halo around it.

3. Propose a phase function yourself (don’t have to describe the exact mathematical form). How would you design the shape of the phase function? What parameter would you set to control it?

    Suppose there is a phase function with a main axis. Lights that are parallel to the main axis can mostly pass through and light that are perpendicular to the main axis will mostly be blocked. Two parameters: a direction vector of the main axis, and a g similar to the Henyey-Greenstein phase function.

## 5
1. When will next event estimation be more efficient than phase function sampling? In our test scenes, which one is more efficient? Why?

    When the phase function sampling struggles to arrive at a light source. Next event estimation is more efficient since it alone can produce almost noise-free results.

2. In scenes/volpath_test/volpath_test4_2.xml, we render a scene with an object composed of dense volume. How does it compare to rendering the object directly with a Lambertian material? Why are they alike or different?

    They do share similarities when the volume is very dense. However, I think it is easier to get a correct result of diffuse transparent surface and subsurface scattering does not need special treatment.

3. Jim Kajiya famously has predicted in 1991 that in 10 years, all rendering will be volume rendering. What do you think that makes him think so? Why hasn’t it happened yet?

    I think this is because volumetric rendering is a more general model that includes both transparent and non-transparent objects. The burden is mostly computational cost as real-time high quaility path tracing is still a difficult thing to do.

## 6
1. Play with the index of refraction parameter of the dielectric interface in scenes/volpath_test/volpath_test5_2.xml. How does that affect appearance? Why?

    The higher the IOR, the surface of the medium appears more transparent because less light is reflected and more light is refracted.

2. In the scene scenes/volpath_test/vol_cbox_teapot.xml, we model the glass teapot as a transparent glass with blue homogeneous medium inside. What is the difference in terms of appearance between this approach and just making the color of the glass blue without any medium inside?

    Transparent glass only reflect and refract light, and the color is only determined by light and view direction. A colored media exhibits color by scatter and absorb certain wavelength. So, the actual color also depends on the volume itself. Such as a thicker volume is more diffuse and has a more saturated color.

## 7
1.  For heterogeneous volumes, what kind of distribution of the volume density makes the null scattering efficient/inefficient? Can you think of a way to improve our current sampling scheme in the inefficient case?

    Since the majorant is the upper bound of $\sigma_t$, if $\sigma_t$ is mostly small but is big only at one place, the majorant will also be very big for the entire medium, causing the sampling of $t$ to be very inefficient. One solution is to subdivide the medium so that the $\sigma$ does not change too much within the region.

2. How do we make the null-scattering work for emissive volumes? Briefly describe a solution.

    Similar to transmittance, we first consider the emission from a homogeneous media. Then, the energy emitted by the the null particle needs to be removed.

3. Why is it important to have an unbiased solution for volume rendering? Would it be sensible to have something that is biased but faster? How would you do it?

    I think it is important to have an unbiased solution when the results of rendering needs to be connected to real world, such as simulation or differentiable rendering. A faster but biased solution is needed because human eyes cannot distingush the results when the visual fidelity reaches certain level.