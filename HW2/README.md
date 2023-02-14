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

    Equation 13 is used to compute the trans_pdf so that when the emission from a surfce is included, the transmittance is devided by the trans_pdf. This is because the case includes all the possibilities of $t\geq t_{hit}$, so the overall probability is the integral on $t\in(t_{hit},\infty)$.

3. Play with the parameters σs and σa, how do they affect the final image? Why? (we assume monochromatic volumes, so don’t set the parameters to be different for each channel.)

    Both parameters will increase the energy loss as they get higher. $\sigma_a$ direct affects the absorption so the image gets darker significantly as $\sigma_a$ increases.

    $\sigma_s$ both contributes the energy loss and increase. So, a higher $\sigma_s$ does not dim the image as much. However, as the scattering increases, the image becomes more blurry as the light bounces more in the medium.

4. Change the phase function from isotropic (the default one) to Henyey-Greenstein by uncommenting
the phase function specification in the scene scenes/volpath_test/volpath_test2.xml. Play with the
parameter g (the valid range is (−1, 1)). What does g mean? How does the g parameter change the
appearance? Why?

    g means the intensity of forward scattering (positive) and back scattering (negative). As the absolute value of g gets bigger, the scattering intensity along the incident axis is stonger and scattering from the sides becomes weaker, which means the purple light becomes less visible. Back scattering also leads to less light entering the camera, making the image darker overall.
