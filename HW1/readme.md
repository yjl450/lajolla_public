Yijian Liu
yil212

# 1

1.1 Both BRDFs are mostly similar. However, the gievn BRDF provides stronger grazing retroreflection response. which makes the edge on the darker side appears lighter.

1.2 The higher the roughness, the surface appears to be lighter, which is more evident on the darker side.

1.3  Base diffuse BRDF is very similar to the lambertion, creating a nice gradient on the surface. The subsurface appears mostly flat and only shows higher response in the grazing angle area. If this area is not in shadow, subsurface is more visible.

1.4 Oren-Nayar BRDF is a microfact BRDF, which is more physically based and thus more accurate. Disney BRDF modifies the Schlick Fresnel, providing a good appearance in general and offering the artists more control. I prefer the Disney BRDF as it looks decent and is easy to manipulate.

# 2
2.1 As the name suggests, DisneyMetal is more like smooth metal and features strong mirror reflection. Roughplastic only combines diffuse and specular component while metal uses a more accurate microfacet model GGX. Each microfacet is treated as a mirror.

2.2 As the roughness gets higher, the surface gets dimmer. This is because as the surface become smoother, more lights are scarterred in a way thay does not contribute to the camera. 

2.3 Beckmann distribution has a shorter tail than Trowbridge-Reitz distribution. Trowbridge-Reitz distribution is used at Disney to more closely match the long tail specular in the MERL dataset.

2.4 Schlick approximation is both easy to compute and implement. However, the actual Fresnel is more accurate and does not have the weird effect of much higher specular in certain conditions.

# 3
3.1 Under similar roughness, the specular of DisneyMetal appears to be more directiona. Clearcoat is more plasticky and its specular is achromatic. This may be because Clearcoat is modeled as a transparent layer and it is designed to have longer tail than the GGX distribution.

3.2 Autodesk Standard Surface's coat layer does not use a long-tail NDF. This can be useful when adding achromatic specular as it will not affect general appearance as much. It is also theoretically more physcially based. However, the missing long tail need to be compensated to fit the MERL dataset and the overall composition is more complex.

3.3 Clearcoat is designed to represent a thin transparent layer that won't typically be anisotropic in real life. Also, it can be difficult to model a transparent while anisotropic material.

# 4
4.1 The square root models two refractions. When light travels through the bounday twice, the refraction results will be combined and get the correct value.

4.2 When the index of refraction is higher, the lights that pass through the material get distorted futher. Thus, material with high $\eta$ looks less transparent and the reflection is much more visible.

4.3 If the object is not closed object, which makes the object appears to have no thickness, the BRDF will not work. This is because there is no space for the light to refract. To make it work, use the thin layer BRDF.

# 5
5.1 It will be completely dark. When the light is at an tilted angle, there is a reflection, representing the grazing reflection.

5.2 The higher the sheentint, the more tinted the reflection tinted by the bass color. When calculating $C_{sheen}$, when shentint is higher, the result is more dominated by $C_{tint}$, which adds more of the base color.

5.3 A microfact BRDF is more accurate but is also more complicated. However, to accurately model the appearance of fabrics, things such as yarns and threads need to be considered on top of microfacet. Disney sheen BRDF provides decent results, which is the one I prefer. 

# 6
6.1 Specular controls the intensity of the specular highlight, and metallic controls how much the surface reflects the light in a mirror-like way. The higher the specular, the highlight becomes brighter. The higher the metallic, the clearer the reflection.

6.2 Roughness controls the general appearance of the surface, mostly the diffusion. On top of that, clearcoat_gloss controls only the concentration of the specular without affecting the diffusion.

6.3 Similar to SheenTint, The higher the specularTint, the more saturated the color of the specular. For metals like gold, the specular is strongly affected by the base color.

6.4 I can see that roughness is a very important surface parameter that has different impact to different materials. If I were an artist, it is better if the parameter has a clearer meaning when it comes to the macro surfaces.

