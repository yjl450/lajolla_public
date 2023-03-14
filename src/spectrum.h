#pragma once

#include "vector.h"
#include <vector>

/// For now, lajolla assumes we are operating in the linear and trimulus RGB color space.
/// In the future we might implement a proper spectral renderer.
using Spectrum = Vector3;

inline Spectrum make_zero_spectrum() {
    return Vector3{0, 0, 0};
}

inline Spectrum make_const_spectrum(Real v) {
    return Vector3{v, v, v};
}

inline Spectrum fromRGB(const Vector3 &rgb) {
    return rgb;
}

inline Spectrum sqrt(const Spectrum &s) {
    return Vector3{sqrt(max(s[0], Real(0))),
                   sqrt(max(s[1], Real(0))),
                   sqrt(max(s[2], Real(0)))};
}

inline Spectrum exp(const Spectrum &s) {
    return Vector3{exp(s[0]), exp(s[1]), exp(s[2])};
}

inline Real luminance(const Spectrum &s) {
    return s.x * Real(0.212671) + s.y * Real(0.715160) + s.z * Real(0.072169);
}

inline Real avg(const Spectrum &s) {
    return (s.x + s.y + s.z) / 3;
}

inline Vector3 toRGB(const Spectrum &s) {
    return s;
}

/// To support spectral data, we need to convert spectral measurements (how much energy at each wavelength) to
/// RGB. To do this, we first convert the spectral data to CIE XYZ, by
/// integrating over the XYZ response curve. Here we use an analytical response
/// curve proposed by Wyman et al.: https://jcgt.org/published/0002/02/01/
inline Real xFit_1931(Real wavelength) {
    Real t1 = (wavelength - Real(442.0)) * ((wavelength < Real(442.0)) ? Real(0.0624) : Real(0.0374));
    Real t2 = (wavelength - Real(599.8)) * ((wavelength < Real(599.8)) ? Real(0.0264) : Real(0.0323));
    Real t3 = (wavelength - Real(501.1)) * ((wavelength < Real(501.1)) ? Real(0.0490) : Real(0.0382));
    return Real(0.362) * exp(-Real(0.5) * t1 * t1) + 
           Real(1.056) * exp(-Real(0.5) * t2 * t2) -
           Real(0.065) * exp(-Real(0.5) * t3 * t3);
}
inline Real yFit_1931(Real wavelength) {
    Real t1 = (wavelength - Real(568.8)) * ((wavelength < Real(568.8)) ? Real(0.0213) : Real(0.0247));
    Real t2 = (wavelength - Real(530.9)) * ((wavelength < Real(530.9)) ? Real(0.0613) : Real(0.0322));
    return Real(0.821) * exp(-Real(0.5) * t1 * t1) +
           Real(0.286) * exp(-Real(0.5) * t2 * t2);
}
inline Real zFit_1931(Real wavelength) {
    Real t1 = (wavelength - Real(437.0)) * ((wavelength < Real(437.0)) ? Real(0.0845) : Real(0.0278));
    Real t2 = (wavelength - Real(459.0)) * ((wavelength < Real(459.0)) ? Real(0.0385) : Real(0.0725));
    return Real(1.217) * exp(-Real(0.5) * t1 * t1) +
           Real(0.681) * exp(-Real(0.5) * t2 * t2);
}
inline Vector3 XYZintegral_coeff(Real wavelength) {
    return Vector3{xFit_1931(wavelength), yFit_1931(wavelength), zFit_1931(wavelength)};
}

inline Vector3 integrate_XYZ(const std::vector<std::pair<Real, Real>> &data) {
    static const Real CIE_Y_integral = 106.856895;
    static const Real wavelength_beg = 400;
    static const Real wavelength_end = 700;
    if (data.size() == 0) {
        return Vector3{0, 0, 0};
    }
    Vector3 ret = Vector3{0, 0, 0};
    int data_pos = 0;
    // integrate from wavelength 400 nm to 700 nm, increment by 1nm at a time
    // linearly interpolate from the data
    for (Real wavelength = wavelength_beg; wavelength <= wavelength_end; wavelength += Real(1)) {
        // assume the spectrum data is sorted by wavelength
        // move data_pos such that wavelength is between two data or at one end
        while(data_pos < (int)data.size() - 1 &&
               !((data[data_pos].first <= wavelength &&
                  data[data_pos + 1].first > wavelength) ||
                 data[0].first > wavelength)) {
            data_pos += 1;
        }
        Real measurement = 0;
        if (data_pos < (int)data.size() - 1 && data[0].first <= wavelength) {
            Real curr_data = data[data_pos].second;
            Real next_data = data[std::min(data_pos + 1, (int)data.size() - 1)].second;
            Real curr_wave = data[data_pos].first;
            Real next_wave = data[std::min(data_pos + 1, (int)data.size() - 1)].first;
            // linearly interpolate
            measurement = curr_data * (next_wave - wavelength) / (next_wave - curr_wave) +
                          next_data * (wavelength - curr_wave) / (next_wave - curr_wave);
        } else {
            // assign the endpoint
            measurement = data[data_pos].second;
        }
        Vector3 coeff = XYZintegral_coeff(wavelength);
        ret += coeff * measurement;
    }
    Real wavelength_span = wavelength_end - wavelength_beg;
    ret *= (wavelength_span / (CIE_Y_integral * (wavelength_end - wavelength_beg)));
    return ret;
}

inline Vector3 XYZ_to_RGB(const Vector3 &xyz) {
    return Vector3{
        Real( 3.240479) * xyz[0] - Real(1.537150) * xyz[1] - Real(0.498535) * xyz[2],
        Real(-0.969256) * xyz[0] + Real(1.875991) * xyz[1] + Real(0.041556) * xyz[2],
        Real( 0.055648) * xyz[0] - Real(0.204043) * xyz[1] + Real(1.057311) * xyz[2]};
}

inline Vector3 sRGB_to_RGB(const Vector3 &srgb) {
    // https://en.wikipedia.org/wiki/SRGB#From_sRGB_to_CIE_XYZ
    Vector3 rgb = srgb;
    for (int i = 0; i < 3; i++) {
        rgb[i] = rgb[i] <= Real(0.04045) ?
            rgb[i] / Real(12.92) :
            pow((rgb[i] + Real(0.055)) / Real(1.055), Real(2.4));
    }
    return rgb;
}
///////////////////////////////////////////////////////////////////////
// Spectral rendering from pbrt-v3
///////////////////////////////////////////////////////////////////////
inline Real Lerp(Real t, Real v1, Real v2) { return (1 - t) * v1 + t * v2; }

inline Real AverageSpectrumSamples(const Real* lambda, const Real* vals, int n,
    Real lambdaStart, Real lambdaEnd) {
    // Handle cases with out-of-bounds range or single sample only
    if (lambdaEnd <= lambda[0]) return vals[0];
    if (lambdaStart >= lambda[n - 1]) return vals[n - 1];
    if (n == 1) return vals[0];
    Real sum = 0;
    // Add contributions of constant segments before/after samples
    if (lambdaStart < lambda[0]) sum += vals[0] * (lambda[0] - lambdaStart);
    if (lambdaEnd > lambda[n - 1])
        sum += vals[n - 1] * (lambdaEnd - lambda[n - 1]);

    // Advance to first relevant wavelength segment
    int i = 0;
    while (lambdaStart > lambda[i + 1]) ++i;

    // Loop over wavelength sample segments and add contributions
    auto interp = [lambda, vals](Real w, int i) {
        return Lerp((w - lambda[i]) / (lambda[i + 1] - lambda[i]), vals[i],
            vals[i + 1]);
    };
    for (; i + 1 < n && lambdaEnd >= lambda[i]; ++i) {
        Real segLambdaStart = std::max(lambdaStart, lambda[i]);
        Real segLambdaEnd = std::min(lambdaEnd, lambda[i + 1]);
        sum += 0.5 * (interp(segLambdaStart, i) + interp(segLambdaEnd, i)) *
            (segLambdaEnd - segLambdaStart);
    }
    return sum / (lambdaEnd - lambdaStart);
}

static const int nCIESamples = 471;
static const Real CIE_Y_integral = 106.856895;
static const int nRGB2SpectSamples = 32;
static const int nSpectralSamples = 60;
extern const Real CIE_X[nCIESamples];
extern const Real CIE_Y[nCIESamples];
extern const Real CIE_Z[nCIESamples];
extern const Real CIE_lambda[nCIESamples];
extern const Real RGB2SpectLambda[nRGB2SpectSamples];
extern const Real RGBRefl2SpectWhite[nRGB2SpectSamples];
extern const Real RGBRefl2SpectCyan[nRGB2SpectSamples];
extern const Real RGBRefl2SpectMagenta[nRGB2SpectSamples];
extern const Real RGBRefl2SpectYellow[nRGB2SpectSamples];
extern const Real RGBRefl2SpectRed[nRGB2SpectSamples];
extern const Real RGBRefl2SpectGreen[nRGB2SpectSamples];
extern const Real RGBRefl2SpectBlue[nRGB2SpectSamples];
extern const Real RGBIllum2SpectWhite[nRGB2SpectSamples];
extern const Real RGBIllum2SpectCyan[nRGB2SpectSamples];
extern const Real RGBIllum2SpectMagenta[nRGB2SpectSamples];
extern const Real RGBIllum2SpectYellow[nRGB2SpectSamples];
extern const Real RGBIllum2SpectRed[nRGB2SpectSamples];
extern const Real RGBIllum2SpectGreen[nRGB2SpectSamples];
extern const Real RGBIllum2SpectBlue[nRGB2SpectSamples];

class WaveSpectrum {
private:
    Real c[nSpectralSamples];
    static WaveSpectrum X, Y, Z;
    static WaveSpectrum rgbRefl2SpectWhite, rgbRefl2SpectCyan;
    static WaveSpectrum rgbRefl2SpectMagenta, rgbRefl2SpectYellow;
    static WaveSpectrum rgbRefl2SpectRed, rgbRefl2SpectGreen;
    static WaveSpectrum rgbRefl2SpectBlue;
    static WaveSpectrum rgbIllum2SpectWhite, rgbIllum2SpectCyan;
    static WaveSpectrum rgbIllum2SpectMagenta, rgbIllum2SpectYellow;
    static WaveSpectrum rgbIllum2SpectRed, rgbIllum2SpectGreen;
    static WaveSpectrum rgbIllum2SpectBlue;
    static WaveSpectrum centeral_wavelength;
public:
    static const int wavelength_beg = 400;
    static const int wavelength_end = 700;
    static Real interval;

    WaveSpectrum() {
        memset(c, Real(0), nSpectralSamples * sizeof(Real));
    }

    WaveSpectrum(Real n) {
        for (int i = 0; i < nSpectralSamples; i++)
        {
            c[i] = n;
        }
    }

    static inline WaveSpectrum fromRGB(Spectrum color, bool is_reflectance = true);

    static void Init() {
        interval = (Real(wavelength_end) - wavelength_beg) / nSpectralSamples;
        // Compute XYZ matching functions for _SampledSpectrum_
        for (int i = 0; i < nSpectralSamples; ++i) {
            Real wl0 = Lerp(Real(i) / Real(nSpectralSamples),
                wavelength_beg, wavelength_end);
            Real wl1 = Lerp(Real(i + 1) / Real(nSpectralSamples),
                wavelength_beg, wavelength_end);
            X.c[i] = AverageSpectrumSamples(CIE_lambda, CIE_X, nCIESamples, wl0,
                wl1);
            Y.c[i] = AverageSpectrumSamples(CIE_lambda, CIE_Y, nCIESamples, wl0,
                wl1);
            Z.c[i] = AverageSpectrumSamples(CIE_lambda, CIE_Z, nCIESamples, wl0,
                wl1);
        }

        // Compute RGB to spectrum functions for _SampledSpectrum_
        for (int i = 0; i < nSpectralSamples; ++i) {
            Real wl0 = Lerp(Real(i) / Real(nSpectralSamples),
                wavelength_beg, wavelength_end);
            Real wl1 = Lerp(Real(i + 1) / Real(nSpectralSamples),
                wavelength_beg, wavelength_end);
            rgbRefl2SpectWhite.c[i] =
                AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectWhite,
                    nRGB2SpectSamples, wl0, wl1);
            rgbRefl2SpectCyan.c[i] =
                AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectCyan,
                    nRGB2SpectSamples, wl0, wl1);
            rgbRefl2SpectMagenta.c[i] =
                AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectMagenta,
                    nRGB2SpectSamples, wl0, wl1);
            rgbRefl2SpectYellow.c[i] =
                AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectYellow,
                    nRGB2SpectSamples, wl0, wl1);
            rgbRefl2SpectRed.c[i] = AverageSpectrumSamples(
                RGB2SpectLambda, RGBRefl2SpectRed, nRGB2SpectSamples, wl0, wl1);
            rgbRefl2SpectGreen.c[i] =
                AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectGreen,
                    nRGB2SpectSamples, wl0, wl1);
            rgbRefl2SpectBlue.c[i] =
                AverageSpectrumSamples(RGB2SpectLambda, RGBRefl2SpectBlue,
                    nRGB2SpectSamples, wl0, wl1);

            rgbIllum2SpectWhite.c[i] =
                AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectWhite,
                    nRGB2SpectSamples, wl0, wl1);
            rgbIllum2SpectCyan.c[i] =
                AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectCyan,
                    nRGB2SpectSamples, wl0, wl1);
            rgbIllum2SpectMagenta.c[i] =
                AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectMagenta,
                    nRGB2SpectSamples, wl0, wl1);
            rgbIllum2SpectYellow.c[i] =
                AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectYellow,
                    nRGB2SpectSamples, wl0, wl1);
            rgbIllum2SpectRed.c[i] =
                AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectRed,
                    nRGB2SpectSamples, wl0, wl1);
            rgbIllum2SpectGreen.c[i] =
                AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectGreen,
                    nRGB2SpectSamples, wl0, wl1);
            rgbIllum2SpectBlue.c[i] =
                AverageSpectrumSamples(RGB2SpectLambda, RGBIllum2SpectBlue,
                    nRGB2SpectSamples, wl0, wl1);

            centeral_wavelength[i] = 0.5 * (wavelength_beg + i * interval) + 0.5 * (wavelength_beg + (i + 1) * interval);
        }
    }
    Real operator [](int i) const { return c[i]; }

    Real& operator [](int i) { return c[i]; }

    WaveSpectrum operator+ (Real n) {
        WaveSpectrum wave = WaveSpectrum();
        for (int i = 0; i < nSpectralSamples; i++)
        {
            wave[i] = c[i] + n;
        }
        return wave;
    }

    WaveSpectrum operator+ (WaveSpectrum& n) {
        WaveSpectrum wave = WaveSpectrum();
        for (int i = 0; i < nSpectralSamples; i++)
        {
            wave[i] = c[i] + n[i];
        }
        return wave;
    }

    WaveSpectrum operator* (Real n) {
        WaveSpectrum wave = WaveSpectrum();
        for (int i = 0; i < nSpectralSamples; i++)
        {
            wave[i] = c[i] * n;
        }
        return wave;
    }

    WaveSpectrum operator* (WaveSpectrum& n) {
        WaveSpectrum wave = WaveSpectrum();
        for (int i = 0; i < nSpectralSamples; i++)
        {
            wave[i] = c[i] * n[i];
        }
        return wave;
    }

    WaveSpectrum operator- (Real n) {
        WaveSpectrum wave = WaveSpectrum();
        for (int i = 0; i < nSpectralSamples; i++)
        {
            wave[i] = c[i] - n;
        }
        return wave;
    }

    WaveSpectrum operator- (WaveSpectrum& n) {
        WaveSpectrum wave = WaveSpectrum();
        for (int i = 0; i < nSpectralSamples; i++)
        {
            wave[i] = c[i] - n[i];
        }
        return wave;
    }

    WaveSpectrum operator/ (Real n) {
        WaveSpectrum wave = WaveSpectrum();
        for (int i = 0; i < nSpectralSamples; i++)
        {
            wave[i] = c[i] / n;
        }
        return wave;
    }

    WaveSpectrum operator/ (WaveSpectrum& n) {
        WaveSpectrum wave = WaveSpectrum();
        for (int i = 0; i < nSpectralSamples; i++)
        {
            wave[i] = c[i] / n[i];
        }
        return wave;
    }

    void operator+= (Real n) {
        for (int i = 0; i < nSpectralSamples; i++)
        {
            c[i] += n;
        }
    }

    void operator+= (WaveSpectrum& n) {
        for (int i = 0; i < nSpectralSamples; i++)
        {
            c[i] += n[i];
        }
    }

    void operator-= (Real n) {
        for (int i = 0; i < nSpectralSamples; i++)
        {
            c[i] -= n;
        }
    }

    void operator-= (WaveSpectrum& n) {
        for (int i = 0; i < nSpectralSamples; i++)
        {
            c[i] -= n[i];
        }
    }

    void operator*= (Real n) {
        for (int i = 0; i < nSpectralSamples; i++)
        {
            c[i] *= n;
        }
    }

    void operator*= (WaveSpectrum& n) {
        for (int i = 0; i < nSpectralSamples; i++)
        {
            c[i] *= n[i];
        }
    }

    void operator/= (Real n) {
        for (int i = 0; i < nSpectralSamples; i++)
        {
            c[i] /= n;
        }
    }

    void operator/= (WaveSpectrum& n) {
        for (int i = 0; i < nSpectralSamples; i++)
        {
            c[i] /= n[i];
        }
    }

    void Clamp(Real low = 0, Real high = infinity<Real>()) {
        for (int i = 0; i < nSpectralSamples; i++)
        {
            c[i] = max(min(c[i], high), low);
        }
    }

    void ToXYZ(Real xyz[3]) const {
        xyz[0] = xyz[1] = xyz[2] = Real(0);
        for (int i = 0; i < nSpectralSamples; ++i) {
            xyz[0] += X.c[i] * c[i];
            xyz[1] += Y.c[i] * c[i];
            xyz[2] += Z.c[i] * c[i];
        }
        Real scale = Real(wavelength_end - wavelength_beg) /
            Real(CIE_Y_integral * nSpectralSamples);
        xyz[0] *= scale;
        xyz[1] *= scale;
        xyz[2] *= scale;
    }

    inline Spectrum XYZToRGB(const Real xyz[3]) {
        Spectrum rgb;
        rgb[0] = 3.240479f * xyz[0] - 1.537150f * xyz[1] - 0.498535f * xyz[2];
        rgb[1] = -0.969256f * xyz[0] + 1.875991f * xyz[1] + 0.041556f * xyz[2];
        rgb[2] = 0.055648f * xyz[0] - 0.204043f * xyz[1] + 1.057311f * xyz[2];
        return rgb;
    }

    Spectrum to_RGB() {
        Real xyz[3];
        ToXYZ(xyz);
        return XYZToRGB(xyz);
    }
};

WaveSpectrum WaveSpectrum::fromRGB(Spectrum rgb, bool is_reflectance) {
    WaveSpectrum r;
    if (is_reflectance) {
        // Convert reflectance spectrum to RGB
    if (rgb[0] <= rgb[1] && rgb[0] <= rgb[2]) {
        // Compute reflectance _SampledSpectrum_ with _rgb[0]_ as minimum
        r += rgbRefl2SpectWhite * rgb[0];
        if (rgb[1] <= rgb[2]) {
            r += rgbRefl2SpectCyan * (rgb[1] - rgb[0]);
            r += rgbRefl2SpectBlue * (rgb[2] - rgb[1]);
        }
        else {
            r += rgbRefl2SpectCyan * (rgb[2] - rgb[0]);
            r += rgbRefl2SpectGreen * (rgb[1] - rgb[2]);
        }
    }
    else if (rgb[1] <= rgb[0] && rgb[1] <= rgb[2]) {
        // Compute reflectance _SampledSpectrum_ with _rgb[1]_ as minimum
        r += rgbRefl2SpectWhite * rgb[1];
        if (rgb[0] <= rgb[2]) {
            r += rgbRefl2SpectMagenta * (rgb[0] - rgb[1]);
            r += rgbRefl2SpectBlue * (rgb[2] - rgb[0]);
        }
        else {
            r += rgbRefl2SpectMagenta * (rgb[2] - rgb[1]);
            r += rgbRefl2SpectRed * (rgb[0] - rgb[2]);
        }
    }
    else {
        // Compute reflectance _SampledSpectrum_ with _rgb[2]_ as minimum
        r += rgbRefl2SpectWhite * rgb[2];
        if (rgb[0] <= rgb[1]) {
            r += rgbRefl2SpectYellow * (rgb[0] - rgb[2]);
            r += rgbRefl2SpectGreen * (rgb[1] - rgb[0]);
        }
        else {
            r += rgbRefl2SpectYellow * (rgb[1] - rgb[2]);
            r += rgbRefl2SpectRed * (rgb[0] - rgb[1]);
        }
    }
    r *= .94;
    }
    else {
        // Convert illuminant spectrum to RGB
        if (rgb[0] <= rgb[1] && rgb[0] <= rgb[2]) {
            // Compute illuminant _SampledSpectrum_ with _rgb[0]_ as minimum
            r += rgbIllum2SpectWhite * rgb[0];
            if (rgb[1] <= rgb[2]) {
                r += rgbIllum2SpectCyan * (rgb[1] - rgb[0]);
                r += rgbIllum2SpectBlue * (rgb[2] - rgb[1]);
            }
            else {
                r += rgbIllum2SpectCyan * (rgb[2] - rgb[0]);
                r += rgbIllum2SpectGreen * (rgb[1] - rgb[2]);
            }
        }
        else if (rgb[1] <= rgb[0] && rgb[1] <= rgb[2]) {
            // Compute illuminant _SampledSpectrum_ with _rgb[1]_ as minimum
            r += rgbIllum2SpectWhite * rgb[1];
            if (rgb[0] <= rgb[2]) {
                r += rgbIllum2SpectMagenta * (rgb[0] - rgb[1]);
                r += rgbIllum2SpectBlue * (rgb[2] - rgb[0]);
            }
            else {
                r += rgbIllum2SpectMagenta * (rgb[2] - rgb[1]);
                r += rgbIllum2SpectRed * (rgb[0] - rgb[2]);
            }
        }
        else {
            // Compute illuminant _SampledSpectrum_ with _rgb[2]_ as minimum
            r += rgbIllum2SpectWhite * rgb[2];
            if (rgb[0] <= rgb[1]) {
                r += rgbIllum2SpectYellow * (rgb[0] - rgb[2]);
                r += rgbIllum2SpectGreen * (rgb[1] - rgb[0]);
            }
            else {
                r += rgbIllum2SpectYellow * (rgb[1] - rgb[2]);
                r += rgbIllum2SpectRed * (rgb[0] - rgb[1]);
            }
        }
        r *= .86445f;
    }
    r.Clamp();
    return r;
}

