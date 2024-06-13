use std::f32::consts::PI;

use rustfft::num_traits::float::FloatCore;

use crate::linalg::*;

use crate::{Rvec, DIM};

#[inline]
fn harmonic(k: f32, x0: f32, x: f32) -> (f32, f32) {
    let dx = x - x0;
    let kdx = k * dx;
    (0.5 * kdx * dx, -kdx)
}

#[inline]
fn periodic(k: f32, n: f32, x0: f32, x: f32) -> (f32, f32) {
    let dx = x * n - x0;
    (k * (1.0 + dx.cos()), k * n * dx.sin())
}

pub fn bond_harm(k: f32, r0: f32, ri: &Rvec, rj: &Rvec) -> (f32, [Rvec; 2]) {
    let rij = displace_vec(ri, rj);
    let norm2_rij = norm2(&rij);
    let norm_rij = norm2_rij.sqrt();

    let (u, mut f) = harmonic(k, r0, norm_rij);
    f *= norm_rij.recip();
    // f *= norm2_rij.sqrt().recip(); // will this be faster?

    let mut fvec = [[0.0; DIM]; 2];
    for i in 0..DIM {
        let dri = f * -rij[i];
        fvec[0][i] = dri;
        fvec[1][i] = -dri;
    }
    (u, fvec)
}

pub fn angle_harm(k: f32, t0: f32, ri: &Rvec, rj: &Rvec, rk: &Rvec) -> (f32, [Rvec; 3]) {
    let (t, mut fvec) = dthetadr(ri, rj, rk);
    let (u, f) = harmonic(k, t0, t);
    for i in 0..DIM {
        fvec[0][i] *= f;
        fvec[1][i] *= f;
        fvec[2][i] *= f;
    }
    (u, fvec)
}

#[inline]
pub fn dthetadr(ri: &Rvec, rj: &Rvec, rk: &Rvec) -> (f32, [Rvec; 3]) {
    let rji = displace_vec(rj, ri);
    let rjk = displace_vec(rj, rk);

    let norm2_rji = norm2(&rji);
    let norm2_rjk = norm2(&rjk);

    let norm_rji = norm2_rji.sqrt();
    let norm_rjk = norm2_rjk.sqrt();

    let cos_t = dot(&rji, &rjk) * (norm_rji * norm_rjk).recip();
    let t = cos_t.acos();

    let f = (1.0 - cos_t.powi(2)).sqrt().recip();

    let mut fvec = [[0.0; DIM]; 3];
    for i in 0..DIM {
        let uji = rji[i] / norm_rji;
        let ujk = rjk[i] / norm_rjk;

        let dri = f * norm_rji.recip() * (uji * cos_t - ujk);
        let drk = f * norm_rjk.recip() * (ujk * cos_t - uji);

        let drj = -dri - drk;
        fvec[0][i] = dri;
        fvec[1][i] = drj;
        fvec[2][i] = drk;
    }
    (t, fvec)
}

#[inline]
pub fn dphidr(ri: &Rvec, rj: &Rvec, rk: &Rvec, rl: &Rvec) -> (f32, [Rvec; 4]) {
    let rij = displace_vec(ri, rj);
    let rjk = displace_vec(rj, rk);
    let rkl = displace_vec(rk, rl);

    let nijk = cross(&rij, &rjk);
    let njkl = cross(&rjk, &rkl);

    let s = dot(&nijk, &rkl).signum();

    let norm2_nijk = norm2(&nijk);
    let norm2_njkl = norm2(&njkl);
    let norm2_rjk = norm2(&rjk);

    let cos_phi = dot(&nijk, &njkl) * (norm2_nijk * norm2_njkl).sqrt().recip();
    let phi = s * cos_phi.acos();

    let mut fvec = [[0.0; DIM]; 4];

    for i in 0..DIM {
        let dri = rjk[i] / norm2_nijk * nijk[i];
        let drl = -rjk[i] / norm2_njkl * njkl[i];
        let a = rij[i] * rjk[i] / norm2_rjk;
        let b = rkl[i] * rjk[i] / norm2_rjk;
        fvec[0][i] = dri;
        fvec[1][i] = (a - 1f32) * dri - b * drl;
        fvec[2][i] = (b - 1f32) * drl - a * dri;
        fvec[3][i] = drl;
    }
    (phi, fvec)
}

pub fn pdih(
    k: f32,
    n: f32,
    p0: f32,
    ri: &Rvec,
    rj: &Rvec,
    rk: &Rvec,
    rl: &Rvec,
) -> (f32, [Rvec; 4]) {
    let (phi, mut fvec) = dphidr(ri, rj, rk, rl);
    let (u, f) = periodic(k, n, p0, phi);
    // fvec.flatten_mut().iter_mut().for_each(|x| *x *= f);
    for i in 0..DIM {
        fvec[0][i] *= f;
        fvec[1][i] *= f;
        fvec[2][i] *= f;
        fvec[3][i] *= f;
    }
    (u, fvec)
}

pub fn idih_harm(k: f32, p0: f32, ri: &Rvec, rj: &Rvec, rk: &Rvec, rl: &Rvec) -> (f32, [Rvec; 4]) {
    let (phi, mut fvec) = dphidr(ri, rj, rk, rl);
    let (u, f) = harmonic(k, p0, phi);
    for i in 0..DIM {
        fvec[0][i] *= f;
        fvec[1][i] *= f;
        fvec[2][i] *= f;
        fvec[3][i] *= f;
    }
    (u, fvec)
}

pub fn rbdih(
    c0: f32,
    c1: f32,
    c2: f32,
    c3: f32,
    c4: f32,
    c5: f32,
    ri: &Rvec,
    rj: &Rvec,
    rk: &Rvec,
    rl: &Rvec,
) -> (f32, [Rvec; 4]) {
    let (phi, mut fvec) = dphidr(ri, rj, rk, rl);
    let cos_psi = (phi - PI).cos();
    let cos_psi_2 = cos_psi.powi(2);
    let cos_psi_3 = cos_psi_2 * cos_psi;
    let cos_psi_4 = cos_psi_2.powi(2);
    let cos_psi_5 = cos_psi_4 * cos_psi;
    let u = c0 + c1 * cos_psi + c2 * cos_psi_2 + c3 * cos_psi_3 + c4 * cos_psi_4 + c5 * cos_psi_5;
    let f = phi.sin()
        * (c1
            + 2.0 * c2 * cos_psi
            + 3.0 * c3 * cos_psi_2
            + 4.0 * c4 * cos_psi_3
            + 5.0 * c5 * cos_psi_4);
    for i in 0..DIM {
        fvec[0][i] *= f;
        fvec[1][i] *= f;
        fvec[2][i] *= f;
        fvec[3][i] *= f;
    }
    (u, fvec)
}

pub fn lj(c12: f32, c6: f32, ri: &Rvec, rj: &Rvec) -> (f32, [Rvec; 2]) {
    let rij = displace_vec(ri, rj);
    let norm2_rij = norm2(&rij);
    let norm_rij = norm2_rij.sqrt();
    let norm6_rij = norm2_rij.powi(3);
    let norm12_rij = norm6_rij.powi(2);
    let u = c12 * norm12_rij.recip() - c6 * norm6_rij.recip();
    let f = (-12.0 * c12 / (norm_rij * norm12_rij) + 6.0 * c6 / (norm_rij * norm6_rij)) / norm_rij;
    let mut fvec = [[0.0; DIM]; 2];
    for (i, r) in rij.iter().enumerate() {
        let dri = f * -r;
        fvec[0][i] = dri;
        fvec[1][i] = -dri;
    }
    (u, fvec)
}

pub fn buckingham(a: f32, b: f32, c: f32, ri: &Rvec, rj: &Rvec) -> (f32, [Rvec; 2]) {
    let rij = displace_vec(ri, rj);
    let norm2_rij = norm2(&rij);
    let norm_rij = norm2_rij.sqrt();
    let norm6_rij = norm2_rij.powi(3);
    let exp = (-b * norm2_rij).exp();
    let u = a * exp - c / norm6_rij;
    let f = (a * b * exp - 6.0 * c / (norm6_rij * norm_rij)) / norm_rij;
    let mut fvec = [[0.0; DIM]; 2];
    for (i, r) in rij.iter().enumerate() {
        let dri = f * -r;
        fvec[0][i] = dri;
        fvec[1][i] = -dri;
    }
    (u, fvec)
}

pub fn coulomb(qi: f32, qj: f32, ri: &Rvec, rj: &Rvec) -> (f32, [Rvec; 2]) {
    let rij = displace_vec(ri, rj);
    let norm2_rij = norm2(&rij);
    let norm_rij = norm2_rij.sqrt();
    let qiqj = 138.935_49 * qi * qj;
    let u = qiqj / norm_rij;
    let f = -u / norm2_rij / norm_rij;
    let mut fvec = [[0.0; DIM]; 2];
    for (i, r) in rij.iter().enumerate() {
        let dri = f * -r;
        fvec[0][i] = dri;
        fvec[1][i] = -dri;
    }
    (u, fvec)
}

pub fn comb_rule_geom(vi: f32, wi: f32, vj: f32, wj: f32) -> (f32, f32) {
    let v = (vi * vj).sqrt();
    let w = (wi * wj).sqrt();
    (v, w)
}

pub fn comb_rule_LB(vi: f32, wi: f32, vj: f32, wj: f32) -> (f32, f32) {
    let v = 0.5 * (vi + vj);
    let w = (wi * wj).sqrt();
    (v, w)
}
