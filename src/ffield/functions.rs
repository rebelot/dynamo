use crate::linalg::*;

use crate::{Rvec, DIM};

#[inline]
fn harmonic(k: f32, x0: f32, x: f32) -> (f32, f32) {
    let dx = x - x0;
    let kdx = k * dx;
    return (0.5 * kdx * dx, -kdx);
}

#[inline]
fn periodic(k: f32, n: f32, x0: f32, x: f32) -> (f32, f32) {
    let dx = x * n - x0;
    return (k * (1.0 + dx.cos()), k * n * dx.sin());
}

pub fn bond_harm(k: f32, r0: f32, ri: &Rvec, rj: &Rvec) -> (f32, [Rvec; 2]) {
    let rij = displace_vec(ri, rj);
    let norm2_rij = norm2(&rij);
    let norm_rij = norm2_rij.sqrt();

    let (u, mut f) = harmonic(k, r0, norm_rij);
    f *= norm_rij.recip();
    // *f *= norm2_rij.sqrt().recip(); // will this be faster?

    let mut fvec = [[0.0; DIM]; 2];
    let mut dri: f32;
    for i in 0..DIM {
        dri = f * -rij[i];
        fvec[0][i] = dri;
        fvec[1][i] = -dri;
    }
    return (u, fvec);
}

pub fn angle_harm(k: f32, t0: f32, ri: &Rvec, rj: &Rvec, rk: &Rvec) -> (f32, [Rvec; 3]) {
    let rji = displace_vec(rj, ri);
    let rjk = displace_vec(rj, rk);

    let norm_rji_rjk_1 = (norm2(&rji) * norm2(&rjk)).sqrt().recip();

    let cos_t = dot(&rji, &rjk) * norm_rji_rjk_1;
    let t = cos_t.acos();

    let (u, mut f) = harmonic(k, t0, t);

    f *= -(1.0 - cos_t * cos_t).sqrt().recip() * norm_rji_rjk_1;

    let mut fvec = [[0.0; DIM]; 3];
    let [mut dri, mut drj, mut drk]: [f32; 3];
    for i in 0..DIM {
        dri = f * rjk[i];
        drk = f * rji[i];
        drj = -dri - drk;
        fvec[0][i] = dri;
        fvec[1][i] = drj;
        fvec[2][i] = drk;
    }
    return (u, fvec);
}

#[inline]
pub fn drdpsi(ri: &Rvec, rj: &Rvec, rk: &Rvec, rl: &Rvec) -> (f32, [Rvec; 4]) {
    let rij = displace_vec(ri, rj);
    let rjk = displace_vec(rj, rk);
    let rkl = displace_vec(rk, rl);

    let nijk = cross(&rij, &rjk);
    let njkl = cross(&rjk, &rkl);

    let y = dot(&cross(&nijk, &njkl), &rjk) * norm2(&rjk).sqrt().recip();
    let x = dot(&nijk, &njkl);
    let psi = y.atan2(x);

    let f = (1.0 + (y / x).powi(2)) * y * x.sqrt().recip();

    let dnijk_drij = [rjk[2] - rjk[1], rjk[0] - rjk[2], rjk[1] - rjk[0]];
    let dnijk_drjk = [rij[1] - rij[2], rij[2] - rij[0], rij[0] - rij[1]];
    let dnjkl_drjk = [rkl[2] - rkl[1], rkl[0] - rkl[2], rkl[1] - rkl[0]];
    let dnjkl_drkl = [rjk[1] - rjk[2], rjk[2] - rjk[0], rjk[0] - rjk[1]];

    let mut fvec = [[0.0; DIM]; 4];
    let [mut drij, mut drjk, mut drkl]: [f32; 3];
    for i in 0..DIM {
        drij = f * njkl[i] * dnijk_drij[i];
        drjk = f * (njkl[i] * dnijk_drjk[i] + nijk[i] * dnjkl_drjk[i]);
        drkl = f * nijk[i] * dnjkl_drkl[i];
        fvec[0][i] = -drij;
        fvec[1][i] = drij - drjk;
        fvec[2][i] = drjk - drkl;
        fvec[3][i] = drkl;
    }
    return (psi, fvec);
}

pub fn pdih(k: f32, n: f32, p0: f32, ri: &Rvec, rj: &Rvec, rk: &Rvec, rl: &Rvec) -> (f32, [Rvec; 4]) {
    let (psi, mut fvec) = drdpsi(ri, rj, rk, rl);
    let (u, f) = periodic(k, n, p0, psi);
    for i in 0..DIM {
        fvec[0][i] *= f;
        fvec[1][i] *= f;
        fvec[2][i] *= f;
        fvec[3][i] *= f;
    }
    return (u, fvec);
}

pub fn idih_harm(k: f32, p0: f32, ri: &Rvec, rj: &Rvec, rk: &Rvec, rl: &Rvec) -> (f32, [Rvec; 4]) {
    let (psi, mut fvec) = drdpsi(ri, rj, rk, rl);
    let (u, f) = harmonic(k, p0, psi);
    for i in 0..DIM {
        fvec[0][i] *= f;
        fvec[1][i] *= f;
        fvec[2][i] *= f;
        fvec[3][i] *= f;
    }
    return (u, fvec);
}
