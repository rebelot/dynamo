// use std::iter::Sum;
// use std::ops::Mul;

use crate::{Rvec, DIM};

#[inline]
pub fn dot(avec: &Rvec, bvec: &Rvec) -> f32 {
    return avec[0] * bvec[0] + avec[1] * bvec[1] + avec[2] * bvec[2];
    // return avec.iter().zip(bvec.iter()).map(|(a, b)| a * b).sum();
}

#[inline]
pub fn cross(avec: &Rvec, bvec: &Rvec) -> Rvec {
    return [
        avec[1] * bvec[2] - avec[2] * bvec[1],
        avec[2] * bvec[0] - avec[0] * bvec[2],
        avec[0] * bvec[1] - avec[1] * bvec[0],
    ];
}

// pub fn angle(avec: &Rvec, bvec: &Rvec) -> f32 {
//     let x = dot(avec, bvec);
//     let y = norm2(&cross(avec, bvec)).sqrt();
//     return y.atan2(x);
// }

#[inline]
pub fn norm2(vec: &[f32; DIM]) -> f32 {
    return dot(vec, vec);
}

#[inline]
pub fn min_image<'a>(vec: &'a mut Rvec, pbc: &Rvec) -> &'a Rvec {
    for i in 0..DIM {
        let half_box = 0.5 * pbc[i];
        if vec[i] > half_box {
            vec[i] -= pbc[i];
        } else if vec[i] < -half_box {
            vec[i] += pbc[i];
        }
    }
    return vec;
}

#[inline]
pub fn wrap<'a>(vec: &'a mut Rvec, pbc: &Rvec) -> &'a Rvec {
    for i in 0..DIM {
        vec[i] %= pbc[i];
    }
    return vec;
}

#[inline]
pub fn displace_vec(avec: &Rvec, bvec: &Rvec) -> Rvec {
    let mut cvec = [0.0; DIM];
    for i in 0..DIM {
        cvec[i] = bvec[i] - avec[i];
    }
    return cvec;
}

// Iterator that return the indices (i, j) corresponding to elements
// of a lower-triangular matrix with shape n x n.
pub fn tril_indices_from(n: usize) -> impl Iterator<Item = (usize, usize)> {
    (0..n).flat_map(move |i| (i..n).map(move |j| (i, j)))
}
